#pragma once

#include "openfhe.h"

// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/bgvrns/bgvrns-ser.h"

#include "emp-tool/emp-tool.h"
#include <math.h>
#include <future>
#include "openfhecore.h"
#include "emp-aby/io/mp_io_channel.h"

// Required to compile on mac, remove on ubuntu
#ifdef __APPLE__
    std::shared_ptr<lbcrypto::PRNG> lbcrypto::PseudoRandomNumberGenerator::m_prng = nullptr;
#endif

namespace emp {

#define MAX_MULT_DEPTH 10

template <typename IO>
class HE {
private:
    ThreadPool* pool;
    lbcrypto::KeyPair<lbcrypto::DCRTPoly> kp;

public:
    PRG prg;
    int party, mult_depth = 3, add_count = 100;
    MPIOChannel<IO>* io;
    int num_party;
    long long int q                                = (1L << 16) + 1;
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc = {};
    std::shared_ptr<std::map<usint, lbcrypto::EvalKey<lbcrypto::DCRTPoly>>> evalAtIndexKeys;
    lbcrypto::PublicKey<lbcrypto::DCRTPoly> pk;

    HE(int num_party, MPIOChannel<IO>* io, ThreadPool* pool, int party, long long int plaintext_mod = (1L << 16) + 1,
       int mult_depth = -1, bool keygen = true, bool YAO = false, bool mult = false, int add_count = 100) {
        this->io         = io;
        this->party      = party;
        this->num_party  = num_party;
        this->pool       = pool;
        this->q          = plaintext_mod;
        this->mult_depth = mult_depth;
        this->add_count = add_count;
        if (!YAO && mult_depth == -1) {
            if (num_party <= MAX_MULT_DEPTH) {
                this->mult_depth = num_party;
            }
            else {
                this->mult_depth              = MAX_MULT_DEPTH;
                int num_bootstrapping_parties = floor((double)(num_party - 1) / (double)this->mult_depth);
                this->mult_depth              = ceil((double)num_party / (double)(num_bootstrapping_parties + 1));
            }
        }
        else {
            if (mult) {
                this->mult_depth = 1;
            }
        }
        if (lbcrypto::Serial::DeserializeFromFile(
                "data/cc_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", cc,
                lbcrypto::SerType::BINARY)) {
            lbcrypto::Serial::DeserializeFromFile(
                "data/pk_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", pk,
                lbcrypto::SerType::BINARY);
            lbcrypto::Serial::DeserializeFromFile(
                "data/kp_public_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", kp.publicKey,
                lbcrypto::SerType::BINARY);
            lbcrypto::Serial::DeserializeFromFile(
                "data/kp_secret_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", kp.secretKey,
                lbcrypto::SerType::BINARY);
            // lbcrypto::Serial::DeserializeFromFile("data/evalAtIndexKeys_"+std::to_string(party)+"_"+std::to_string(num_party)+".txt", evalAtIndexKeys, lbcrypto::SerType::BINARY);
        }
        else {
            setup_cryptocontext(YAO);

            if (keygen)
                multiparty_keygen();
        }
    }

    ~HE() {
        lbcrypto::Serial::SerializeToFile("data/cc_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                                          cc, lbcrypto::SerType::BINARY);
        lbcrypto::Serial::SerializeToFile("data/pk_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                                          pk, lbcrypto::SerType::BINARY);
        // lbcrypto::Serial::SerializeToFile("data/evalAtIndexKeys_"+std::to_string(party)+"_"+std::to_string(num_party)+".txt", evalAtIndexKeys, lbcrypto::SerType::BINARY);
        lbcrypto::Serial::SerializeToFile(
            "data/kp_secret_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", kp.secretKey,
            lbcrypto::SerType::BINARY);
        lbcrypto::Serial::SerializeToFile(
            "data/kp_public_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt", kp.publicKey,
            lbcrypto::SerType::BINARY);
        // lbcrypto::Serial::SerializeToFile("data/evalMultKey"+)
        std::ofstream file("data/eval_mult_key_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                           std::ios::out | std::ios::binary);
        if (file.is_open()) {
            cc->SerializeEvalMultKey(file, lbcrypto::SerType::BINARY, cc);
            file.close();
        }
        std::ofstream file1("data/eval_auto_key_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                            std::ios::out | std::ios::binary);
        if (file1.is_open()) {
            cc->SerializeEvalAutomorphismKey(file1, lbcrypto::SerType::BINARY, cc);
            file1.close();
        }
    }

    template <typename T>
    void serialize_send(T& obj, int i, int j = 0, MESSAGE_TYPE mt = NORM_MSG) {
        std::stringstream s;

        lbcrypto::Serial::Serialize(obj, s, lbcrypto::SerType::BINARY);
        string str      = s.str();
        int string_size = str.size();
        char* c         = (char*)malloc(string_size);
        s.read(c, string_size);
        io->send_data(i, c, string_size, j, mt);
        io->flush(i, j);
        // free(c);
    }

    template <typename T>
    void serialize_sendall(T& obj, int j = 0, MESSAGE_TYPE mt = NORM_MSG) {
        std::stringstream s;

        lbcrypto::Serial::Serialize(obj, s, lbcrypto::SerType::BINARY);
        string str      = s.str();
        int string_size = str.size();
        char* c         = (char*)malloc(string_size);
        s.read(c, string_size);
        std::vector<std::future<void>> res;
        for (int i = 1; i <= num_party; ++i) {
            if (i != party) {
                res.push_back(std::async([this, i, c, string_size, j, mt]() {
                    io->send_data(i, c, string_size, j, mt);
                    io->flush(i, j);
                }));
            }
        }
        for (auto& fut : res)
            fut.get();
        res.clear();
        // free(c);
    }

    template <typename T>
    void deserialize_recv(T& obj, int i, int j = 0, MESSAGE_TYPE mt = NORM_MSG) {
        std::stringstream s;
        int string_size = 0;
        char* c         = (char*)io->recv_data(i, string_size, j, mt);
        s.write(c, string_size);
        // free(c);
        lbcrypto::Serial::Deserialize(obj, s, lbcrypto::SerType::BINARY);
    }

    void setup_cryptocontext(bool YAO = false) {
        if (party == ALICE) {
            lbcrypto::CCParams<lbcrypto::CryptoContextBGVRNS> parameters;
            parameters.SetPlaintextModulus(q);
            parameters.SetSecretKeyDist(UNIFORM_TERNARY);
            parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
            parameters.SetMultipartyMode(FIXED_NOISE_MULTIPARTY);
            if (YAO) {
                parameters.SetKeySwitchCount(0);
                parameters.SetFirstModSize(17);
                parameters.SetMultiplicativeDepth(0);
                parameters.SetScalingTechnique(FLEXIBLEAUTO);
            }
            else {
                parameters.SetMultiplicativeDepth(mult_depth);
                parameters.SetKeySwitchCount(2);
                parameters.SetEvalAddCount(add_count);
                parameters.SetScalingTechnique(FLEXIBLEAUTOEXT);
                // parameters.SetRingDim(65536);
            }
            this->cc = lbcrypto::GenCryptoContext(parameters);
            cc->Enable(PKE);
            if (!YAO)
                cc->Enable(KEYSWITCH);
            cc->Enable(LEVELEDSHE);
            cc->Enable(MULTIPARTY);

            serialize_sendall(cc);
        }
        else {
            deserialize_recv(cc, ALICE);
        }
    }

    void multiparty_keygen() {
        if (party == ALICE) {
            kp = cc->KeyGen();
            if (!kp.good())
                error("Bad key pair");
            serialize_send(kp.publicKey, party + 1);
        }
        else {
            deserialize_recv(pk, party - 1);

            kp = cc->MultipartyKeyGen(pk);
            if (!kp.good())
                error("Bad key pair");
            if (party != num_party)
                serialize_send(kp.publicKey, party + 1);
        }
        if (party == num_party) {
            pk = kp.publicKey;
            serialize_sendall(pk);
        }
        else {
            deserialize_recv(pk, num_party);
        }
    }

    void rotation_keygen() {
        //Rotation keys assuming tables are only size 2
        std::ifstream file("data/eval_auto_key_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                           std::ios::out | std::ios::binary);
        if (file.is_open()) {
            cc->DeserializeEvalAutomorphismKey(file, lbcrypto::SerType::BINARY);
            file.close();
            return;
        }
        std::vector<int32_t> indices = {1, -1};

        if (party == ALICE) {
            cc->EvalAtIndexKeyGen(kp.secretKey, indices);

            evalAtIndexKeys = std::make_shared<std::map<usint, lbcrypto::EvalKey<lbcrypto::DCRTPoly>>>(
                cc->GetEvalAutomorphismKeyMap(kp.secretKey->GetKeyTag()));
        }
        else {
            deserialize_recv(evalAtIndexKeys, party - 1);

            auto evalAtIndexKeys2 =
                cc->MultiEvalAtIndexKeyGen(kp.secretKey, evalAtIndexKeys, indices, kp.publicKey->GetKeyTag());
            evalAtIndexKeys =
                cc->MultiAddEvalAutomorphismKeys(evalAtIndexKeys, evalAtIndexKeys2, kp.publicKey->GetKeyTag());
        }

        if (party != num_party) {
            serialize_send(evalAtIndexKeys, party + 1);
            deserialize_recv(evalAtIndexKeys, num_party);
        }
        else {
            serialize_sendall(evalAtIndexKeys);
        }
        cc->InsertEvalAutomorphismKey(evalAtIndexKeys);
    }

    void multiplication_keygen() {
        std::ifstream file("data/eval_mult_key_" + std::to_string(party) + "_" + std::to_string(num_party) + ".txt",
                           std::ios::out | std::ios::binary);
        if (file.is_open()) {
            cc->DeserializeEvalMultKey(file, lbcrypto::SerType::BINARY);
            file.close();
            return;
        }
        lbcrypto::EvalKey<lbcrypto::DCRTPoly> evalMultKey, evalMultKeyShare, evalMultFinal;
        if (party == ALICE) {
            evalMultKeyShare = cc->KeySwitchGen(kp.secretKey, kp.secretKey);
            evalMultKey      = evalMultKeyShare;
            this->serialize_send(evalMultKey, party + 1);

            this->deserialize_recv(evalMultKey, num_party);
            evalMultKeyShare = cc->MultiMultEvalKey(kp.secretKey, evalMultKey, pk->GetKeyTag());
            evalMultFinal    = evalMultKeyShare;
            this->serialize_send(evalMultFinal, party + 1);
            this->deserialize_recv(evalMultFinal, num_party);
        }
        else {
            this->deserialize_recv(evalMultKey, party - 1);
            evalMultKeyShare = cc->MultiKeySwitchGen(kp.secretKey, kp.secretKey, evalMultKey);
            evalMultKey      = cc->MultiAddEvalKeys(evalMultKey, evalMultKeyShare, pk->GetKeyTag());

            if (party != num_party) {
                this->serialize_send(evalMultKey, party + 1);
                this->deserialize_recv(evalMultKey, num_party);
            }
            else
                this->serialize_sendall(evalMultKey);

            this->deserialize_recv(evalMultFinal, party - 1);
            evalMultKeyShare = cc->MultiMultEvalKey(kp.secretKey, evalMultKey, pk->GetKeyTag());
            evalMultFinal    = cc->MultiAddEvalMultKeys(evalMultFinal, evalMultKeyShare, evalMultKey->GetKeyTag());

            if (party != num_party) {
                this->serialize_send(evalMultFinal, party + 1);
                this->deserialize_recv(evalMultFinal, num_party);
            }
            else
                this->serialize_sendall(evalMultFinal);
        }

        cc->InsertEvalMultKey({evalMultFinal});
    }

    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> decrypt_partial(
        const std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& ciphertextVec) {
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> part_dec;
        if (party == ALICE) {
            part_dec = cc->MultipartyDecryptLead(ciphertextVec, kp.secretKey);
        }
        else {
            part_dec = cc->MultipartyDecryptMain(ciphertextVec, kp.secretKey);
        }
        return part_dec;
    }

    void enc_to_share(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& ciphertext, int64_t* share, uint n,
                      PlaintextEncodings encoding = PACKED_ENCODING) {
        uint batch_size = this->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
        int threads     = this->pool->size();
        if (ceil((double)n / (double)batch_size) != ciphertext.size()) {
            std::cout << "ciphertext size" + std::to_string(ciphertext.size()) + "is smaller than" +
                             std::to_string(ceil((double)n / (double)batch_size))
                      << std::endl;
            exit(1);
        }

        if (party == ALICE) {
            auto pt = cc->MultipartyDecryptLead(ciphertext, kp.secretKey);
            std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>> partial_decs;

            for (size_t i = 0; i < ciphertext.size(); ++i) {
                partial_decs.push_back({pt[i]});
                partial_decs[i].resize(num_party, pt[i]);
            }
            vector<std::future<void>> res;
            for (int i = 2; i <= num_party; ++i) {
                res.push_back(pool->enqueue([this, &partial_decs, i, x = ciphertext.size()]() {
                    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> partial_decs_i;
                    deserialize_recv(partial_decs_i, i);

                    for (int j = 0; j < x; ++j) {
                        partial_decs[j][i - 1] = partial_decs_i[j];
                    }
                    partial_decs_i.clear();
                }));
            }

            for (auto& v : res)
                v.get();
            res.clear();

            int num_steps = ceil((double)ciphertext.size() / ((double)threads));
            for (int i = 0; i < threads; ++i) {
                res.push_back(pool->enqueue([this, i, &partial_decs, x = ciphertext.size(), num_steps, batch_size,
                                             share, encoding, n]() {
                    for (int j = 0; j < num_steps; ++j) {
                        if (i * num_steps + j >= x)
                            break;
                        lbcrypto::Plaintext ptxt;
                        cc->MultipartyDecryptFusion(partial_decs[i * num_steps + j], &ptxt);
                        // std::cout << ptxt << std::endl;
                        vector<int64_t> tmp;
                        if (encoding == COEF_PACKED_ENCODING)
                            tmp = ptxt->GetCoefPackedValue();
                        else
                            tmp = ptxt->GetPackedValue();

                        if ((i + 1) * batch_size <= n)
                            memcpy(share + (i * num_steps + j) * batch_size, tmp.data(), batch_size * sizeof(int64_t));
                        else
                            memcpy(share + (i * num_steps + j) * batch_size, tmp.data(),
                                   (n % batch_size) * sizeof(int64_t));

                        tmp.clear();
                    }
                }));
            }

            for (auto& v : res)
                v.get();
            res.clear();
        }
        else {
            this->prg.random_data(share, n * sizeof(int64_t));
            for (int i = 0; i < n; ++i) {
                share[i] = ((share[i] % this->q) + this->q) % this->q;
                if (encoding == COEF_PACKED_ENCODING) {
                    if (share[i] > this->q / 2)
                        share[i] -= this->q;
                }
            }
            auto partial_dec_i = cc->MultipartyDecryptMain(ciphertext, kp.secretKey);
            vector<std::future<void>> res;
            uint64_t num_steps = ceil((double)ciphertext.size() / ((double)threads));
            for (uint64_t i = 0; i < threads; ++i) {
                res.push_back(pool->enqueue([this, i, &partial_dec_i, x = ciphertext.size(), num_steps, batch_size,
                                             share, encoding, n]() {
                    for (uint64_t j = 0; j < num_steps; ++j) {
                        if (i * num_steps + j >= x)
                            break;

                        vector<int64_t> tmp;
                        tmp.resize(batch_size);
                        if ((i * num_steps + j + 1) * batch_size <= n)
                            memcpy(tmp.data(), share + (i * num_steps + j) * batch_size, batch_size * sizeof(int64_t));
                        else
                            memcpy(tmp.data(), share + (i * num_steps + j) * batch_size,
                                   (n % batch_size) * sizeof(int64_t));

                        lbcrypto::Plaintext ptxt;
                        if (encoding == COEF_PACKED_ENCODING)
                            ptxt = cc->MakeCoefPackedPlaintext(tmp);
                        else
                            ptxt = cc->MakePackedPlaintext(tmp);
                        // std::cout << ptxt << std::endl;
                        partial_dec_i[i * num_steps + j] = cc->EvalSub(partial_dec_i[i * num_steps + j], ptxt);
                        tmp.clear();
                    }
                }));
            }

            for (auto& v : res)
                v.get();
            res.clear();
            serialize_send(partial_dec_i, ALICE);
        }
    }

    void bootstrap(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ciphertext, uint holding_party, uint final_party,
                   PlaintextEncodings encoding = PACKED_ENCODING) {
        vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> cts;
        cts.push_back(ciphertext);
        bootstrap(cts, holding_party, final_party, encoding);
        ciphertext = cts[0];
    }

    void bootstrap(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& ciphertext, uint holding_party,
                   uint final_party, PlaintextEncodings encoding = PACKED_ENCODING) {
        uint batch_size = this->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
        // std::cout << "bootstrap \t" << party << ": " << holding_party << " " << final_party << std::endl;
        if (this->party == holding_party) {
            this->serialize_sendall(ciphertext, 1, BOOT_REQ_MSG);
        }
        else {
            this->deserialize_recv(ciphertext, holding_party, 1, BOOT_REQ_MSG);
        }
        // std::cout <<  party << "\t request received \n";
        int n = ciphertext.size();
        if (party != final_party) {
            int64_t* share = new int64_t[n * batch_size];
            this->prg.random_data(share, n * batch_size * sizeof(int64_t));
            for (int i = 0; i < n * batch_size; ++i) {
                share[i] = ((share[i] % this->q) + this->q) % this->q;
                if (encoding == COEF_PACKED_ENCODING) {
                    if (share[i] > this->q / 2)
                        share[i] -= this->q;
                }
            }

            std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> partial_decs_encs = this->decrypt_partial(ciphertext);

            vector<int64_t> tmp;
            tmp.resize(batch_size);
            for (int j = 0; j < n; ++j) {
                memcpy(tmp.data(), share + (j)*batch_size, batch_size * sizeof(int64_t));
                lbcrypto::Plaintext ptxt;
                if (encoding == COEF_PACKED_ENCODING)
                    ptxt = cc->MakeCoefPackedPlaintext(tmp);
                else
                    ptxt = cc->MakePackedPlaintext(tmp);
                partial_decs_encs[j] = cc->EvalSub(partial_decs_encs[j], ptxt);

                // // Append encryption of the share to end
                lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c = cc->Encrypt(this->pk, ptxt);
                partial_decs_encs.push_back(c);
            }
            tmp.clear();

            this->serialize_send(partial_decs_encs, final_party, 0, BOOT_RSP_MSG);
            // std::cout << party << "\t partial_decs_encs sent \n";
            delete[] share;
        }
        else {
            auto pt = decrypt_partial(ciphertext);
            std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>> partial_decs;
            std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> new_ciphertext;
            for (int i = 0; i < n; ++i) {
                partial_decs.push_back({pt[i]});
                partial_decs[i].resize(num_party, pt[i]);
            }
            std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> partial_decs_encs;

            for (int i = 1; i <= num_party; ++i) {
                if (i != final_party) {
                    deserialize_recv(partial_decs_encs, i, 0, BOOT_RSP_MSG);

                    for (int j = 0; j < n; ++j) {
                        partial_decs[j][i - 1] = partial_decs_encs[j];
                        if (i == 1)
                            new_ciphertext.push_back(partial_decs_encs[ciphertext.size() + j]);
                        else
                            cc->EvalAddInPlace(new_ciphertext[j], partial_decs_encs[ciphertext.size() + j]);
                    }
                }
            }
            // std::cout << party << "\t partial_decs_encs received \n";
            for (int j = 0; j < n; ++j) {
                lbcrypto::Plaintext ptxt;
                cc->MultipartyDecryptFusion(partial_decs[j], &ptxt);
                vector<int64_t> tmp;
                if (encoding == COEF_PACKED_ENCODING)
                    tmp = ptxt->GetCoefPackedValue();
                else
                    tmp = ptxt->GetPackedValue();

                if (encoding == COEF_PACKED_ENCODING)
                    ptxt = cc->MakeCoefPackedPlaintext(tmp);
                else
                    ptxt = cc->MakePackedPlaintext(tmp);
                ciphertext[j] = cc->Encrypt(this->pk, ptxt);
                cc->EvalAddInPlace(ciphertext[j], new_ciphertext[j]);
            }
        }
    }
};
}  // namespace emp
