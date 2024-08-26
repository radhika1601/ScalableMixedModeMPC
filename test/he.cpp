#include "emp-aby/he_interface.hpp"
#include "emp-aby/io/multi-io.hpp"
using namespace emp;

#include <typeinfo>
#include <math.h>
int party, port;

const static int threads = 4;

int num_party;

template <typename IO>
void test_decrypt(HE<IO>* he, MPIOChannel<IO>* io) {
    int n = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;

    if (party == ALICE) {
        std::vector<int64_t> original(n);
        PRG prg;
        prg.random_data(original.data(), n * sizeof(int64_t));
        for (int i = 0; i < n; ++i) {
            original[i] = ((original[i] % he->q) + he->q) % 2;
        }

        lbcrypto::Plaintext plaintext = he->cc->MakePackedPlaintext(original);

        auto ciphertext = he->cc->Encrypt(he->pk, plaintext);

        he->serialize_sendall(ciphertext);
        auto partial_dec = he->decrypt_partial({ciphertext});

        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> partialCiphertextVec;
        partialCiphertextVec.push_back(partial_dec[0]);
        for (int i = 2; i <= num_party; ++i) {
            he->deserialize_recv(partial_dec, i);
            partialCiphertextVec.push_back(partial_dec[0]);
        }

        lbcrypto::Plaintext new_ptxt;
        he->cc->MultipartyDecryptFusion(partialCiphertextVec, &new_ptxt);
        // new_ptxt->SetLength(plaintext->GetLength());

        vector<int64_t> dec_vals = new_ptxt->GetPackedValue();
        for (int i = 0; i < n; ++i) {
            if (dec_vals[i] != original[i])
                if ((dec_vals[i] + he->q) % he->q != original[i]) {
                    std::cout << (dec_vals[i] + he->q) % he->q << " " << original[i] << "\n";
                    error("Decryption not working!");
                }
        }
        std::cout << "Decryption test passed!" << std::endl;
    }
    else {
        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> ciphertext;
        he->deserialize_recv(ciphertext, ALICE);
        auto partial_dec = he->decrypt_partial({ciphertext});

        he->serialize_send(partial_dec, ALICE);
    }
}

template <typename IO>
void test_rotation(HE<IO>* he, MPIOChannel<IO>* io) {
    int n = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    if (party == ALICE) {
        std::vector<int64_t> original(n);
        PRG prg;
        prg.random_data(original.data(), n * sizeof(int64_t));
        for (int i = 0; i < n; ++i) {
            original[i] = ((original[i] % he->q) + he->q) % 2;
        }

        lbcrypto::Plaintext plaintext = he->cc->MakePackedPlaintext(original);
        auto ciphertext               = he->cc->Encrypt(he->pk, plaintext);
        auto permutedCiphertext       = he->cc->EvalRotate(ciphertext, 1);
        std::cout << "send to all\n";
        he->serialize_sendall(permutedCiphertext);

        int64_t* share                                            = new int64_t[n];
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> vec = {ciphertext};
        he->enc_to_share(vec, share, n);

        for (int i = 2; i <= he->num_party; ++i) {
            std::vector<int64_t> tmp;
            tmp.resize(n);
            io->recv_data(i, (int64_t*)tmp.data(), tmp.size() * sizeof(int64_t));
            for (int j = 0; j < n; ++j) {
                share[j] = (he->q + share[j] + tmp[j]) % he->q;
            }
        }
        std::cout << "original\t";
        for (int i = 0; i < 5; ++i) {
            std::cout << original[i] << " ";
        }
        std::cout << "\n";

        std::cout << "rotated\t";
        for (int i = 0; i < 5; ++i) {
            std::cout << share[i] << " ";
        }
        std::cout << "\n";
        std::cout << "Rotation test passed!" << std::endl;
    }
    else {
        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> ciphertext;
        he->deserialize_recv(ciphertext, ALICE);
        std::cout << "receive ciphertext\n";
        int64_t* share                                            = new int64_t[n];
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> vec = {ciphertext};
        he->enc_to_share(vec, share, n);
        io->send_data(ALICE, (int64_t*)share, n * sizeof(int64_t));
    }
}

template <typename IO>
void test_enc_to_share(HE<IO>* he, MPIOChannel<IO>* io, int n = 100) {
    PRG prg;
    int batch_size    = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    n                 = batch_size;
    int64_t* original = new int64_t[n];
    prg.random_data(original, n * sizeof(int64_t));
    for (int i = 0; i < n; ++i) {
        original[i] = ((original[i] % he->q) + he->q) % he->q;
    }

    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> ciphertext;

    if (party == ALICE) {
        for (int i = 0; i < ceil((double)n / (double)batch_size); ++i) {
            vector<int64_t> tmp;
            if ((i + 1) * batch_size <= n) {
                tmp.resize(batch_size);
                memcpy(tmp.data(), original + i * batch_size, batch_size * sizeof(int64_t));
            }
            else {
                tmp.resize(n % batch_size);
                memcpy(tmp.data(), original + i * batch_size, (n % batch_size) * sizeof(int64_t));
            }

            auto plaintext = he->cc->MakePackedPlaintext(tmp);
            // std::cout << "Plaintext: " << plaintext << std::endl;
            ciphertext.push_back(he->cc->Encrypt(he->pk, plaintext));
        }
        he->serialize_sendall(ciphertext);

        int64_t* share = new int64_t[n];
        he->enc_to_share(ciphertext, share, n);

        for (int i = 2; i <= he->num_party; ++i) {
            std::vector<int64_t> tmp;
            tmp.resize(n);
            io->recv_data(i, (int64_t*)tmp.data(), tmp.size() * sizeof(int64_t));
            for (int j = 0; j < n; ++j) {
                share[j] = (he->q + share[j] + tmp[j]) % he->q;
            }
        }

        for (int i = 0; i < n; ++i) {
            if (original[i] != share[i]) {
                std::cout << i << " " << original[i] << " " << share[i] << std::endl;
                error("enc_to_share Failed!");
            }
        }
        std::cout << "enc_to_share test passed!" << std::endl;
    }
    else {
        he->deserialize_recv(ciphertext, ALICE);

        he->enc_to_share(ciphertext, original, n);
        io->send_data(ALICE, (int64_t*)original, n * sizeof(int64_t));
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Format: test_mp_bit_triple PartyID Port num_parties" << std::endl;
        exit(0);
    }
    parse_party_and_port(argv, &party, &port);
    num_party = atoi(argv[3]);

    std::vector<std::pair<std::string, unsigned short>> net_config;

    if (argc == 5) {
        const char* file = argv[4];
        FILE* f          = fopen(file, "r");
        for (int i = 0; i < num_party; ++i) {
            char* c = (char*)malloc(15 * sizeof(char));
            uint p;
            fscanf(f, "%s %d\n", c, &p);
            std::string s(c);
            net_config.push_back(std::make_pair(s, p));
            fflush(f);
        }
        fclose(f);
    }
    else {
        for (int i = 0; i < num_party; ++i) {
            std::string s = "127.0.0.1";
            uint p        = (port + 4 * num_party * i);
            net_config.push_back(std::make_pair(s, p));
        }
    }

    MultiIO* io = new MultiIO(party, num_party, net_config);
    std::cout << party << " connected \n";
    ThreadPool pool(threads);

    io->flush();
    const long long int modulus = (1L << 32) - (1 << 30) + 1;
    HE<MultiIOBase>* he         = new HE<MultiIOBase>(num_party, io, &pool, party, modulus);
    he->multiplication_keygen();
    he->rotation_keygen();

    std::cout << party << " p = " << he->cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << party << " n = " << he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    std::cout << party
              << " log2 q = " << log2(he->cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    vector<std::future<void>> res;
    // uint batch_size = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    test_decrypt(he, io);
    test_enc_to_share(he, io);
    test_rotation(he, io);
    delete io;
}
