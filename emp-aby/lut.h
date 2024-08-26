#pragma once

#include "emp-aby/he_interface.hpp"
#include "emp-aby/utils.h"

#include <cmath>
#include <poll.h>
namespace emp {

template <typename IO>
class LUT {
private:
    int rotated_pool_size;
    int num_used = 0;
    ThreadPool* pool;
    bool* rotation;
//    int64_t* lut_share;
    HE<IO>* he;
    MPIOChannel<IO>* io;
    PRG prg;

    void shuffle(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& c, bool* rotation, size_t batch_size, size_t i);

public:
    int64_t *lut_share;
    int num_party;
    int party;
    int64_t table[2];
    LUT(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int rot_pool_size = 20);
    LUT(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int64_t table[2],
        int rot_pool_size = 20);
    ~LUT();
    void generate_shares(int64_t* lut_share, bool* rotation, int num_shares, int64_t table[2]);
    void lookup(int64_t* out, bool* in, size_t length);
};

template <typename IO>
LUT<IO>::LUT(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int rot_pool_size) {
    this->io        = io;
    this->party     = party;
    this->num_party = num_party;
    this->pool      = pool;

    this->he = he;
    this->rotated_pool_size =
        rot_pool_size * (he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2);
}

template <typename IO>
LUT<IO>::LUT(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int64_t table[2],
             int rot_pool_size)
    : LUT(num_party, party, io, pool, he, rot_pool_size) {
    this->table[0]  = table[0];
    this->table[1]  = table[1];
    this->rotation  = new bool[rotated_pool_size];
    this->lut_share = new int64_t[2 * rotated_pool_size];
    // std::cout << "go in gen shares \n";
    this->generate_shares(this->lut_share, this->rotation, this->rotated_pool_size, this->table);
}

template <typename IO>
void LUT<IO>::shuffle(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& c, bool* rotation, size_t batch_size, size_t rot_idx) {
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> rot_1, rot_2;
    rot_1 = he->cc->EvalRotate(c, 1);
    rot_2 = he->cc->EvalRotate(c, -1);
    vector<int64_t> tmp;
    vector<int64_t> mult1, mult2, mult3;
    mult1.resize(batch_size);
    mult2.resize(batch_size);
    mult3.resize(batch_size);
    for (int j = 0; j < batch_size / 2; ++j) {
        if (rotation[(rot_idx)*batch_size / 2 + j] == true) {
            mult1[2 * j]     = 1;
            mult2[2 * j + 1] = 1;
        }
        else {
            mult3[2 * j]     = 1;
            mult3[2 * j + 1] = 1;
        }
    }

    auto plain1 = he->cc->MakePackedPlaintext(mult3);
    // auto tmp1 = he->cc->Encrypt(he->pk, plain1);
    auto tmp1 = he->cc->EvalMult(c, plain1);
    c         = tmp1;

    plain1 = he->cc->MakePackedPlaintext(mult1);
    // tmp1 = he->cc->Encrypt(he->pk, plain1);
    tmp1 = he->cc->EvalMult(rot_1, plain1);
    he->cc->EvalAddInPlace(c, tmp1);

    plain1 = he->cc->MakePackedPlaintext(mult2);
    // tmp1 = he->cc->Encrypt(he->pk, plain1);
    tmp1 = he->cc->EvalMult(rot_2, plain1);
    he->cc->EvalAddInPlace(c, tmp1);
}

template <typename IO>
void LUT<IO>::generate_shares(int64_t* lut_share, bool* rotation, int num_shares, int64_t* table) {
    int batch_size = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    int n          = 2 * num_shares;
    prg.random_bool((bool*)rotation, num_shares);
    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> ciphertext;
    vector<std::future<void>> res;
    res.push_back(pool->enqueue([this, n, rotation, batch_size, &ciphertext, table]() {
        // auto start = clock_start();
        for (int i = 0; i < ceil((double)n / (double)batch_size); ++i) {
            lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c;

            // for (int j = 1; j < party; ++j) {
            if (((party - 1) >= he->mult_depth) && (((party - 1) % he->mult_depth) == 0)) {
                he->bootstrap(c, party - 1, party);
            }
            // }
            if (party == ALICE) {
                vector<int64_t> tmp;
                tmp.resize(batch_size);
                for (int j = 0; j < batch_size / 2; ++j) {
                    if (rotation[(i * batch_size) / 2 + j] == true) {
                        tmp[2 * j]     = table[1];
                        tmp[2 * j + 1] = table[0];
                    }
                    else {
                        tmp[2 * j]     = table[0];
                        tmp[2 * j + 1] = table[1];
                    }
                }

                lbcrypto::Plaintext plaintext = he->cc->MakePackedPlaintext(tmp);
                c                             = he->cc->Encrypt(he->pk, plaintext);
            }
            else {
                if (!(party > he->mult_depth) || !((party - 1) % he->mult_depth == 0)) {
                    he->deserialize_recv(c, party - 1);
                }
                shuffle(c, rotation, batch_size, i);
                // Refresh ciphertext
                const std::vector<lbcrypto::DCRTPoly>& cv = c->GetElements();
                const auto cryptoParams = std::dynamic_pointer_cast<lbcrypto::CryptoParametersRLWE<lbcrypto::DCRTPoly>>(
                    he->pk->GetCryptoParameters());
                const auto ns = cryptoParams->GetNoiseScale();

                lbcrypto::DCRTPoly::DggType dgg(NOISE_FLOODING::MP_SD);
                lbcrypto::DCRTPoly e(dgg, cv[0].GetParams(), Format::EVALUATION);
                lbcrypto::DCRTPoly b = cv[0] + ns * e;
                c->SetElements({std::move(b), std::move(cv[1])});
            }

            if (party != num_party) {
                if ((party < he->mult_depth) || (party % he->mult_depth != 0)) {
                    he->serialize_send(c, party + 1);
                }
                else {
                    he->bootstrap(c, party, party + 1);
                }
            }
            if (party == num_party)
                ciphertext.push_back(c);
        }
        // double timeused = time_from(start);
        // std::cout << party << "\tprocessing time\t" << timeused / 1000 << std::endl;
    }));

    for (int j = 1; j < num_party; ++j) {
        if ((j != party) && (j != party - 1)) {
            if ((j >= he->mult_depth) && ((j % he->mult_depth) == 0)) {
                res.push_back(pool->enqueue([this, j, batch_size, n]() {
                    for (int i = 0; i < ceil((double)n / (double)batch_size); ++i) {
                        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c;
                        he->bootstrap(c, j, j + 1);
                    }
                }));
            }
        }
    }

    for (auto& v : res)
        v.get();
    res.clear();

    if (party != num_party)
        he->deserialize_recv(ciphertext, num_party);
    else
        he->serialize_sendall(ciphertext);

    he->enc_to_share(ciphertext, lut_share, n);
}

template <typename IO>
void LUT<IO>::lookup(int64_t* out, bool* in, size_t length) {
    bool* r           = nullptr;
    int64_t* t        = nullptr;
    bool delete_array = false;
    if (length > rotated_pool_size) {
        r            = new bool[(length + rotated_pool_size - 1) / rotated_pool_size * rotated_pool_size];
        t            = new int64_t[(length + rotated_pool_size - 1) / rotated_pool_size * rotated_pool_size * 2];
        delete_array = true;
        for (uint i = 0; i < (length + rotated_pool_size - 1) / rotated_pool_size; ++i)
            this->generate_shares(t + 2 * i * rotated_pool_size, r + i * rotated_pool_size, rotated_pool_size,
                                  this->table);
        size_t tocp = std::min((int)((length + rotated_pool_size - 1) / rotated_pool_size * rotated_pool_size - length),
                               num_used);
        memcpy(this->rotation, r + (length + rotated_pool_size - 1) / rotated_pool_size * rotated_pool_size - length,
               tocp);
        memcpy(this->lut_share,
               t + 2 * ((length + rotated_pool_size - 1) / rotated_pool_size * rotated_pool_size - length),
               2 * tocp * sizeof(int64_t));
        num_used     = 0;
        delete_array = true;
    }
    else if (length > rotated_pool_size - num_used) {
        r            = new bool[length];
        t            = new int64_t[2 * length];
        delete_array = true;
        memcpy(r, rotation + num_used, rotated_pool_size - num_used);
        memcpy(t, lut_share + 2 * num_used, 2 * (rotated_pool_size - num_used) * sizeof(int64_t));
        this->generate_shares(this->lut_share, this->rotation, this->rotated_pool_size, this->table);
        memcpy(r + rotated_pool_size - num_used, this->rotation, length - (rotated_pool_size - num_used));
        memcpy(t + 2 * (rotated_pool_size - num_used), this->lut_share,
               2 * (length - (rotated_pool_size - num_used)) * sizeof(int64_t));
        num_used = length - (rotated_pool_size - num_used);
    }
    else {
        r = rotation + num_used;
        t = lut_share + 2 * num_used;
        num_used += length;
    }
    bool* e = new bool[length];
    xorBools_arr(e, r, in, length);
    if (party == ALICE) {
        bool* tmp = new bool[length];
        for (uint i = 2; i <= num_party; ++i) {
            io->recv_bool(i, tmp, length);
            xorBools_arr(e, e, tmp, length);
        }
        vector<std::future<void>> res;
        for (uint i = 2; i <= num_party; ++i) {
            res.push_back(pool->enqueue([this, i, e, length]() {
                this->io->send_bool(i, e, length);
                io->flush(i);
            }));
        }
        for (auto& v : res)
            v.get();
        res.clear();

        delete[] tmp;
    }
    else {
        io->send_bool(ALICE, e, length);
        io->flush(ALICE);
        io->recv_bool(ALICE, e, length);
    }
    for (int i = 0; i < length; ++i) {
        out[i] = t[i * 2 + (int)e[i]];
    }

    delete[] e;
    if (delete_array) {
        delete[] r;
        delete[] t;
    }
}

template <typename IO>
LUT<IO>::~LUT() {
    delete[] rotation;
    delete[] lut_share;
}

}  // namespace emp
