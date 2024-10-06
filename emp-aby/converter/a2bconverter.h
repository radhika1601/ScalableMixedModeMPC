#pragma once

#include "emp-aby/lut.h"
#include "emp-aby/mp-circuit.hpp"
#include "emp-aby/simd_interface/arithmetic-circ.h"

namespace emp {

template <typename IO>
class A2BConverter {
private:
    LUT<IO>* bit_to_a;
    MPIOChannel<IO>* io;
    int num_party, party;
    ThreadPool* pool;
    HE<IO>* he;
    PRG prg;
    Circuit<MPSIMDCircExec<IO>>* circuit;
    ArithmeticCirc<IO>* arithmetic_circ;
    bool* b_share;
    int64_t* a_share;

public:
    size_t ab_share_pool, num_used_shares = 0, num_rejected = 0;
    A2BConverter(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he,
                 MPSIMDCircExec<IO>* simd_circ, int pool_size = 20);

    void convert(bool* out, int64_t* in, size_t length);

    size_t rand_ab_shares(int64_t* a_share, bool* b_share, const size_t length);
};

template <typename IO>
A2BConverter<IO>::A2BConverter(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he,
                               MPSIMDCircExec<IO>* simd_circ, int pool_size) {
    int64_t table[2] = {0, 1};
    this->pool       = pool;
    this->io         = io;
    this->num_party  = num_party;
    this->party      = party;
    this->he         = he;
    this->bit_to_a   = new LUT<IO>(num_party, party, io, pool, he, table, pool_size);
    // std::cout << "lut setup done \n";
    this->circuit         = new Circuit<MPSIMDCircExec<IO>>("emp-aby/modsum.txt", party, simd_circ);
    this->arithmetic_circ = new ArithmeticCirc<IO>(num_party, party, io, he);
    // std::cout << "arithmetic circ setup done \n";
    int l               = ceil(log2(he->q));
    this->ab_share_pool = pool_size * (he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2) / l;
    this->a_share       = new int64_t[ab_share_pool];
    this->b_share       = new bool[ab_share_pool * (int)ceil(log2(he->q))];
    this->num_rejected  = this->rand_ab_shares(this->a_share, this->b_share, this->ab_share_pool);
    // std::cout << "random ab shares setup done \n";
}

template <typename IO>
size_t A2BConverter<IO>::rand_ab_shares(int64_t* a_share, bool* b_share, const size_t length) {
    int64_t l    = ceil(log2(he->q));
    bool* r_b    = new bool[length * l];
    int64_t* r_a = new int64_t[length * l];
    prg.random_bool(r_b, length * l);

    this->bit_to_a->lookup(r_a, r_b, length * l);
    // std::cout << "lookup done\n";
    for (int i = 0; i < length * l; ++i) {
        r_a[i] = (r_a[i] % he->q + he->q) % he->q;
    }

    bool* bits = new bool[l];
    int x      = he->q - 1;
    for (int i = 0; i < l; ++i) {
        bits[i] = ((x & 1) == 1);
        x >>= 1;
    }
    int64_t* zero_sum = new int64_t[length];
    memset(zero_sum, 0, length * sizeof(length));
    for (int i = 0; i < l; ++i) {
        if (!bits[i]) {
            for (int j = 0; j < length; ++j) {
                zero_sum[j] = (zero_sum[j] + r_a[j * l + i]) % he->q;
            }
        }
    }
    int64_t* tmp = new int64_t[length];
    for (int i = 0; i < l; ++i) {
        if (bits[i]) {
            for (int j = 0; j < length; ++j) {
                tmp[j] = r_a[j * l + i];
            }
            this->arithmetic_circ->mult(zero_sum, zero_sum, tmp, length);
        }
    }
    if (party == ALICE) {
        for (int i = 2; i <= num_party; ++i) {
            io->recv_data(i, tmp, length * sizeof(int64_t));
            for (int j = 0; j < length; ++j)
                zero_sum[j] = (zero_sum[j] + tmp[j]) % he->q;
        }

        for (int i = 2; i <= num_party; ++i) {
            io->send_data(i, zero_sum, length * sizeof(int64_t));
            io->flush(i);
        }
    }
    else {
        io->send_data(ALICE, zero_sum, length * sizeof(int64_t));
        io->flush(ALICE);

        io->recv_data(ALICE, zero_sum, length * sizeof(int64_t));
    }
    memset(a_share, 0, length * sizeof(int64_t));

    size_t waste = 0;
    uint64_t y   = 0;
    for (int i = 0; i < length; ++i) {
        if (zero_sum[i] == 0) {
            for (int j = 0; j < l; ++j) {
                y = r_a[i * l + j];
                y = y << j;
                y %= he->q;
                a_share[i - waste]           = (a_share[i - waste] + y) % he->q;
                b_share[(i - waste) * l + j] = r_b[i * l + j];
            }
        }
        else
            waste += 1;
    }

    delete[] tmp;
    return waste;
}

template <typename IO>
void A2BConverter<IO>::convert(bool* out, int64_t* in, size_t length) {
    int64_t l = ceil(log2(he->q));
    bool* r_b;
    int64_t* r_a;
    bool delete_array = false;

    if (length > ab_share_pool - num_used_shares - num_rejected) {
        r_a = new int64_t[length];
        r_b = new bool[length * l];
        memcpy(r_a, this->a_share + num_used_shares,
               (ab_share_pool - num_used_shares - num_rejected) * sizeof(int64_t));
        memcpy(r_b, this->b_share + num_used_shares * l, (ab_share_pool - num_used_shares - num_rejected) * l);
        int accepted = ab_share_pool - num_used_shares - num_rejected;
        while (accepted < length) {
            num_rejected = this->rand_ab_shares(a_share, b_share, ab_share_pool);
            uint tocp    = min(ab_share_pool - num_rejected, length - accepted);
            memcpy(r_a + accepted, this->a_share, tocp * sizeof(int64_t));
            memcpy(r_b + accepted * l, this->b_share, tocp * l);
            accepted += ab_share_pool - num_rejected;
        }
        num_used_shares = ab_share_pool - (accepted - length) - num_rejected;
        delete_array    = true;
    }
    else {
        r_a = a_share + num_used_shares;
        r_b = b_share + num_used_shares * l;
        num_used_shares += length;
    }
    // std::cout << "AB shares generated" << std::endl;
    bool* y_b = new bool[length * (circuit->n1 + circuit->n2)];
    memset(y_b, 0, length * (circuit->n1 + circuit->n2));  // x-r (length * n), r (length * n)

    for (size_t i = 0; i < length; ++i) {
        for (int j = 0; j < l; ++j) {
            y_b[i * circuit->n2 + j + length * circuit->n1] = r_b[i * l + j];
        }
    }

    for (size_t i = 0; i < length; ++i) {
        r_a[i] = (in[i] - r_a[i]) % he->q;
        r_a[i] = (he->q + r_a[i]) % he->q;
    }

    // std::cout << "Do in - r_a" << std::endl;
    if (party == ALICE) {
        if (pool->size() == 1) {
            int64_t* tmp = new int64_t[length];
            for (int i = 2; i <= num_party; ++i) {
                io->recv_data(i, tmp, length * sizeof(int64_t));
                io->flush(i);
                for (size_t j = 0; j < length; ++j)
                    r_a[j] = (tmp[j] + r_a[j]) % he->q;
            }
        }
        else {
            int threads = pool->size();
            int64_t* tmp[threads];
            for (int i = 0; i < threads; ++i) {
                tmp[i] = new int64_t[length];
                memset(tmp[i], 0, length * sizeof(int64_t));
            }
            vector<std::future<void>> res;
            int num_steps = ceil((double)(num_party) / (double)threads);
            for (int i = 0; i < threads; ++i) {
                res.push_back(pool->enqueue([this, i, num_steps, length, t = tmp[i]]() {
                    int64_t* tmp_i = new int64_t[length];

                    for (int j = 0; j < num_steps; ++j) {
                        if (i * num_steps + j + 1 > num_party)
                            break;
                        if (i * num_steps + j + 1 == party)
                            continue;

                        io->recv_data(i * num_steps + j + 1, tmp_i, length * sizeof(int64_t));
                        io->flush(i * num_steps + j + 1);
                        for (size_t k = 0; k < length; ++k)
                            t[k] = (t[k] + tmp_i[k]) % he->q;
                    }
                    delete[] tmp_i;
                }));
            }

            for (auto& v : res)
                v.get();
            res.clear();

            for (int i = 0; i < threads; ++i) {
                for (size_t j = 0; j < length; ++j)
                    r_a[j] = (tmp[i][j] + r_a[j]) % he->q;
                delete[] tmp[i];
            }
        }

        for (size_t i = 0; i < length; ++i) {
            int64_t x = r_a[i];
            for (int j = 0; j < circuit->n1; ++j) {
                y_b[i * circuit->n1 + j] = ((x & 1) == 1);
                x >>= 1;
            }
        }

        //            delete[] tmp;
    }
    else {
        this->io->send_data(ALICE, r_a, length * sizeof(int64_t));
        io->flush(ALICE);
    }
    // std::cout << "Circuit to be computed \n";
    // std::cout << length *circuit->n3 << std::endl;
    bool* tmp_out = new bool[length * circuit->n3];
    circuit->template compute<IO>(tmp_out, y_b, length);
    // std::cout << "Circuit computed \n";

    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < l; ++j) {
            out[i * l + j] = tmp_out[i * circuit->n3 + j];
        }
    }
    // std::cout << "Circuit computed \n";

    delete[] y_b;
    delete[] tmp_out;
    if (delete_array) {
        delete[] r_a;
        delete[] r_b;
    }
}

}  //namespace emp
