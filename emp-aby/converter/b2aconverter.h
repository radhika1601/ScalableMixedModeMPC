#pragma once

#include "emp-aby/lut.h"
#include "emp-aby/mp-circuit.hpp"
#include "emp-aby/simd_interface/arithmetic-circ.h"

namespace emp {

template <typename IO>
class B2AConverter {
private:
    LUT<IO>* bit_to_a;
    MPIOChannel<IO>* io;
    int num_party, party;
    ThreadPool* pool;
    HE<IO>* he;

public:
    B2AConverter(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int pool_size);

    /**
             * @brief This function converts bit vector for boolean shares to arithmetic shares
             * 
             * @param out Output arithmetic share array
             * @param in Input array of boolean shares
             * @param length length of boolean array
             * @param l length of each bit-vector
             */
    void convert(int64_t* out, bool* in, size_t length, int l);
};

template <typename IO>
B2AConverter<IO>::B2AConverter(int num_party, int party, MPIOChannel<IO>* io, ThreadPool* pool, HE<IO>* he, int pool_size) {
    int64_t table[2] = {0, 1};
    this->pool       = pool;
    this->io         = io;
    this->num_party  = num_party;
    this->party      = party;
    this->he         = he;
    // std::cout << "go in lut constructor \n";
    this->bit_to_a = new LUT<IO>(num_party, party, io, pool, he, table, pool_size);
}

template <typename IO>
void B2AConverter<IO>::convert(int64_t* out, bool* in, size_t length, int l) {
    if (length % l != 0)
        error("Length of boolean array is not divisible by length of each bit vector.");
    size_t n           = length / l;
    int64_t* in_ashare = new int64_t[length];
    this->bit_to_a->lookup(in_ashare, in, length);
    memset(out, 0, n * sizeof(int64_t));

    size_t threads   = pool->size();
    size_t num_steps = ceil((double)n / (double)threads);
    vector<std::future<void>> res;
    for (size_t t = 0; t < threads; ++t) {
        res.push_back(pool->enqueue([this, n, t, num_steps, in_ashare, out, l]() {
            for (size_t step = 0; step < num_steps; ++step) {
                size_t j = t * num_steps + step;
                if (j >= n)
                    break;

                for (size_t i = 0; i < l; ++i) {
                    uint64_t x;
                    in_ashare[j * l + i] = (in_ashare[j * l + i] % he->q + he->q) % he->q;
                    x                    = in_ashare[j * l + i];
                    x                    = x << i;
                    x %= he->q;
                    out[j] = (out[j] + x) % he->q;
                }
            }
        }));
    }
    for (auto& fut : res)
        fut.get();
    res.clear();
}

}  // namespace emp
