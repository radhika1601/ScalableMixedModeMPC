#pragma once

#include "simd.h"
#include "emp-aby/triple-providers/mp-bit-triple.h"

namespace emp {

template <typename IO>
class MPSIMDCircExec: SIMDCircuitExecution {
private:
    MPBitTripleProvider<IO>* btp = nullptr;
    MPIOChannel<IO>* io;
    ThreadPool* pool;
    bool *bit_triple_a, *bit_triple_b, *bit_triple_c;
    block *block_triple_a, *block_triple_b, *block_triple_c;
    size_t num_triples = 0, num_block_triples = 0;
    size_t num_triples_pool;
    size_t num_block_triples_pool;
    int num_party;
    double total_time = 0;
    vector<block> d, d1, e, e1;

public:
    int cur_party;

    MPSIMDCircExec(int num_party, int party, ThreadPool* pool, MPIOChannel<IO>* io) {
        this->cur_party = party;
        this->num_party = num_party;
        this->pool      = pool;
        btp             = new MPBitTripleProvider<IO>(num_party, party, pool, io);
        num_triples_pool       = btp->BUFFER_SZ;
        num_block_triples_pool = btp->BUFFER_SZ / 128;
        bit_triple_a           = new bool[num_triples_pool];
        bit_triple_b           = new bool[num_triples_pool];
        bit_triple_c           = new bool[num_triples_pool];
        btp->get_triple(bit_triple_a, bit_triple_b, bit_triple_c);
        block_triple_a = new block[num_block_triples_pool];
        block_triple_b = new block[num_block_triples_pool];
        block_triple_c = new block[num_block_triples_pool];
        btp->get_block_triple(block_triple_a, block_triple_b, block_triple_c);
        this->io = io;
    }

    ~MPSIMDCircExec() {
        delete btp;
        delete[] bit_triple_a;
        delete[] bit_triple_b;
        delete[] bit_triple_c;

        delete[] block_triple_a;
        delete[] block_triple_b;
        delete[] block_triple_c;
        std::cout << "Total time in SIMD: " << total_time << " us\n";
    };

    void and_gate(bool* out1, bool* in1, bool* in2, size_t length) {
        auto t  = clock_start();
        bool *a = nullptr, *b = nullptr, *c = nullptr;
        bool delete_array = false;
        if (length > num_triples_pool) {
            a = new bool[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
            b = new bool[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
            c = new bool[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
            for (uint i = 0; i < (length + num_triples_pool - 1) / num_triples_pool; ++i)
                btp->get_triple(a + i * num_triples_pool, b + i * num_triples_pool, c + i * num_triples_pool);
            size_t tocp =
                min((length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length, num_triples);
            memcpy(bit_triple_a, a + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
                   tocp);
            memcpy(bit_triple_b, b + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
                   tocp);
            memcpy(bit_triple_c, c + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
                   tocp);
            num_triples  = 0;
            delete_array = true;
        }
        else if (length > num_triples_pool - num_triples) {  // buffer is not long enough
            a            = new bool[length];
            b            = new bool[length];
            c            = new bool[length];
            delete_array = true;
            memcpy(a, bit_triple_a + num_triples, num_triples_pool - num_triples);
            memcpy(b, bit_triple_b + num_triples, num_triples_pool - num_triples);
            memcpy(c, bit_triple_c + num_triples, num_triples_pool - num_triples);
            btp->get_triple(bit_triple_a, bit_triple_b, bit_triple_c);
            memcpy(a + num_triples_pool - num_triples, bit_triple_a, length - (num_triples_pool - num_triples));
            memcpy(b + num_triples_pool - num_triples, bit_triple_b, length - (num_triples_pool - num_triples));
            memcpy(c + num_triples_pool - num_triples, bit_triple_c, length - (num_triples_pool - num_triples));
            num_triples = length - (num_triples_pool - num_triples);
        }
        else {
            a = bit_triple_a + num_triples;
            b = bit_triple_b + num_triples;
            c = bit_triple_c + num_triples;
            num_triples += length;
        }

        bool *d = new bool[length], *e = new bool[length];

        for (uint i = 0; i < length; ++i) {
            d[i] = in1[i] ^ a[i];
            e[i] = in2[i] ^ b[i];
        }

        io->sync();
        if (cur_party == ALICE) {
            bool *d0 = new bool[length], *e0 = new bool[length];

            for (int i = 2; i <= num_party; ++i) {
                io->recv_bool(i, d0, length);
                io->recv_bool(i, e0, length);
                xorBools_arr(d, d, d0, length);
                xorBools_arr(e, e, e0, length);
            }

            vector<future<void>> res;
            for (int i = 2; i <= num_party; ++i) {
                res.push_back(pool->enqueue([this, i, d, length]() { this->io->send_bool(i, d, length); }));
                res.push_back(pool->enqueue([this, i, e, length]() { this->io->send_bool(i, e, length); }));
            }

            for (auto& v : res)
                v.get();
            res.clear();

            delete[] d0;
            delete[] e0;
        }
        else {
            io->send_bool(ALICE, d, length);
            io->send_bool(ALICE, e, length);

            io->recv_bool(ALICE, d, length);
            io->recv_bool(ALICE, e, length);
        }
        io->flush();
        if (cur_party == ALICE) {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i] ^ (d[i] & e[i]);
        }
        else {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i];
        }
        delete[] d;
        delete[] e;
        if (delete_array) {
            delete[] a;
            delete[] b;
            delete[] c;
        }
        total_time += time_from(t);
    }

    void and_gate(block* out1, block* in1, block* in2, size_t length) {
        auto t   = clock_start();
        block *a = nullptr, *b = nullptr, *c = nullptr;
        bool delete_array = false;
        if (length > num_block_triples_pool) {
            size_t alloc_length =
                (length + num_block_triples_pool - 1) / num_block_triples_pool * num_block_triples_pool;
            a = new block[alloc_length];
            b = new block[alloc_length];
            c = new block[alloc_length];
            for (uint i = 0; i < alloc_length / num_block_triples_pool; ++i) {
                btp->get_block_triple(a + i * num_block_triples_pool, b + i * num_block_triples_pool,
                                      c + i * num_block_triples_pool);
            }
            size_t tocp = min(alloc_length - length, num_block_triples);
            memcpy(block_triple_a, a + alloc_length - length, tocp);
            memcpy(block_triple_b, b + alloc_length - length, tocp);
            memcpy(block_triple_c, c + alloc_length - length, tocp);
            num_block_triples = 0;
            delete_array      = true;
        }
        else if (length > (num_block_triples_pool - num_block_triples)) {  // buffer is not long enough
            a            = new block[length];
            b            = new block[length];
            c            = new block[length];
            delete_array = true;
            memcpy(a, block_triple_a + num_block_triples, 16 * (num_block_triples_pool - num_block_triples));
            memcpy(b, block_triple_b + num_block_triples, 16 * (num_block_triples_pool - num_block_triples));
            memcpy(c, block_triple_c + num_block_triples, 16 * (num_block_triples_pool - num_block_triples));
            btp->get_block_triple(block_triple_a, block_triple_b, block_triple_c);
            memcpy(a + (num_block_triples_pool - num_block_triples), block_triple_a,
                   16 * (length - (num_block_triples_pool - num_block_triples)));
            memcpy(b + (num_block_triples_pool - num_block_triples), block_triple_b,
                   16 * (length - (num_block_triples_pool - num_block_triples)));
            memcpy(c + (num_block_triples_pool - num_block_triples), block_triple_c,
                   16 * (length - (num_block_triples_pool - num_block_triples)));
            num_block_triples = length - (num_block_triples_pool - num_block_triples);
        }
        else {
            a = block_triple_a + num_block_triples;
            b = block_triple_b + num_block_triples;
            c = block_triple_c + num_block_triples;
            num_block_triples += length;
        }

        if (length > d.size()) {
            d.resize(length);
            d1.resize(length);
            e.resize(length);
            e1.resize(length);
        }

        for (uint i = 0; i < length; ++i) {
            d[i] = in1[i] ^ a[i];
            e[i] = in2[i] ^ b[i];
        }

        io->sync();
        if (cur_party == ALICE) {
            for (int i = 2; i <= num_party; ++i) {
                io->recv_block(i, d1.data(), length);
                io->recv_block(i, e1.data(), length);
                xorBlocks_arr(d.data(), d.data(), d1.data(), length);
                xorBlocks_arr(e.data(), e.data(), e1.data(), length);
            }

            for (int i = 2; i <= num_party; ++i) {
                io->send_block(i, d.data(), length);
                io->send_block(i, e.data(), length);
            }
        }
        else {
            io->send_block(ALICE, d.data(), length);
            io->send_block(ALICE, e.data(), length);

            io->recv_block(ALICE, d.data(), length);
            io->recv_block(ALICE, e.data(), length);
        }
        io->flush();

        if (cur_party == ALICE) {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i] ^ (d[i] & e[i]);
        }
        else {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i];
        }

        if (delete_array) {
            delete[] a;
            delete[] b;
            delete[] c;
        }
        total_time += time_from(t);
    }

    void xor_gate(bool* out1, bool* in1, bool* in2, size_t length) {
        //auto t = clock_start();
        xorBools_arr(out1, in1, in2, length);
        //total_time += time_from(t);
    }

    void xor_gate(block* out1, block* in1, block* in2, size_t length) {
        //auto t = clock_start();
        for (uint i = 0; i < length; ++i) {
            out1[i] = in1[i] ^ in2[i];
        }
        //total_time += time_from(t);
    }

    void not_gate(bool* out1, bool* in1, size_t length) {
        //auto t = clock_start();
        bool bit_to_xor = (cur_party == ALICE);
        for (uint i = 0; i < length; ++i) {
            (out1[i]) = (in1[i]) ^ bit_to_xor;
        }
        //total_time += time_from(t);
    }

    void not_gate(block* out1, block* in1, size_t length) {
        //    auto t = clock_start();
        block bit_to_xor = (cur_party == ALICE) ? all_one_block : zero_block;
        for (uint i = 0; i < length; ++i) {
            out1[i] = in1[i] ^ bit_to_xor;
        }
        //  total_time += time_from(t);
    }
};

}  // namespace emp
