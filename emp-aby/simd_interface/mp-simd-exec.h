#pragma once

#include "simd.h"
#include "emp-aby/triple-providers/mp-bit-triple.h"

namespace emp {

template <typename IO>
class MPSIMDCircExec : SIMDCircuitExecution<MPBitTripleProvider<IO>> {
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
        this->cur_party        = party;
        this->num_party        = num_party;
        this->pool             = pool;
        btp                    = new MPBitTripleProvider<IO>(num_party, party, pool, io);
        num_triples_pool       = btp->BUFFER_SZ;
        num_block_triples_pool = btp->BUFFER_SZ / 128;
        bit_triple_a           = new bool[num_triples_pool];
        bit_triple_b           = new bool[num_triples_pool];
        bit_triple_c           = new bool[num_triples_pool];
        btp->get_triple(bit_triple_a, bit_triple_b, bit_triple_c);
        block_triple_a = new block[num_block_triples_pool];
        block_triple_b = new block[num_block_triples_pool];
        block_triple_c = new block[num_block_triples_pool];
        btp->get_triple(block_triple_a, block_triple_b, block_triple_c);
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

    MPBitTripleProvider<IO> * getBtp(){
        return btp;
    }

    void and_gate(bool* out1, bool* in1, bool* in2, size_t length) {
        auto t  = clock_start();
        bool *a = nullptr, *b = nullptr, *c = nullptr;
        bool delete_array = false;
        this->template and_helper<bool>(a, b, c, length, delete_array, bit_triple_a, bit_triple_b, bit_triple_c, num_triples_pool,
                         num_triples);
        bool *d = new bool[length], *e = new bool[length];

        for (uint i = 0; i < length; ++i) {
            d[i] = in1[i] ^ a[i];
            e[i] = in2[i] ^ b[i];
        }

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
                res.push_back(pool->enqueue([this, i, d, e, length]() {
                    this->io->send_bool(i, d, length);
                    this->io->send_bool(i, e, length);
                }));
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
        this->template and_helper<block>(a, b, c, length, delete_array, block_triple_a, block_triple_b, block_triple_c,
                          num_block_triples_pool, num_block_triples);
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

    void and_gate(bool* out, bool* in1, bool* in2, size_t bool_length, block* block_out, block* block_in1,
                  block* block_in2, size_t length) {
        auto t  = clock_start();
        bool *a = nullptr, *b = nullptr, *c = nullptr;
        bool delete_array = false;
        this->template and_helper<bool>(a, b, c, bool_length, delete_array, bit_triple_a, bit_triple_b, bit_triple_c, num_triples_pool,
                         num_triples);
        bool *d = new bool[bool_length], *e = new bool[bool_length];

        for (uint i = 0; i < bool_length; ++i) {
            d[i] = in1[i] ^ a[i];
            e[i] = in2[i] ^ b[i];
        }

        block *block_a = nullptr, *block_b = nullptr, *block_c = nullptr;
        bool delete_block_array = false;
        this->template and_helper<block>(block_a, block_b, block_c, length, delete_block_array, block_triple_a, block_triple_b,
                          block_triple_c, num_block_triples_pool, num_block_triples);

        if (length > this->d.size()) {
            this->d.resize(length);
            this->d1.resize(length);
            this->e.resize(length);
            this->e1.resize(length);
        }

        for (uint i = 0; i < length; ++i) {
            this->d[i] = block_in1[i] ^ block_a[i];
            this->e[i] = block_in2[i] ^ block_b[i];
        }

        if (cur_party == ALICE) {
            bool *d0 = new bool[bool_length], *e0 = new bool[bool_length];

            for (int i = 2; i <= num_party; ++i) {
                io->recv_bool(i, d0, bool_length);
                io->recv_bool(i, e0, bool_length);
                xorBools_arr(d, d, d0, bool_length);
                xorBools_arr(e, e, e0, bool_length);
                io->recv_block(i, this->d1.data(), length);
                io->recv_block(i, this->e1.data(), length);
                xorBlocks_arr(this->d.data(), this->d.data(), this->d1.data(), length);
                xorBlocks_arr(this->e.data(), this->e.data(), this->e1.data(), length);
            }

            vector<future<void>> res;
            for (int i = 2; i <= num_party; ++i) {
                res.push_back(pool->enqueue([this, i, d, e, bool_length, length]() {
                    this->io->send_bool(i, d, bool_length);
                    this->io->send_bool(i, e, bool_length);
                    io->send_block(i, this->d.data(), length);
                    io->send_block(i, this->e.data(), length);
                }));
            }

            for (auto& v : res)
                v.get();
            res.clear();

            delete[] d0;
            delete[] e0;
        }
        else {
            io->send_bool(ALICE, d, bool_length);
            io->send_bool(ALICE, e, bool_length);
            io->send_block(ALICE, this->d.data(), length);
            io->send_block(ALICE, this->e.data(), length);

            io->recv_bool(ALICE, d, bool_length);
            io->recv_bool(ALICE, e, bool_length);
            io->recv_block(ALICE, this->d.data(), length);
            io->recv_block(ALICE, this->e.data(), length);
        }
        io->flush();
        if (cur_party == ALICE) {
            for (uint i = 0; i < bool_length; ++i)
                out[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i] ^ (d[i] & e[i]);
            for (uint i = 0; i < length; ++i)
                block_out[i] =
                    (this->d[i] & block_b[i]) ^ (this->e[i] & block_a[i]) ^ block_c[i] ^ (this->d[i] & this->e[i]);
        }
        else {
            for (uint i = 0; i < bool_length; ++i)
                out[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i];
            for (uint i = 0; i < length; ++i)
                block_out[i] = (this->d[i] & block_b[i]) ^ (this->e[i] & block_a[i]) ^ block_c[i];
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
