#ifndef EMP_SIMD_EXEC_H__
#define EMP_SIMD_EXEC_H__
#include "emp-aby/simd_interface/simd.h"

namespace emp {

template <typename IO>
class SIMDCircExec: SIMDCircuitExecution {
private:
    BitTripleProvider<IO>* btp;
    IO* io;
    bool *bit_triple_a, *bit_triple_b, *bit_triple_c, *sel_mux;
    block *block_triple_a, *block_triple_b, *block_triple_c;
    block *A_mux, *B_mux;
    int length_mux, width_mux;
    size_t num_triples = 0, num_block_triples = 0;
    size_t num_triples_pool;
    size_t num_block_triples_pool;

    template <typename T>
    void trim_arr(T*& a, int prev_length, int length_cut) {
        T* chunk = new T[prev_length - length_cut];
        memcpy(chunk, a + length_cut, (prev_length - length_cut) * sizeof(T));
        delete[] a;
        a = chunk;
    }
    double total_time = 0;
    vector<block> d, d1, e, e1;

public:
    int cur_party;
    size_t num_simd_and_gates = 0, num_and_gates = 0;

    SIMDCircExec(int party, int threads, IO** ios, int length_mux = 0, int width_mux = 0) {
        this->cur_party        = party;
        btp                    = new BitTripleProvider<IO>(party, threads, ios);
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
        /*/    if (length_mux > 0 && width_mux > 0) {
      A_mux = new block[length_mux * width_mux];
      B_mux = new block[length_mux * width_mux];
      sel_mux = new bool[length_mux];
      btp->get_mux_triple(A_mux, sel_mux, B_mux, width_mux, length_mux);
      }
      this->length_mux = length_mux;
      this->width_mux = width_mux;*/
        io = ios[0];
    }

    ~SIMDCircExec() {
        delete btp;
        delete[] bit_triple_a;
        delete[] bit_triple_b;
        delete[] bit_triple_c;

        delete[] block_triple_a;
        delete[] block_triple_b;
        delete[] block_triple_c;
        /*    delete[] A_mux;
                    delete[] B_mux;
                    delete[] sel_mux;*/
        std::cout << "Total time in SIMD: " << total_time << " us\n";
    };
    void and_gate(bool* out1, bool* in1, bool* in2, size_t length) {
        num_simd_and_gates += 1;
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

        bool *d = new bool[length], *d1 = new bool[length];
        bool *e = new bool[length], *e1 = new bool[length];

        for (uint i = 0; i < length; ++i) {
            d[i] = in1[i] ^ a[i];
            e[i] = in2[i] ^ b[i];
        }
        if (cur_party == ALICE) {
            io->send_bool(d, length);
            io->send_bool(e, length);

            io->recv_bool(d1, length);
            io->recv_bool(e1, length);
        }
        else if (cur_party == BOB) {
            io->recv_bool(d1, length);
            io->recv_bool(e1, length);

            io->send_bool(d, length);
            io->send_bool(e, length);
            io->flush();
        }
        xorBools_arr(d, d, d1, length);
        xorBools_arr(e, e, e1, length);

        if (cur_party == ALICE) {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i] ^ (d[i] & e[i]);
        }
        else if (cur_party == BOB) {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i];
        }
        delete[] d;
        delete[] e;
        delete[] d1;
        delete[] e1;
        if (delete_array) {
            delete[] a;
            delete[] b;
            delete[] c;
        }
        total_time += time_from(t);
    }

    void and_gate(block* out1, block* in1, block* in2, size_t length) {
        num_simd_and_gates += 1;
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

        //    block *d = new block[length], *d1 = new block[length];
        //    block *e = new block[length], *e1 = new block[length];
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
            io->send_block(d.data(), length);
            io->send_block(e.data(), length);

            io->recv_block(d1.data(), length);
            io->recv_block(e1.data(), length);
        }
        else if (cur_party == BOB) {
            io->recv_block(d1.data(), length);
            io->recv_block(e1.data(), length);

            io->send_block(d.data(), length);
            io->send_block(e.data(), length);
            io->flush();
        }
        xorBlocks_arr(d.data(), d.data(), d1.data(), length);
        xorBlocks_arr(e.data(), e.data(), e1.data(), length);

        if (cur_party == ALICE) {
            for (uint i = 0; i < length; ++i)
                out1[i] = (d[i] & b[i]) ^ (e[i] & a[i]) ^ c[i] ^ (d[i] & e[i]);
        }
        else if (cur_party == BOB) {
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
        bool bit_to_xor = (cur_party == ALICE) ? true : false;
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

    void mux_gate(block** out1, block** in1, block** in2, bool** select, int width, size_t length) {
        if (width_mux < width) {
            this->length_mux = length;
            this->width_mux  = width;
            delete[] A_mux;
            delete[] B_mux;
            delete[] sel_mux;
            A_mux   = new block[length_mux * width_mux];
            B_mux   = new block[length_mux * width_mux];
            sel_mux = new bool[length_mux];
            btp->get_mux_triple(A_mux, sel_mux, B_mux, width_mux, length_mux);
        }
        if (length_mux < length) {
            block* x = new block[(length + length_mux) * width_mux];
            block* y = new block[(length + length_mux) * width_mux];
            bool* z  = new bool[length + length_mux];
            btp->get_mux_triple(x, z, y, width_mux, length);
            memcpy(x + length * width_mux, A_mux, length_mux * width_mux);
            memcpy(y + length * width_mux, B_mux, length_mux * width_mux);
            memcpy(z + length, sel_mux, length_mux);
            delete[] A_mux;
            delete[] B_mux;
            delete[] sel_mux;
            A_mux   = x;
            B_mux   = y;
            sel_mux = z;
            length_mux += length;
        }
        bool* d  = new bool[length];
        block *v = new block[length * width], *e = new block[length * width];

        bool* d1  = new bool[length];
        block* e1 = new block[length * width];

        for (int i = 0; i < length * width; ++i) {
            v[i] = *(in1[i]) ^ *(in2[i]);
        }
        for (uint i = 0; i < length; ++i) {
            d[i] = *(select[i]) ^ sel_mux[i];
        }

        for (uint i = 0; i < length; ++i) {
            xorBlocks_arr(e + i * width, v + i * width, A_mux + i * width_mux, width);
        }
        if (cur_party == ALICE) {
            io->send_bool(d, length);
            io->send_block(e, length * width);
            io->recv_bool(d1, length);
            io->recv_block(e1, length * width);
        }
        else if (cur_party == BOB) {
            io->recv_bool(d1, length);
            io->recv_block(e1, length * width);
            io->send_bool(d, length);
            io->send_block(e, length * width);
            io->flush();
        }

        xorBlocks_arr(e, e, e1, length * width);
        xorBools_arr(d, d, d1, length);

        for (uint i = 0; i < length; ++i) {
            for (int j = 0; j < width; ++j) {
                *(out1[i * width + j]) = B_mux[i * width_mux + j] ^ *(in1[i * width + j]);
                if (d[i]) {
                    *(out1[i * width + j]) ^= e[i * width + j] ^ A_mux[i * width_mux + j];
                }
                if (select[i]) {
                    *(out1[i * width + j]) ^= e[i * width + j];
                }
            }
        }

        trim_arr(A_mux, length_mux * width_mux, length * width_mux);
        trim_arr(B_mux, length_mux * width_mux, length * width_mux);
        trim_arr(sel_mux, length_mux, length);
        length_mux -= length;

        delete[] d;
        delete[] e;
        delete[] d1;
        delete[] e1;
        delete[] v;
    }
};

}  // namespace emp

#endif  // EMP_SIMD_EXEC_H__