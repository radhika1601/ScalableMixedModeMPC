#ifndef EMP__SIMD_CIRCUIT_EXECUTION_H__
#define EMP__SIMD_CIRCUIT_EXECUTION_H__
#include "emp-aby/triple-providers/bit-triple.h"

namespace emp {

template <typename BTP>
class SIMDCircuitExecution {
public:
#ifndef THREADING
    static SIMDCircuitExecution* simd_circ_exec;
#else
    static __thread SIMDCircuitExecution* simd_circ_exec;
#endif
    virtual BTP* getBtp() = 0;
    virtual void and_gate(block* out1, block* in1, block* in2, size_t length) = 0;
    virtual void xor_gate(block* out1, block* in1, block* in2, size_t length) = 0;
    virtual void not_gate(block* out1, block* in1, size_t length)             = 0;
    virtual void and_gate(bool* out, bool* in1, bool* in2, size_t bool_length, block* block_out, block* block_in1,
                          block* block_in2, size_t length)                    = 0;
    // virtual void mux_gate(block** out1, block** in1, block** in2, bool** select, int width, size_t length);
    virtual void and_gate(bool* out1, bool* in1, bool* in2, size_t length) = 0;
    virtual void xor_gate(bool* out1, bool* in1, bool* in2, size_t length) = 0;
    virtual void not_gate(bool* out1, bool* in1, size_t length)            = 0;

    virtual ~SIMDCircuitExecution() {}

    protected:
    
    template <typename T>
    inline void and_helper(T*& a, T*& b, T*& c, size_t length, bool& delete_array, T* bit_triple_a, T* bit_triple_b,
                           T* bit_triple_c, size_t num_triples_pool, size_t& num_triples) {
        BTP* btp = getBtp();
        if (length > num_triples_pool) {
            size_t total = (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool;
            a            = new T[total];
            b            = new T[total];
            c            = new T[total];
            for (uint i = 0; i < total / num_triples_pool; ++i)
                btp->get_triple(a + i * num_triples_pool, b + i * num_triples_pool, c + i * num_triples_pool);
            // The unused triples are a[length, total). Move as many as fit into the consumed
            // slots just below num_triples, so the pool stays a contiguous unused suffix.
            size_t tocp = min(total - length, num_triples);
            memcpy(bit_triple_a + num_triples - tocp, a + length, sizeof(T) * tocp);
            memcpy(bit_triple_b + num_triples - tocp, b + length, sizeof(T) * tocp);
            memcpy(bit_triple_c + num_triples - tocp, c + length, sizeof(T) * tocp);
            num_triples -= tocp;
            delete_array = true;
        }
        else if (length > num_triples_pool - num_triples) {  // buffer is not long enough
            a            = new T[length];
            b            = new T[length];
            c            = new T[length];
            delete_array = true;
            memcpy(a, bit_triple_a + num_triples, sizeof(T) * (num_triples_pool - num_triples));
            memcpy(b, bit_triple_b + num_triples, sizeof(T) * (num_triples_pool - num_triples));
            memcpy(c, bit_triple_c + num_triples, sizeof(T) * (num_triples_pool - num_triples));
            btp->get_triple(bit_triple_a, bit_triple_b, bit_triple_c);
            memcpy(a + num_triples_pool - num_triples, bit_triple_a,
                   sizeof(T) * (length - (num_triples_pool - num_triples)));
            memcpy(b + num_triples_pool - num_triples, bit_triple_b,
                   sizeof(T) * (length - (num_triples_pool - num_triples)));
            memcpy(c + num_triples_pool - num_triples, bit_triple_c,
                   sizeof(T) * (length - (num_triples_pool - num_triples)));
            num_triples = length - (num_triples_pool - num_triples);
        }
        else {
            a = bit_triple_a + num_triples;
            b = bit_triple_b + num_triples;
            c = bit_triple_c + num_triples;
            num_triples += length;
        }
    }

private:
};

}  // namespace emp

#endif  // EMP__SIMD_CIRCUIT_EXECUTION_H__