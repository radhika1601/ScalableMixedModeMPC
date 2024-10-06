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
            a = new T[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
            b = new T[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
            c = new T[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
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