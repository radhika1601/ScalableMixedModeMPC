#ifndef EMP__SIMD_CIRCUIT_EXECUTION_H__
#define EMP__SIMD_CIRCUIT_EXECUTION_H__
#include "emp-aby/triple-providers/bit-triple.h"

namespace emp {

class SIMDCircuitExecution {
public:
#ifndef THREADING
    static SIMDCircuitExecution* simd_circ_exec;
#else
    static __thread SIMDCircuitExecution* simd_circ_exec;
#endif
    virtual void and_gate(block* out1, block* in1, block* in2, size_t length) = 0;
    virtual void xor_gate(block* out1, block* in1, block* in2, size_t length) = 0;
    virtual void not_gate(block* out1, block* in1, size_t length)             = 0;
    // virtual void mux_gate(block** out1, block** in1, block** in2, bool** select, int width, size_t length);
    virtual void and_gate(bool* out1, bool* in1, bool* in2, size_t length) = 0;
    virtual void xor_gate(bool* out1, bool* in1, bool* in2, size_t length) = 0;
    virtual void not_gate(bool* out1, bool* in1, size_t length)            = 0;

    virtual ~SIMDCircuitExecution() {}

private:
};

}  // namespace emp

#endif  // EMP__SIMD_CIRCUIT_EXECUTION_H__