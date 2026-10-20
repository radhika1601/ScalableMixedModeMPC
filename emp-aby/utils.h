#pragma once
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include <cstdlib>
#include <future>
#include <string>
#include <vector>
namespace emp {
using std::future;
using std::min;
using std::string;
using std::vector;
using tuning::ferret_b13;

// Ferret b13 RCOTs produced per round, and the triple buffer size derived from it.
inline constexpr int64_t FERRET_N = ferret_b13.t << ferret_b13.tree_depth;
inline constexpr int64_t FERRET_BUFFER_SZ = ((FERRET_N - ferret_b13.M) / 128) * 128;

// argv[1] is the party id, argv[2] the port.
inline void parse_party_and_port(char** argv, int* party, int* port, int max_party = BOB) {
    *party = parse_party(argv, max_party);
    *port  = atoi(argv[2]);
}

inline void andBlocks_arr(block* res, const block* x, const block* y, int nblocks) {
    const block* dest = nblocks + x;
    for (; x != dest;)
        *(res++) = *(x++) & *(y++);
}

inline void andBlocks_arr(block* res, const block* x, block y, int nblocks) {
    const block* dest = nblocks + x;
    for (; x != dest;)
        *(res++) = *(x++) & y;
}

inline void xorBools_arr(bool* res, const bool* x, const bool* y, int nbools) {
    const bool* dest = nbools + x;
    for (; x != dest;)
        *(res++) = *(x++) ^ *(y++);
}

inline void andBools_arr(bool* res, const bool* x, const bool* y, int nbools) {
    const bool* dest = nbools + x;
    // Don't change & to &&
    for (; x != dest;)
        *(res++) = *(x++) & *(y++);
}

inline void bool_to_block_arr(block* a, bool* data, int bool_length) {
    uint64_t* data64 = (uint64_t*)data;
    uint64_t low = 0, high = 0;
    int i = 0;
    for (; i < bool_length / 8; ++i) {
        unsigned long long mask = 0x0101010101010101ULL;
        unsigned long long tmp  = 0;
#if defined(__BMI2__)
        tmp = _pext_u64(data64[i], mask);
#else
        // https://github.com/Forceflow/libmorton/issues/6
        for (unsigned long long bb = 1; mask != 0; bb += bb) {
            if (data64[i] & mask & -mask) {
                tmp |= bb;
            }
            mask &= (mask - 1);
        }
#endif
        if (i > 0 && i % 16 == 0) {
            a[i / 16 - 1] = makeBlock(high, low);
            high          = 0;
            low           = tmp;
        }
        else if (i % 16 < 8) {
            low += tmp << (i % 16) * 8;
        }
        else if (i % 16 >= 8) {
            high += tmp << ((i % 16 - 8) * 8);
        }
    }
    if (i % 16 == 0) {
        a[i / 16 - 1] = makeBlock(high, low);
    }
    else {
        a[i / 16] = makeBlock(high, low);
    }
}
}  // namespace emp
