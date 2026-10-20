#pragma once
#include <cstdint>

// Triple buffer size for tests. Must be a multiple of 128.
// The library default (FERRET_BUFFER_SZ, about 15M) needs several GB per party,
// which does not fit when all parties run on one machine.
inline constexpr int64_t TEST_BUFFER_SZ = int64_t(1) << 20;
