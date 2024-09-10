#ifndef BIT_TRIPLE_H__
#define BIT_TRIPLE_H__
#include "emp-ot/emp-ot.h"
#include "emp-aby/utils.h"

namespace emp {

template <typename IO>
class BitTripleProvider {
public:
    BitTripleProvider(int party, int threads, IO** ios);

    void get_triple(block* a, block* b, block* c);  //length should be the #blocks
    void get_triple(bool* a, bool* b, bool* c);     // length ... # bools

    ~BitTripleProvider();

    static long long int BUFFER_SZ;

private:
    int party, threads;
    IO** ios;
    FerretCOT<IO>*ot0 = nullptr, *ot1 = nullptr;
    CRH crh;
    PRG prg;
    block delta;
    vector<block> r0, r1, scratch, A_hat, A_star, B_hat;
    vector<unsigned char> a_bool, b_bool, c_bool;
    void compute_rcots(bool* b);
};

#include "bit-triple.hpp"
}  // namespace emp

#endif  //BIT_TRIPLE_H__
