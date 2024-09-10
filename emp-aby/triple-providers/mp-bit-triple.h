#pragma once

#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "emp-aby/io/multi-io.hpp"
#include "emp-aby/utils.h"

namespace emp {
template <typename IO>
class MPBitTripleProvider {
public:
    MPBitTripleProvider(
        int num_party, int party, ThreadPool* pool, MPIOChannel<IO>* io,
        int buffer_length = ((ferret_b13.n - ferret_b13.k - ferret_b13.t * ferret_b13.log_bin_sz_pre - 128) / 128) *
                            128);

    void get_triple(block* a, block* b, block* c);  //length should be the #blocks
    void get_triple(bool* a, bool* b, bool* c);     // length ... # bools
    ~MPBitTripleProvider();

    int BUFFER_SZ;

private:
    vector<FerretCOT<IO>*> cot_sender;
    vector<FerretCOT<IO>*> cot_receiver;
    block ch[2];

    void seed_gen() {
        vector<future<void>> res;
        prg.random_block(&seed, 1);
        this->io->flush();
        if (party == ALICE) {
            for (int i = 2; i <= num_party; ++i) {
                res.push_back(pool->enqueue([this, i] { this->io->send_block(i, &seed, 1); }));
            }
        }
        else {
            io->recv_block(ALICE, &seed, 1);
        }
        for (auto& v : res)
            v.get();
        res.clear();
        this->io->flush();
    }

    block gen_delta() {
        block delta_i;

        prg.random_block(&delta_i, 1);
        block one = makeBlock(0xFFFFFFFFFFFFFFFFLL, 0xFFFFFFFFFFFFFFFELL);
        delta_i   = delta_i & one;
        delta_i   = delta_i ^ 0x1;

        return delta_i;
    }

    block seed;
    int party, threads;
    ThreadPool* pool    = nullptr;
    MPIOChannel<IO>* io = nullptr;
    PRG prg;
    vector<unsigned char> a_bool, b_bool, c_bool;
    block **key_star, **xor_mac, **xor_key, **key, **mac;
    bool *sent, *received;
    CRH crh;
    block delta;
    int num_party;
    bool b_set = false;
    void send(int send_to, const bool* a, block* s, int thread_idx) {
        if (send_to == party - 1)
            return;
        if (sent[send_to])
            return;
        sent[send_to] = true;
        cot_sender[send_to]->rcot_inplace(this->key[thread_idx], ferret_b13.n, seed);
        for (int i = 0; i < BUFFER_SZ; i += 128) {
            xorBlocks_arr(key_star[thread_idx] + i, key[thread_idx] + i, delta, 128);
            crh.Hn(key[thread_idx] + i, key[thread_idx] + i, 128);
            crh.Hn(key_star[thread_idx] + i, key_star[thread_idx] + i, 128);
        }
        io->get(send_to + 1, party > (send_to + 1))->flush();
        for (int j = 0; j < BUFFER_SZ; ++j) {
            s[j] = key[thread_idx][j] ^ key_star[thread_idx][j] ^ ch[a[j]];
        }
        io->get(send_to + 1, party > (send_to + 1))->send_block(s, BUFFER_SZ);
        io->get(send_to + 1, party > (send_to + 1))->flush();
        xorBlocks_arr(xor_key[thread_idx], xor_key[thread_idx], key[thread_idx], BUFFER_SZ);
    }

    void recv(int receive_from, bool* b, block* w, block* s, int thread_idx) {
        if (receive_from == party - 1)
            return;
        if (received[receive_from])
            return;
        received[receive_from] = true;
        cot_receiver[receive_from]->rcot_inplace(this->mac[thread_idx], ferret_b13.n, seed);
        if (!b_set) {
            for (int i = 0; i < BUFFER_SZ; ++i)
                b[i] = getLSB(mac[thread_idx][i]);

            b_set = true;
        }
        io->get(receive_from + 1, party < (receive_from + 1))->flush();
        io->get(receive_from + 1, party < (receive_from + 1))->recv_block(s, BUFFER_SZ);
        io->get(receive_from + 1, party < (receive_from + 1))->flush();

        for (int i = 0; i < BUFFER_SZ; i += 128) {
            crh.Hn(mac[thread_idx] + i, mac[thread_idx] + i, 128);
            xorBlocks_arr(w + i, s + i, w + i, 128);
            xorBlocks_arr(xor_mac[thread_idx] + i, xor_mac[thread_idx] + i, mac[thread_idx] + i, 128);
        }
    }
};

#include "emp-aby/triple-providers/mp-bit-triple.hpp"

}  // namespace emp
