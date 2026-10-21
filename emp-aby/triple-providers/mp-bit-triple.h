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
        int buffer_length = FERRET_BUFFER_SZ);

    void get_triple(block* a, block* b, block* c);  //length should be the #blocks
    void get_triple(bool* a, bool* b, bool* c);     // length ... # bools
    ~MPBitTripleProvider();

    int BUFFER_SZ;

private:
    vector<Ferret*> cot_sender;
    vector<Ferret*> cot_receiver;
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
        delta_i = delta_i | makeBlock(0, 1);

        return delta_i;
    }

    block seed;
    int party, threads;
    ThreadPool* pool    = nullptr;
    MPIOChannel<IO>* io = nullptr;
    PRG prg;
    vector<unsigned char> a_bool, b_bool, c_bool;
    block **key_star, **xor_mac, **xor_key, **key, **mac;
    CCRH crh;
    block delta;
    int num_party;
    void send(int send_to, const bool* a, block* s, int thread_idx) {
        cot_sender[send_to]->rcot(this->key[thread_idx], BUFFER_SZ);
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

    // All receivers share one choice seed, so any one of them can write b. Pass b to
    // exactly one call and nullptr to the rest.
    void recv(int receive_from, bool* b, block* w, block* s, int thread_idx) {
        cot_receiver[receive_from]->rcot(this->mac[thread_idx], BUFFER_SZ);
        if (b != nullptr) {
            for (int i = 0; i < BUFFER_SZ; ++i)
                b[i] = getLSB(mac[thread_idx][i]);
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
