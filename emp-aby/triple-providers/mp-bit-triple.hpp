template <typename IO>
MPBitTripleProvider<IO>::MPBitTripleProvider(int num_party, int party, ThreadPool* pool, MPIOChannel<IO>* io,
                                             const int buffer_length) {
    this->party   = party;
    this->threads = pool->size();
    if (this->threads % 2)
        error("MPBitTripleProvider needs even number of threads!");
    this->io        = io;
    this->num_party = num_party;
    this->BUFFER_SZ = buffer_length;
    this->pool      = pool;
    a_bool.resize(BUFFER_SZ);
    b_bool.resize(BUFFER_SZ);
    c_bool.resize(BUFFER_SZ);
    mac      = (block**)(malloc((threads / 2) * sizeof(block*)));
    key      = (block**)(malloc((threads / 2) * sizeof(block*)));
    key_star = (block**)(malloc((threads / 2) * sizeof(block*)));
    xor_mac  = (block**)(malloc((threads / 2) * sizeof(block*)));
    xor_key  = (block**)(malloc((threads / 2) * sizeof(block*)));
    for (int i = 0; i < threads / 2; ++i) {
        mac[i]      = new block[BUFFER_SZ];
        key[i]      = new block[BUFFER_SZ];
        xor_mac[i]  = new block[BUFFER_SZ];
        xor_key[i]  = new block[BUFFER_SZ];
        key_star[i] = new block[BUFFER_SZ];
    }
    delta = gen_delta();
    bool delta_bits[128];
    bits_to_bools(delta_bits, &delta, 128);
    // All receivers of this party must output the same choice bits, so they share one choice seed.
    block choice_seed;
    prg.random_block(&choice_seed, 1);
    // ALICE broadcasts a public LPN matrix seed shared by every Ferret instance.
    seed_gen();
    cot_sender.assign(num_party, nullptr);
    cot_receiver.assign(num_party, nullptr);
    // Create the semi-honest Ferret instances. Bootstrap runs lazily on the first rcot.
    for (int i = 0; i < num_party; ++i) {
        int peer = i + 1;
        if (party == peer)
            continue;
        cot_sender[i] = new Ferret(ALICE, io->get(peer, party > peer), false);
        cot_sender[i]->set_sid(makeBlock(party, peer));
        cot_sender[i]->set_lpn_matrix_seed(seed);
        cot_sender[i]->set_delta(delta_bits);

        cot_receiver[i] = new Ferret(BOB, io->get(peer, party < peer), false);
        cot_receiver[i]->set_sid(makeBlock(peer, party));
        cot_receiver[i]->set_lpn_matrix_seed(seed);
        cot_receiver[i]->set_choice_seed(choice_seed);
    }
}

template <typename IO>
MPBitTripleProvider<IO>::~MPBitTripleProvider() {
    for (int i = 0; i < threads / 2; ++i) {
        delete[] mac[i];
        delete[] key[i];
        delete[] key_star[i];
        delete[] xor_mac[i];
        delete[] xor_key[i];
    }
    free(mac);
    free(key);
    free(key_star);
    free(xor_mac);
    free(xor_key);
    for (int i = 0; i < num_party; ++i) {
        delete cot_sender[i];
        delete cot_receiver[i];
    }
}

template <typename IO>
void MPBitTripleProvider<IO>::get_triple(block* a, block* b, block* c) {
    this->get_triple((bool*)a_bool.data(), (bool*)b_bool.data(), (bool*)c_bool.data());
    bool_to_block_arr(a, (bool*)a_bool.data(), BUFFER_SZ);
    bool_to_block_arr(b, (bool*)b_bool.data(), BUFFER_SZ);
    bool_to_block_arr(c, (bool*)c_bool.data(), BUFFER_SZ);
}

template <typename IO>
void MPBitTripleProvider<IO>::get_triple(bool* a, bool* b, bool* c) {
    prg.random_bool(a, BUFFER_SZ);
    ch[0] = zero_block;
    ch[1] = all_one_block;
    block **w, **s;
    w = (block**)malloc((threads / 2) * sizeof(block*));
    s = (block**)malloc(threads * sizeof(block*));
    for (int i = 0; i < threads / 2; ++i) {
        w[i] = new block[BUFFER_SZ];
        memset(w[i], 0, BUFFER_SZ * sizeof(block));
        memset(xor_mac[i], 0, BUFFER_SZ * sizeof(block));
        memset(xor_key[i], 0, BUFFER_SZ * sizeof(block));
    }
    for (int i = 0; i < threads; ++i) {
        s[i] = new block[BUFFER_SZ];
    }
    // Peer offsets 1 .. num_party-1 are split into contiguous runs, one per worker.
    // Each party sends to party+offset and receives from party-offset, so the two ends
    // of every pair handle it in the same worker at the same step.
    int num_steps = ceil((double)(num_party - 1) / ((double)threads / 2));
    vector<future<void>> res;
    for (int i = 0; i < threads / 2; ++i) {
        res.push_back(pool->enqueue([this, i, num_steps, a, s] {
            for (int offset = num_steps * i + 1; offset <= min(num_steps * (i + 1), num_party - 1); ++offset) {
                int send_to = ((party - 1) + offset) % num_party;
                send(send_to, a, s[i], i);
            }
        }));
    }
    for (int i = 0; i < threads / 2; ++i) {
        res.push_back(pool->enqueue([this, i, num_steps, b, w, s] {
            for (int offset = num_steps * i + 1; offset <= min(num_steps * (i + 1), num_party - 1); ++offset) {
                int receive_from = ((party - 1) + num_party - offset) % num_party;
                recv(receive_from, offset == 1 ? b : nullptr, w[i], s[i + threads / 2], i);
            }
        }));
    }

    for (auto& v : res)
        v.get();
    res.clear();

    for (int i = 1; i < threads / 2; ++i) {
        xorBlocks_arr(w[0], w[0], w[i], BUFFER_SZ);
    }

    for (int i = 0; i < BUFFER_SZ; ++i) {
        w[0][i] = w[0][i] & ch[b[i]];
    }

    for (int i = 0; i < threads / 2; ++i) {
        xorBlocks_arr(w[0], xor_mac[i], w[0], BUFFER_SZ);
        xorBlocks_arr(w[0], xor_key[i], w[0], BUFFER_SZ);
    }
    for (int i = 0; i < BUFFER_SZ; ++i) {
        c[i] = (a[i] & b[i]) ^ getLSB(w[0][i]);
    }

    for (int i = 0; i < threads / 2; ++i)
        delete[] w[i];
    for (int i = 0; i < threads; ++i)
        delete[] s[i];
    free(s);
    free(w);
}
