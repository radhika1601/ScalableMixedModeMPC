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
    sent     = new bool[num_party];
    received = new bool[num_party];
    mac      = (block**)(malloc((threads / 2) * sizeof(block*)));
    key      = (block**)(malloc((threads / 2) * sizeof(block*)));
    key_star = (block**)(malloc((threads / 2) * sizeof(block*)));
    xor_mac  = (block**)(malloc((threads / 2) * sizeof(block*)));
    xor_key  = (block**)(malloc((threads / 2) * sizeof(block*)));
    for (int i = 0; i < threads / 2; ++i) {
        mac[i]      = new block[emp::ferret_b13.n];
        key[i]      = new block[emp::ferret_b13.n];
        xor_mac[i]  = new block[BUFFER_SZ];
        xor_key[i]  = new block[BUFFER_SZ];
        key_star[i] = new block[BUFFER_SZ];
    }
    delta = gen_delta();
    cot_sender.resize((num_party) * sizeof(FerretCOT<IO>*));
    cot_receiver.resize((num_party) * sizeof(FerretCOT<IO>*));
    // Create the semi-honest Ferret instances.
    for (int i = 0; i < num_party; ++i) {
        if (party != i + 1) {
            cot_sender[i]   = new FerretCOT<IO>(ALICE, 1, &(io->get(i + 1, party > (i + 1))), false, false);
            cot_receiver[i] = new FerretCOT<IO>(BOB, 1, &(io->get(i + 1, party < (i + 1))), false, false);
        }
    }
    seed_gen();
    int l            = ferret_b13.log_bin_sz * ferret_b13.t + ferret_b13.k + 128;
    bool* pre_choice = new bool[l];
    prg.random_bool(pre_choice, l);
    vector<future<void>> res;
    for (int i = 0; i < num_party; ++i)
        if (party < (i + 1)) {
            res.push_back(pool->enqueue([this, io, i, pre_choice]() {
                cot_sender[i]->setup(this->delta,
                                     "./data/" + std::to_string(this->party) + "to" + std::to_string(i + 1) + ".txt",
                                     pre_choice, seed);
                io->flush(i + 1);
            }));
            res.push_back(pool->enqueue([this, io, i, pre_choice]() {
                cot_receiver[i]->setup(
                    "./data/" + std::to_string(this->party) + "from" + std::to_string(i + 1) + ".txt", pre_choice,
                    seed);
                io->flush(i + 1);
            }));
        }
        else if (party > (i + 1)) {
            res.push_back(pool->enqueue([this, io, i, pre_choice]() {
                cot_receiver[i]->setup(
                    "./data/" + std::to_string(this->party) + "from" + std::to_string(i + 1) + ".txt", pre_choice,
                    seed);
                io->flush(i + 1);
            }));
            res.push_back(pool->enqueue([this, io, i, pre_choice]() {
                cot_sender[i]->setup(this->delta,
                                     "./data/" + std::to_string(this->party) + "to" + std::to_string(i + 1) + ".txt",
                                     pre_choice, seed);
                io->flush(i + 1);
            }));
        }

    for (auto& v : res)
        v.get();
    res.clear();

    delete[] pre_choice;
}

template <typename IO>
MPBitTripleProvider<IO>::~MPBitTripleProvider() {
    delete io;
    delete pool;
    delete[] mac;
    delete[] key;
    delete[] key_star;
    delete[] xor_mac;
    delete[] xor_key;
    delete[] sent;
    delete[] received;
    for (int i = 0; i < num_party; ++i) {
        if (i + 1 == party)
            continue;
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
    memset(sent, 0, num_party * sizeof(bool));
    memset(received, 0, num_party * sizeof(bool));
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
    seed_gen();
    int num_steps = ceil((double)(num_party - 1) / ((double)threads / 2));
    vector<future<void>> res;
    for (int i = 0; i < threads / 2; ++i) {
        res.push_back(pool->enqueue([this, i, num_steps, a, s] {
            for (int step = 1; step <= num_steps; ++step) {
                int send_to = ((party - 1) + step + num_steps * i) % num_party;
                send(send_to, a, s[i], i);
            }
        }));
    }
    for (int i = 0; i < threads / 2; ++i) {
        res.push_back(pool->enqueue([this, i, num_steps, b, w, s] {
            for (int step = 1; step <= num_steps; ++step) {
                int receive_from = ((party - 1) + num_party - step - num_steps * i) % num_party;
                recv(receive_from, b, w[i], s[i + threads / 2], i);
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

    free(s);
    free(w);
    b_set = false;
}
