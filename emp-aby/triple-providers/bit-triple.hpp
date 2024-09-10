template <typename IO>
long long int BitTripleProvider<IO>::BUFFER_SZ =
    ((ferret_b13.n - ferret_b13.k - ferret_b13.t * ferret_b13.log_bin_sz - 128) / 128) * 128;

template <typename IO>
void BitTripleProvider<IO>::compute_rcots(bool* b) {
    auto t = clock_start();
    ot0->rcot_inplace(r0.data(), ferret_b13.n);
    ot1->rcot_inplace(r1.data(), ferret_b13.n);
    //	duplex_rcot_inplace(ot0, ot1, r0.data(),r1.data(), N_REG, N_REG);
    // std::cout << "1: " << time_from(t) << "\n";
    t = clock_start();
    if (party == ALICE) {
        for (int i = 0; i < BUFFER_SZ; i += 128) {
            this->crh.Hn(A_hat.data() + i, r0.data() + i, 128, scratch.data());
            xorBlocks_arr(A_star.data() + i, r0.data() + i, this->delta, 128);
            this->crh.Hn(A_star.data() + i, A_star.data() + i, 128, scratch.data());
            this->crh.Hn(B_hat.data() + i, r1.data() + i, 128, scratch.data());
            for (int j = 0; j < 128; ++j)
                b[j + i] = getLSB(r1[j + i]);
        }
    }
    else {  // BOB
        for (int i = 0; i < BUFFER_SZ; i += 128) {
            this->crh.Hn(B_hat.data() + i, r0.data() + i, 128, scratch.data());
            for (int j = 0; j < 128; ++j)
                b[j + i] = getLSB(r0[j + i]);
            this->crh.Hn(A_hat.data() + i, r1.data() + i, 128, scratch.data());
            xorBlocks_arr(A_star.data() + i, r1.data() + i, this->delta, 128);
            this->crh.Hn(A_star.data() + i, A_star.data() + i, 128, scratch.data());
        }
    }
    // std::cout << "2: " << time_from(t) << "\n";
}

template <typename IO>
BitTripleProvider<IO>::BitTripleProvider(int party, int threads, IO** ios) {
    this->party                            = party;
    this->threads                          = threads;
    this->ios                              = ios;
    static std::string ot0_pre_ot_filename = (party == ALICE ? "./data/ALICE_sender" : "./data/BOB_receiver");
    this->ot0     = new FerretCOT<IO>(party, threads, this->ios, false, true, ferret_b13, ot0_pre_ot_filename);
    int rev_party = party == ALICE ? BOB : ALICE;
    static std::string ot1_pre_ot_filename = (party == ALICE ? "./data/ALICE_receiver" : "./data/BOB_sender");
    this->ot1   = new FerretCOT<IO>(rev_party, threads, this->ios, false, true, ferret_b13, ot1_pre_ot_filename);
    BUFFER_SZ   = ((ferret_b13.n - ferret_b13.k - ferret_b13.t * ferret_b13.log_bin_sz - 128) / 128) * 128;
    this->delta = (party == ALICE) ? this->ot0->Delta : this->ot1->Delta;
    r0.resize(emp::ferret_b13.n);
    r1.resize(emp::ferret_b13.n);
    scratch.resize(BUFFER_SZ);
    A_hat.resize(BUFFER_SZ);
    A_star.resize(BUFFER_SZ);
    B_hat.resize(BUFFER_SZ);
    a_bool.resize(BUFFER_SZ);
    b_bool.resize(BUFFER_SZ);
    c_bool.resize(BUFFER_SZ);
}

template <typename IO>
BitTripleProvider<IO>::~BitTripleProvider() {
    delete ot0;
    delete ot1;
}

template <typename IO>
void BitTripleProvider<IO>::get_triple(block* a, block* b, block* c) {
    this->get_triple((bool*)a_bool.data(), (bool*)b_bool.data(), (bool*)c_bool.data());
    bool_to_block_arr(a, (bool*)a_bool.data(), BUFFER_SZ);
    bool_to_block_arr(b, (bool*)b_bool.data(), BUFFER_SZ);
    bool_to_block_arr(c, (bool*)c_bool.data(), BUFFER_SZ);
}

template <typename IO>
void BitTripleProvider<IO>::get_triple(bool* a, bool* b, bool* c) {
    this->compute_rcots(b);
    //	auto t = clock_start();
    for (int i = 0; i < BUFFER_SZ; ++i) {
        a[i] = getLSB(A_hat[i]) ^ getLSB(A_star[i]);
        c[i] = (a[i] & b[i]) ^ getLSB(B_hat[i]) ^ getLSB(A_hat[i]);
    }
    //	std::cout<<"3: "<< time_from(t)<<"\n";
}

// bit-vector multiplication. AKA multiplexier
// \Xor C[i,...,i+wdith-1] = (\xor b[i])(\Xor A [i,...,i+width-1]) for i in[0, length)
/* template <typename IO>
void BitTripleProvider<IO>::get_mux_triple(block* A, bool* b, block* C, int width, int length) {
    int num      = length * width;
    block *A_hat = new block[num], *A_star = new block[num];
    block* B_hat = new block[num];
    bool* b_all  = new bool[num];
    this->compute_rcots(b_all);
    block* r = new block[num];
    for (int i = 0; i < num; ++i)
        if (b_all[i])
            r[i] = A_star[i];
        else
            r[i] = A_hat[i];

    xorBlocks_arr(A, A_hat, A_star, num);
    xorBlocks_arr(C, B_hat, r, num);

    bool* xors = new bool[num];

    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < width; ++j) {
            xors[i * width + j] = b_all[i * width] ^ b_all[i * width + j];
        }
    }
    if (party == ALICE) {
        this->ios[0]->send_bool(xors, num);
        this->ios[0]->recv_bool(xors, num);
    }
    else if (party == BOB) {
        bool* xors_a = new bool[num];
        this->ios[0]->recv_bool(xors_a, num);
        xorBools_arr(xors, xors, xors_a, num);
        delete[] xors_a;
        this->ios[0]->send_bool(xors, num);
    }
    this->ios[0]->sync();
    for (int i = 0; i < length; ++i) {
        b[i] = b_all[i * width];
        for (int j = 0; j < width; ++j) {
            if (xors[i * width + j]) {
                C[i * width + j] = C[i * width + j] ^ A[i * width + j];
            }
        }
    }

    delete[] A_hat;
    delete[] B_hat;
    delete[] b_all;
    delete[] A_star;
    delete[] xors;
    delete[] r;
}
*/
