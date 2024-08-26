#include "emp-aby/simd_interface/mp-simd-exec.h"
#include "fstream"
#include "emp-aby/io/multi-io.hpp"

int party, port;

const static int threads = 4;
int num_party;

template <typename IO>
double test_and(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 2000) {
    bool *A = new bool[length], *B = new bool[length], *C = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    prg.random_bool(B, length);
    auto start = clock_start();
    //    std::cout << "start multiplying" << std::endl;
    simd_circ->and_gate(C, A, B, length);
    //    for(int i = 0; i < length; ++i){
    //        std::cout  << i << " A: " << A[i] << " B: " << B[i] << " C: " << C[i] << std::endl;
    //    }
    long long t = time_from(start);
    io->flush();
    //    std::cout << "Now checking" << std::endl;
    io->sync();
    if (party == ALICE) {
        bool *A0 = new bool[length], *B0 = new bool[length], *C0 = new bool[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, A0, length);
            io->recv_bool(i, B0, length);
            io->recv_bool(i, C0, length);
            xorBools_arr(A, A, A0, length);
            xorBools_arr(B, B, B0, length);
            xorBools_arr(C, C, C0, length);
        }
        for (int i = 0; i < length; ++i) {
            if ((A[i] & B[i]) != C[i]) {
                std::cout << i << " ";
                error(" Bool AND Failed!");
            }
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    else {
        io->send_bool(ALICE, A, length);
        io->send_bool(ALICE, B, length);
        io->send_bool(ALICE, C, length);
    }
    io->flush();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Bool AND passed \n";
}

template <typename IO>
double test_block_and(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 2000) {
    block *A = new block[length], *B = new block[length], *C = new block[length];
    PRG prg;
    prg.random_block(A, length);
    prg.random_block(B, length);

    auto start = clock_start();
    simd_circ->and_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    //    std::cout << "Now checking" << std::endl;
    if (party == ALICE) {
        block *A0 = new block[length], *B0 = new block[length], *C0 = new block[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_block(i, A0, length);
            io->recv_block(i, B0, length);
            io->recv_block(i, C0, length);
            xorBlocks_arr(A, A, A0, length);
            xorBlocks_arr(B, B, B0, length);
            xorBlocks_arr(C, C, C0, length);
        }
        andBlocks_arr(A, A, B, length);

        if (!cmpBlock(A, C, length)) {
            error(" Bool AND Failed!");
        }

        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    else {
        io->send_block(ALICE, A, length);
        io->send_block(ALICE, B, length);
        io->send_block(ALICE, C, length);
    }
    io->flush();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Block AND passed \n";
}

template <typename IO>
double test_not(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 2000) {
    bool *A = new bool[length], *B = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    auto start = clock_start();
    simd_circ->not_gate(B, A, length);
    long long t = time_from(start);
    io->flush();
    io->sync();
    if (party == ALICE) {
        bool *A0 = new bool[length], *B0 = new bool[length], *C0 = new bool[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, A0, length);
            io->recv_bool(i, B0, length);
            xorBools_arr(A, A, A0, length);
            xorBools_arr(B, B, B0, length);
        }
        for (int i = 0; i < length; ++i) {
            if (A[i] != !B[i]) {
                std::cout << i << " ";
                error(" Bool NOT Failed!");
            }
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    else {
        io->send_bool(ALICE, A, length);
        io->send_bool(ALICE, B, length);
    }
    io->flush();
    delete[] A;
    delete[] B;
    return t;
    // std::cout << "Bool NOT gate test Passed \n";
}

template <typename IO>
double test_block_not(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 2000) {
    block *A = new block[length], *B = new block[length];
    PRG prg;
    prg.random_block(A, length);
    auto start = clock_start();
    simd_circ->not_gate(B, A, length);
    long long t = time_from(start);
    if (party == ALICE) {
        block *A0 = new block[length], *B0 = new block[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_block(i, A0, length);
            io->recv_block(i, B0, length);
            xorBlocks_arr(A, A, A0, length);
            xorBlocks_arr(B, B, B0, length);
        }
        //        andBlocks_arr(A, A, B, length);
        xorBlocks_arr(A, A, all_one_block, length);
        if (!cmpBlock(A, B, length)) {
            error(" Bool AND Failed!");
        }

        delete[] A0;
        delete[] B0;
    }
    else {
        io->send_block(ALICE, A, length);
        io->send_block(ALICE, B, length);
    }
    io->flush();
    delete[] A;
    delete[] B;
    return t;
    // std::cout << "Block NOT Gate tests Passed \n";
}

template <typename IO>
double test_xor(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 10000) {
    bool *A = new bool[length], *B = new bool[length], *C = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    prg.random_bool(B, length);

    auto start = clock_start();
    simd_circ->xor_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    io->sync();
    if (party == ALICE) {
        bool *A0 = new bool[length], *B0 = new bool[length], *C0 = new bool[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, A0, length);
            io->recv_bool(i, B0, length);
            io->recv_bool(i, C0, length);
            xorBools_arr(A, A, A0, length);
            xorBools_arr(B, B, B0, length);
            xorBools_arr(C, C, C0, length);
        }
        for (int i = 0; i < length; ++i) {
            if ((A[i] ^ B[i]) != C[i]) {
                std::cout << i << " ";
                error(" Bool AND Failed!");
            }
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    else {
        io->send_bool(ALICE, A, length);
        io->send_bool(ALICE, B, length);
        io->send_bool(ALICE, C, length);
    }
    io->flush();
    io->sync();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Bool XOR passed \n";
}

template <typename IO>
double test_block_xor(MPSIMDCircExec<IO>* simd_circ, MPIOChannel<IO>* io, int length = 2000) {
    block *A = new block[length], *B = new block[length], *C = new block[length];
    PRG prg;
    prg.random_block(A, length);
    prg.random_block(B, length);

    auto start = clock_start();
    simd_circ->xor_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        block *A0 = new block[length], *B0 = new block[length], *C0 = new block[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_block(i, A0, length);
            io->recv_block(i, B0, length);
            io->recv_block(i, C0, length);
            xorBlocks_arr(A, A, A0, length);
            xorBlocks_arr(B, B, B0, length);
            xorBlocks_arr(C, C, C0, length);
        }
        xorBlocks_arr(A, A, B, length);

        if (!cmpBlock(A, C, length)) {
            error(" Bool AND Failed!");
        }

        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    else {
        io->send_block(ALICE, A, length);
        io->send_block(ALICE, B, length);
        io->send_block(ALICE, C, length);
    }
    io->flush();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Block XOR passed \n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Format: test_mp_bit_triple PartyID Port num_parties" << std::endl;
        exit(0);
    }
    parse_party_and_port(argv, &party, &port);
    num_party = atoi(argv[3]);

    std::vector<std::pair<std::string, unsigned short>> net_config;

    if (argc == 5) {
        const char* file = argv[4];
        FILE* f          = fopen(file, "r");
        for (int i = 0; i < num_party; ++i) {
            char* c = (char*)malloc(15 * sizeof(char));
            uint p;
            fscanf(f, "%s %d\n", c, &p);
            std::string s(c);
            net_config.push_back(std::make_pair(s, p));
            fflush(f);
        }
        fclose(f);
    }
    else {
        for (int i = 0; i < num_party; ++i) {
            std::string s = "127.0.0.1";
            uint p        = (port + 4 * num_party * i);
            net_config.push_back(std::make_pair(s, p));
        }
    }

    MultiIO* io = new MultiIO(party, num_party, net_config);
    std::cout << party << " connected \n";

    ThreadPool pool(threads);
    io->setup_ot_ios();
    io->flush();
    std::cout << "io setup" << std::endl;
    long long int buffer_length            =((ferret_b13.n - ferret_b13.k - ferret_b13.t * ferret_b13.log_bin_sz_pre - 128) / 128) * 128;
    auto start                             = clock_start();
    MPSIMDCircExec<MultiIOBase>* simd_circ = new MPSIMDCircExec<MultiIOBase>(num_party, party, &pool, io);
    double timeused                        = time_from(start);
    std::cout << party << "\tsetup\t" << timeused / (2000 * buffer_length) << " ms" << std::endl;
    //    std::cout << party << " BOOL AND EVALUATION\t" << test_and(simd_circ, io, buffer_length) / 1000 << "ms"
    //              << std::endl;
    //    std::cout << party << " BOOL XOR EVALUATION\t" << test_xor(simd_circ, io) / 1000 << "ms"
    //              << std::endl;
    //    std::cout << party << " BOOL NOT EVALUATION\t" << test_not(simd_circ, io) / 1000 << "ms"
    //              << std::endl;

    std::cout << party << " BLOCK AND EVALUATION\t"
              << test_block_and(simd_circ, io, buffer_length / 128) / (buffer_length) << " us" << std::endl;
    //    std::cout << party << " BLOCK XOR EVALUATION\t" << test_block_xor(simd_circ, io) / 1000
    //              << "ms" << std::endl;
    //    std::cout << party << " BLOCK NOT EVALUATION\t" << test_block_not(simd_circ, io) / 1000
    //              << "ms" << std::endl;
    // delete simd_circ;
    delete io;
}
