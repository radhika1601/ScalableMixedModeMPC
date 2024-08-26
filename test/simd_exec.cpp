#include "emp-aby/simd_interface/simd_exec.h"
int party, port;

const static int threads = 1;

double test_and(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    bool *A = new bool[length], *B = new bool[length], *C = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    prg.random_bool(B, length);
    auto start = clock_start();
    simd_circ->and_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        io->send_bool(A, length);
        io->send_bool(B, length);
        io->send_bool(C, length);
    }
    else if (party == BOB) {
        bool *A0 = new bool[length], *B0 = new bool[length], *C0 = new bool[length];
        io->recv_bool(A0, length);
        io->recv_bool(B0, length);
        io->recv_bool(C0, length);
        for (int i = 0; i < length; ++i) {
            if ((C0[i] ^ C[i]) != ((A0[i] ^ A[i]) & (B0[i] ^ B[i]))) {
                error(" Bool AND Failed!");
            };
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    io->sync();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Bool AND passed \n";
}

double test_block_and(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    block *A = new block[length], *B = new block[length], *C = new block[length];
    PRG prg;
    prg.random_block(A, length);
    prg.random_block(B, length);

    auto start = clock_start();
    simd_circ->and_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        io->send_block(A, length);
        io->send_block(B, length);
        io->send_block(C, length);
    }
    else if (party == BOB) {
        block *A0 = new block[length], *B0 = new block[length], *C0 = new block[length];
        io->recv_block(A0, length);
        io->recv_block(B0, length);
        io->recv_block(C0, length);
        xorBlocks_arr(C, C0, C, length);
        xorBlocks_arr(A, A0, A, length);
        xorBlocks_arr(B, B0, B, length);
        andBlocks_arr(A, A, B, length);
        if (!cmpBlock(A, C, length)) {
            error("Block AND Failed");
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    io->sync();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Block AND passed \n";
}

double test_not(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    bool *A = new bool[length], *B = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    auto start = clock_start();
    simd_circ->not_gate(B, A, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        io->send_bool(A, length);
        io->send_bool(B, length);
    }
    else if (party == BOB) {
        bool *A0 = new bool[length], *B0 = new bool[length];
        io->recv_bool(A0, length);
        io->recv_bool(B0, length);

        xorBools_arr(A, A, A0, length);
        xorBools_arr(B, B, B0, length);

        for (int i = 0; i < length; ++i) {
            assert(A[i] == (1 ^ B[i]));
        }
        delete[] A0;
        delete[] B0;
    }
    io->sync();
    delete[] A;
    delete[] B;
    return t;
    // std::cout << "Bool NOT gate test Passed \n";
}
double test_block_not(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    block *A = new block[length], *B = new block[length];
    PRG prg;
    prg.random_block(A, length);
    auto start = clock_start();
    simd_circ->not_gate(B, A, length);
    long long t = time_from(start);
    if (party == ALICE) {
        io->send_block(A, length);
        io->send_block(B, length);
    }
    else if (party == BOB) {
        block *A0 = new block[length], *B0 = new block[length];
        io->recv_block(A0, length);
        io->recv_block(B0, length);

        xorBlocks_arr(A, A, A0, length);
        xorBlocks_arr(B, B, B0, length);
        xorBlocks_arr(A, A, all_one_block, length);

        if (!cmpBlock(A, B, length)) {
            error("Block NOT Failed");
        }
        delete[] A0;
        delete[] B0;
    }
    delete[] A;
    delete[] B;
    return t;
    // std::cout << "Block NOT Gate tests Passed \n";
}
double test_xor(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    bool *A = new bool[length], *B = new bool[length], *C = new bool[length];
    PRG prg;
    prg.random_bool(A, length);
    prg.random_bool(B, length);

    auto start = clock_start();
    simd_circ->xor_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        io->send_bool(A, length);
        io->send_bool(B, length);
        io->send_bool(C, length);
    }
    else if (party == BOB) {
        bool *A0 = new bool[length], *B0 = new bool[length], *C0 = new bool[length];
        io->recv_bool(A0, length);
        io->recv_bool(B0, length);
        io->recv_bool(C0, length);
        for (int i = 0; i < length; ++i) {
            if ((C0[i] ^ C[i]) != ((A0[i] ^ A[i]) ^ (B0[i] ^ B[i]))) {
                error("Bool XOR Failed!");
            };
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    io->sync();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Bool XOR passed \n";
}
double test_block_xor(SIMDCircExec<NetIO>* simd_circ, NetIO* io, int length = 2000) {
    block *A = new block[length], *B = new block[length], *C = new block[length];
    PRG prg;
    prg.random_block(A, length);
    prg.random_block(B, length);

    auto start = clock_start();
    simd_circ->xor_gate(C, A, B, length);
    long long t = time_from(start);
    io->flush();
    if (party == ALICE) {
        io->send_block(A, length);
        io->send_block(B, length);
        io->send_block(C, length);
    }
    else if (party == BOB) {
        block *A0 = new block[length], *B0 = new block[length], *C0 = new block[length];
        io->recv_block(A0, length);
        io->recv_block(B0, length);
        io->recv_block(C0, length);
        xorBlocks_arr(C, C0, C, length);
        xorBlocks_arr(A, A0, A, length);
        xorBlocks_arr(B, B0, B, length);
        xorBlocks_arr(A, A, B, length);
        if (!cmpBlock(A, C, length)) {
            error("Block XOR Failed");
        }
        delete[] A0;
        delete[] B0;
        delete[] C0;
    }
    io->sync();
    delete[] A;
    delete[] B;
    delete[] C;
    return t;
    // std::cout << "Block XOR passed \n";
}

// double test_mux(SIMDCircExec *simd_circ, NetIO *io, int length = 2000,
//                 int width = 120) {
//   block *A = new block[length * width], *B = new block[length * width],
//         *C = new block[length * width];
//   bool *sel = new bool[length];
//   PRG prg;
//   prg.random_block(A, length * width);
//   prg.random_block(B, length * width);
//   prg.random_bool(sel, length);

//   auto start = clock_start();
//   simd_circ->mux_gate(C, A, B, sel, width, length);
//   long long t = time_from(start);
//   io->flush();
//   if (party == ALICE) {
//     io->send_block(A, length * width);
//     io->send_block(B, length * width);
//     io->send_block(C, length * width);
//     io->send_bool(sel, length);
//   } else if (party == BOB) {
//     block *A0 = new block[length * width], *B0 = new block[length * width],
//           *C0 = new block[length * width];
//     bool *sel0 = new bool[length];

//     io->recv_block(A0, length * width);
//     io->recv_block(B0, length * width);
//     io->recv_block(C0, length * width);
//     io->recv_bool(sel0, length);

//     xorBools_arr(sel, sel, sel0, length);
//     xorBlocks_arr(A, A, A0, length * width);
//     xorBlocks_arr(B, B, B0, length * width);
//     xorBlocks_arr(C, C, C0, length * width);
//     for (int i = 0; i < length; ++i) {
//       if (sel[i]) {
//         if (!cmpBlock(B + i * width, C + i * width, width)) {
//           error("MUX Failed");
//         }
//       } else {
//         if (!cmpBlock(A + i * width, C + i * width, width)) {
//           error("MUX Failed");
//         }
//       }
//     }
//     delete[] A0;
//     delete[] B0;
//     delete[] C0;
//     delete[] sel0;
//   }
//   io->sync();
//   delete[] A;
//   delete[] B;
//   delete[] C;
//   delete[] sel;
//   return t;
// }

int main(int argc, char** argv) {
    parse_party_and_port(argv, &party, &port);
    vector<NetIO*> ios;
    for (int i = 0; i < threads; ++i)
        ios.push_back(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port));

    auto start                     = clock_start();
    SIMDCircExec<NetIO>* simd_circ = new SIMDCircExec<NetIO>(party, threads, ios.data());
    double timeused                = time_from(start);
    std::cout << party << "\tsetup\t" << timeused / 1000 << "ms" << std::endl;
    NetIO* io = ios[0];

    // std::cout << "BOOL MUX EVALUATION\t" << test_mux(simd_circ, io)/1000 <<
    // "ms" << std::endl;
    std::cout << "BOOL AND EVALUATION\t" << test_and(simd_circ, io) / 1000 << "ms" << std::endl;
    std::cout << "BOOL XOR EVALUATION\t" << test_xor(simd_circ, io) / 1000 << "ms" << std::endl;
    std::cout << "BOOL NOT EVALUATION\t" << test_not(simd_circ, io) / 1000 << "ms" << std::endl;

    std::cout << "BLOCK AND EVALUATION\t" << test_block_and(simd_circ, io) / 1000 << "ms" << std::endl;
    std::cout << "BLOCK XOR EVALUATION\t" << test_block_xor(simd_circ, io) / 1000 << "ms" << std::endl;
    std::cout << "BLOCK NOT EVALUATION\t" << test_block_not(simd_circ, io) / 1000 << "ms" << std::endl;

    delete simd_circ;
    delete io;
}