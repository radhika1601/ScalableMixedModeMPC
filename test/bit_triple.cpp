#include "emp-aby/triple-providers/bit-triple.h"
#include "emp-aby/io/multi-io.hpp"

#include <typeinfo>
int party, port;

const static int threads = 1;

template <typename IO>
double test_get_block_triple(BitTripleProvider<IO>* bt, IO* io) {
    const int length = BitTripleProvider<IO>::BUFFER_SZ / 128;
    block* a         = new block[length];
    block* b         = new block[length];
    block* c         = new block[length];
    auto start       = clock_start();
    bt->get_triple(a, b, c);
    long long t = time_from(start);
    if (party == ALICE) {
        io->send_block(a, length);
        io->send_block(b, length);
        io->send_block(c, length);
    }
    else {
        block* a0 = new block[length];
        block* b0 = new block[length];
        block* c0 = new block[length];
        io->recv_block(a0, length);
        io->recv_block(b0, length);
        io->recv_block(c0, length);

        block* lhs = new block[length];
        block* rhs = new block[length];
        xorBlocks_arr(lhs, c0, c, length);
        block* a_xor = new block[length];
        block* b_xor = new block[length];
        xorBlocks_arr(a_xor, a0, a, length);
        xorBlocks_arr(b_xor, b0, b, length);
        andBlocks_arr(rhs, a_xor, b_xor, length);
        if (!cmpBlock(lhs, rhs, length)) {
            std::cout << "Bit Triple Failed \n";
            error("Bit Triple Failed");
        }

        delete[] a0;
        delete[] b0;
        delete[] c0;
        delete[] lhs;
        delete[] rhs;
        delete[] a_xor;
        delete[] b_xor;
    }
    delete[] a;
    delete[] b;
    delete[] c;
    return t;
}

template <typename IO>
double test_get_bool_triple(BitTripleProvider<IO>* bt, IO* io) {
    const int length = BitTripleProvider<IO>::BUFFER_SZ;
    bool* a          = new bool[length];
    bool* b          = new bool[length];
    bool* c          = new bool[length];
    memset(a, 0, length);
    memset(b, 0, length);
    memset(c, 0, length);
    auto start = clock_start();
    bt->get_triple(a, b, c);
    long long t = time_from(start);
    if (party == ALICE) {
        io->send_bool(a, length);
        io->send_bool(b, length);
        io->send_bool(c, length);
    }
    else {
        bool* a0 = new bool[length];
        bool* b0 = new bool[length];
        bool* c0 = new bool[length];
        io->recv_bool(a0, length);
        io->recv_bool(b0, length);
        io->recv_bool(c0, length);

        bool* lhs = new bool[length];
        bool* rhs = new bool[length];
        for (int i = 0; i < length; ++i) {
            lhs[i] = c0[i] ^ c[i];
            rhs[i] = (a0[i] ^ a[i]) & (b0[i] ^ b[i]);
            if (lhs[i] != rhs[i]) {
                std::cout << c0[i] << c[i] << std::endl;
                std::cout << a0[i] << a[i] << std::endl;
                std::cout << b0[i] << b[i] << std::endl;
                std::cout << "Bool Triple Failed \n";
            }
        }

        delete[] a0;
        delete[] b0;
        delete[] c0;
        delete[] lhs;
        delete[] rhs;
    }
    delete[] a;
    delete[] b;
    delete[] c;
    return t;
}

/*double test_get_mux_triple(BitTripleProvider *bt, int length = 10, int width = 12)
  {

  block *A = new block[length * width], *C = new block[length * width];
  bool *b = new bool[length];
  auto start = clock_start();
  bt->get_mux_triple(A, b, C, width, length);
  long long t = time_from(start);
  IO *io = new IO(party == ALICE ? nullptr : "127.0.0.1", port + length);
  if (party == ALICE)
  {
  io->send_block(A, length * width);
  io->send_bool(b, length);
  io->send_block(C, length * width);
  }
  else if (party == BOB)
  {
  block *A0 = new block[length * width], *C0 = new block[length * width];
  bool *b0 = new bool[length];

  io->recv_block(A0, length * width);
  io->recv_bool(b0, length);
  io->recv_block(C0, length * width);

  xorBlocks_arr(A, A0, A, length * width);
  xorBlocks_arr(C, C0, C, length * width);

  for (int i = 0; i < length; ++i)
  {
  if (b[i] ^ b0[i])
  {
  if (!cmpBlock(A + i * width, C + i * width, width))
  error("Mux Triple Failed");
  }
  else
  {
  for (int j = 0; j < width; ++j)
  if (!cmpBlock(C + i * width + j, &zero_block, 1))
  error("Mux Triple Failed");
  }
  }
  delete[] A0;
  delete[] b0;
  delete[] C0;
  }
// std::cout << "Mux triple passed. \t";
delete io;
delete[] A;
delete[] b;
delete[] C;
return t;
}*/

int main(int argc, char** argv) {
    parse_party_and_port(argv, &party, &port);

    std::vector<std::pair<std::string, unsigned short>> net_config;

    for (int i = 0; i < 2; ++i) {
        std::string s = "127.0.0.1";
        uint p        = (port + 4 * 2 * i);
        net_config.push_back(std::make_pair(s, p));
    }
    auto start = clock_start();
    vector<NetIO*> ios;

    for (int i = 0; i < threads; ++i)
        ios.push_back(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port));
    BitTripleProvider<NetIO>* bt = new BitTripleProvider<NetIO>(party, threads, ios.data());
    double timeused              = time_from(start);
    const double buffer_size     = BitTripleProvider<NetIO>::BUFFER_SZ;
    std::cout << party << "\tsetup\t" << timeused / 1000 << "ms" << std::endl;
    std::cout << "BLOCK TRIPLE GENERATION\t" << test_get_block_triple<NetIO>(bt, ios[0]) / buffer_size * 128
              << " ns per triple" << std::endl;
    std::cout << "BOOL TRIPLE GENERATION\t" << test_get_bool_triple<NetIO>(bt, ios[0]) / buffer_size << " ns per triple"
              << std::endl;
    //std::cout << "(10000, 100) MUX TRIPLE GENERATION\t" << test_get_mux_triple(bt, 10000, 100)/1000 << "ms" << std::endl;
    delete bt;
    for (auto io : ios)
        delete io;
}
