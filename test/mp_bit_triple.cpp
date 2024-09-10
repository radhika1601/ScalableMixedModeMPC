#include "emp-aby/triple-providers/mp-bit-triple.h"
#include "emp-aby/io/multi-io.hpp"
using namespace emp;

#include <typeinfo>

int party, port;

const static int threads = 4;

int num_party;

template <typename IO>
double test_bool_triple(MPBitTripleProvider<IO>* mp_bit_triple_provider, MPIOChannel<IO>* io, int length = 10) {
    length  = mp_bit_triple_provider->BUFFER_SZ;
    bool *a = new bool[length], *b = new bool[length], *c = new bool[length];
    auto start = clock_start();

    mp_bit_triple_provider->get_triple(a, b, c);
    double t = time_from(start);
    io->sync();

    if (party == ALICE) {
        bool *a0 = new bool[length], *b0 = new bool[length], *c0 = new bool[length];

        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, a0, length);
            io->recv_bool(i, b0, length);
            io->recv_bool(i, c0, length);
            xorBools_arr(a, a, a0, length);
            xorBools_arr(b, b, b0, length);
            xorBools_arr(c, c, c0, length);
        }

        andBools_arr(c0, a, b, length);

        for (int i = 0; i < length; ++i) {
            if (c[i] != c0[i]) {
                std::cout << i << " " << a[i] << " " << b[i] << " " << c[i] << " " << c0[i] << std::endl;
                error("Bool Triple failed!");
            }
        }

        //        std::cout << "Bool Triple Passed" << std::endl;
    }
    else {
        io->send_bool(ALICE, a, length);
        io->send_bool(ALICE, b, length);
        io->send_bool(ALICE, c, length);
    }
    io->flush();
    return t;
}

template <typename IO>
double test_block_triple(MPBitTripleProvider<IO>* mp_bit_triple_provider, MPIOChannel<IO>* io) {
    int length = mp_bit_triple_provider->BUFFER_SZ / 128;
    block *a = new block[length], *b = new block[length], *c = new block[length];
    auto start = clock_start();
    mp_bit_triple_provider->get_triple(a, b, c);
    double t = time_from(start);
    io->sync();

    if (party == ALICE) {
        block *a0 = new block[length], *b0 = new block[length], *c0 = new block[length];

        for (int i = 2; i <= num_party; ++i) {
            io->recv_block(i, a0, length);
            io->recv_block(i, b0, length);
            io->recv_block(i, c0, length);
            xorBlocks_arr(a, a, a0, length);
            xorBlocks_arr(b, b, b0, length);
            xorBlocks_arr(c, c, c0, length);
        }

        andBlocks_arr(c0, a, b, length);

        if (!cmpBlock(c, c0, length)) {
            error("block Triple failed!");
        }

        //        std::cout << "block Triple Passed" << std::endl;
    }
    else {
        io->send_block(ALICE, a, length);
        io->send_block(ALICE, b, length);
        io->send_block(ALICE, c, length);
    }
    io->flush();
    return t;
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
    auto start = clock_start();

    MPBitTripleProvider<MultiIOBase>* mp_bit_triple_provider =
        new MPBitTripleProvider<MultiIOBase>(num_party, party, &pool, io);
    double timeused = time_from(start);
    std::cout << party << "\tsetup\t" << timeused / 1000 << "ms" << std::endl;
    //    const double buffer_size = pow(128, 3);

    std::cout << party << " BOOL TRIPLE GENERATION\t" << test_bool_triple(mp_bit_triple_provider, io) / (1000) << " ms"
              << std::endl;
    std::cout << party << " BLOCK TRIPLE GENERATION\t" << test_block_triple(mp_bit_triple_provider, io) / (1000)
              << " ms" << std::endl;

    delete io;
}
