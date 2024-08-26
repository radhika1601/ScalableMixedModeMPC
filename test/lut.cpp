#include "emp-aby/lut.h"
#include "emp-aby/io/multi-io.hpp"

using namespace emp;

#include <typeinfo>

int party, port;

const static int threads = 4;

int num_party;

template <typename IO>
void test_generate_shares(HE<IO>* he, LUT<IO>* lut, MPIOChannel<IO>* io, int n = 100) {
    n                  = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    int64_t* lut_share = new int64_t[2 * n];
    bool* rotation     = new bool[n];
    int64_t table[2] = {0, 1};
    lut->generate_shares(lut_share, rotation, n, table);
    if (party == ALICE) {
        bool* tmp_rot = new bool[n];
        int64_t* tmp  = new int64_t[2 * n];

        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, tmp_rot, n);
            io->recv_data(i, (int64_t*)tmp, 2 * n * sizeof(int64_t));

            for (int j = 0; j < 2 * n; ++j) {
                lut_share[j] = (he->q + lut_share[j] + tmp[j]) % he->q;
            }
            xorBools_arr(rotation, tmp_rot, rotation, n);
        }
        for (int i = 0; i < n; ++i) {
            // std::cout << rotation[i] << " " << lut_share[2 * i] << " " << lut_share[2 * i + 1] << "\n";
            if (table[0] != lut_share[2 * i + rotation[i] * 1])
                error("Failed lut share generation!");
            if (table[1] != lut_share[2 * i + (!rotation[i]) * 1])
                error("Failed lut share generation!");
        }
        std::cout << "LUT Shares generated!" << std::endl;
    }
    else {
        io->send_bool(ALICE, rotation, n);
        io->send_data(ALICE, lut_share, 2 * n * sizeof(int64_t));
    }
}

template <typename IO>
void test_lookup(MPIOChannel<IO>* io, LUT<IO>* lut, HE<IO>* he, int64_t t[2], int n = 100) {
    bool* in     = new bool[n];
    int64_t* out = new int64_t[n];
    PRG prg;
    prg.random_bool(in, n);
    auto start = clock_start();
    lut->lookup(out, in, n);
    long long timeused = time_from(start);

    if (party == ALICE) {
        bool* tmp        = new bool[n];
        int64_t* tmp_out = new int64_t[n];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, tmp, n);
            io->recv_data(i, tmp_out, n * sizeof(int64_t));
            xorBools_arr(in, in, tmp, n);

            for (size_t j = 0; j < n; ++j) {
                out[j] = (he->q + out[j] + tmp_out[j]) % he->q;
            }
        }

        for (int i = 0; i < n; ++i) {
            if (t[(int)in[i]] != out[i]) {
                std::cout << t[(int)in[i]] << " " << out[i] << "\n";
                error("Lookup failed!");
            }
        }
        // std::cout << "Lookup test passed " << n << std::endl;
    }
    else {
        io->send_bool(ALICE, in, n);
        io->send_data(ALICE, out, n * sizeof(int64_t));
        io->flush(ALICE);
    }
    std::cout << party << " Compute Done " << n << ": " << timeused << " micro sec" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Format: lut PartyID port num_parties" << std::endl;
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

    ThreadPool pool(threads);

    MultiIO* io = new MultiIO(party, num_party, net_config);
    // io->extra_io();
    std::cout << "io setup" << std::endl;

    auto start                  = clock_start();
    const long long int modulus = (1L << 32) - (1 << 30) + 1;

    HE<MultiIOBase>* he = new HE<MultiIOBase>(num_party, io, &pool, party, modulus);
    he->rotation_keygen();
    int n = he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    // std::cout << "log2 q = " << log2(he->cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
    //           << std::endl;
    int64_t table[2] = {0, 1};

    LUT<MultiIOBase>* lut = new LUT<MultiIOBase>(num_party, party, io, &pool, he, table);
    long long timeused    = time_from(start);
    std::cout << party << "\tsetup\t" << timeused / 1000 << "ms" << std::endl;

    test_lookup(io, lut, he, table, 1 * n);
    test_generate_shares(he, lut, io, n);
    delete io;
}
