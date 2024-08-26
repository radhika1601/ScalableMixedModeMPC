#include "emp-aby/converter/b2aconverter.h"
#include "emp-aby/io/multi-io.hpp"

using namespace emp;

#include <typeinfo>

int party, port;

const static int threads = 4;

int num_party;

template <typename IO>
void check(MPIOChannel<IO>* io, bool* b, int64_t* a, int n, long long q, int l, string msg, bool mod = true) {
    int length = l * n;
    // io->sync();
    if (party == ALICE) {
        bool* tmp_b    = new bool[length];
        int64_t* tmp_a = new int64_t[n];

        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, tmp_b, length);
            xorBools_arr(b, b, tmp_b, length);

            io->recv_data(i, tmp_a, n * sizeof(int64_t));
            io->flush(i);
            for (int j = 0; j < n; ++j) {
                a[j] = (a[j] + tmp_a[j]) % q;
            }
        }
        int64_t* check = new int64_t[n];
        memset(check, 0, n * sizeof(int64_t));
        // std::cout << l << " " << q << " " << n << std::endl;
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < n; ++j) {
                if (b[j * l + i]) {
                    check[j] += (1L << i);
                }
                if (mod)
                    check[j] %= q;
            }
        }
        for (int i = 0; i < n; ++i)
            if (check[i] != a[i]) {
                std::cout << i << " " << check[i] << " " << a[i] << std::endl;
                std::cout << msg << " ";
                error("Test failed!");
            }

        // std::cout << msg << " Test passed" << std::endl;

        delete[] tmp_a;
        delete[] tmp_b;
        delete[] check;
    }
    else {
        io->send_bool(ALICE, b, length);
        io->send_data(ALICE, a, n * sizeof(int64_t));
        io->flush(ALICE);
    }
}

template <typename IO>
void test_b2a(MPIOChannel<IO>* io, B2AConverter<IO>* converter, HE<IO>* he, int n = 100, double comm_offset = 0) {
    PRG prg;
    int l         = ceil(log2(he->q));
    int length    = l * n;
    bool* boolean = new bool[length];
    prg.random_bool(boolean, length);
    int64_t* arithmetic = new int64_t[n];
    auto start          = clock_start();
    converter->convert(arithmetic, boolean, length, l);
    double timeused    = time_from(start);
    double online_comm = io->get_total_bytes_sent() - comm_offset;
    check<IO>(io, boolean, arithmetic, n, he->q, l, "B2A");
    std::cout << party << " Online time B2A conversion per 32-bit: " << timeused / (n * 1000) << " ms" << std::endl;
    std::cout << party << " Online comm B2A conversion per 32-bit: " << online_comm / (n) << " KB" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Format: b2aconverter PartyID port num_parties" << std::endl;
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
    std::cout << "io setup" << std::endl;
    ThreadPool* pool            = new ThreadPool(threads);
    const long long int modulus = (1L << 32) - (1L << 30) + 1;
    HE<MultiIOBase>* he         = new HE<MultiIOBase>(num_party, io, pool, party, modulus);
    he->multiplication_keygen();
    he->rotation_keygen();
    std::cout << party << " p = " << he->cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << party << " n = " << he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    std::cout << party
              << " log2 q = " << log2(he->cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;
    //   io->flush();
    int pool_size = std::max(20 / num_party - 1, 0) * num_party + num_party;
    if (num_party > 16)
        pool_size = 16;
    std::cout << "pool size " << pool_size << std::endl;

    // int pool_size = 20;
    // std::cout << "go in b2a constructor \n";
    auto start                           = clock_start();
    B2AConverter<MultiIOBase>* converter = new B2AConverter<MultiIOBase>(num_party, party, io, pool, he, pool_size);
    double timeused                      = time_from(start);
    double offline_comm                  = io->get_total_bytes_sent();
    int n = (he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2);

    std::cout << party << "\tB2A offline time\t" << timeused / ((pool_size * n / 32) * 1000) << " ms" << std::endl;
    std::cout << party << "\tB2A offline comm\t" << offline_comm / ((pool_size * n / 32)) << " KB" << std::endl;
    test_b2a<MultiIOBase>(io, converter, he, (min(20, pool_size) * n) / 32, offline_comm);
    delete he;
    delete io;
}
