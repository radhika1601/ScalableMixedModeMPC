#include "emp-aby/simd_interface/arithmetic-circ.h"
#include "emp-aby/io/multi-io.hpp"

using namespace emp;

#include <typeinfo>

int party, port;

const static int threads = 4;

int num_party;

template <typename IO>
void test_triple(MPIOChannel<IO>* io, HE<IO>* he, ArithmeticCirc<IO>* circ) {
    int n = circ->num_triples_pool;
    int64_t *triple_a, *triple_b, *triple_c;

    triple_a = new int64_t[n];
    triple_b = new int64_t[n];
    triple_c = new int64_t[n];

    circ->get_triples(triple_a, triple_b, triple_c);

    if (party == ALICE) {
        int64_t *a, *b, *c;

        a = new int64_t[n];
        b = new int64_t[n];
        c = new int64_t[n];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_data(i, a, n * sizeof(int64_t));
            io->recv_data(i, b, n * sizeof(int64_t));
            io->recv_data(i, c, n * sizeof(int64_t));

            for (int j = 0; j < n; ++j) {
                triple_a[j] = (he->q + triple_a[j] + a[j]) % he->q;
                triple_b[j] = (he->q + triple_b[j] + b[j]) % he->q;
                triple_c[j] = (he->q + triple_c[j] + c[j]) % he->q;
            }
        }
        for (int i = 0; i < n; ++i) {
            uint64_t x = ((uint64_t)triple_a[i] * (uint64_t)triple_b[i]);
            if (triple_c[i] != x % he->q) {
                std::cout << i << " " << triple_a[i] << " " << triple_b[i] << " " << triple_c[i] << " " << x % he->q
                          << std::endl;
                error("Arithmetic triple failed!");
            }
        }
        std::cout << "Arithmetic triple test passed!" << std::endl;
    }
    else {
        io->send_data(ALICE, triple_a, n * sizeof(int64_t));
        io->send_data(ALICE, triple_b, n * sizeof(int64_t));
        io->send_data(ALICE, triple_c, n * sizeof(int64_t));
    }
    io->flush();
}

template <typename IO>
void test_mul(MPIOChannel<IO>* io, HE<IO>* he, ArithmeticCirc<IO>* circ, int length = 1000) {
    int64_t *in1 = new int64_t[length], *in2 = new int64_t[length], *out = new int64_t[length];
    PRG prg;
    prg.random_data(in1, length * sizeof(int64_t));
    prg.random_data(in2, length * sizeof(int64_t));
    for (size_t i = 0; i < length; ++i) {
        in1[i] %= he->q;
        in2[i] %= he->q;
        in1[i] = (he->q + in1[i]) % he->q;
        in2[i] = (he->q + in2[i]) % he->q;
    }
    auto start = clock_start();
    circ->mult(out, in1, in2, length);
    double timeused = time_from(start);
    if (party == ALICE) {
        int64_t *a, *b, *c;

        a = new int64_t[length];
        b = new int64_t[length];
        c = new int64_t[length];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_data(i, a, length * sizeof(int64_t));
            io->recv_data(i, b, length * sizeof(int64_t));
            io->recv_data(i, c, length * sizeof(int64_t));

            for (int j = 0; j < length; ++j) {
                in1[j] = (he->q + in1[j] + a[j]) % he->q;
                in2[j] = (he->q + in2[j] + b[j]) % he->q;
                out[j] = (he->q + out[j] + c[j]) % he->q;
            }
        }
        for (int i = 0; i < length; ++i) {
            uint64_t x = ((uint64_t)in1[i] * (uint64_t)in2[i]);
            if (out[i] != x % he->q) {
                std::cout << i << " " << out[i] << " " << in1[i] << " " << in2[i] << std::endl;
                error("Arithmetic multiplication failed!");
            }
        }
    }
    else {
        io->send_data(ALICE, in1, length * sizeof(int64_t));
        io->send_data(ALICE, in2, length * sizeof(int64_t));
        io->send_data(ALICE, out, length * sizeof(int64_t));
    }
    std::cout << "Arith. Mult Online:\t" << timeused / (1000 * length) << " ms" << std::endl;
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
    io->flush();
    std::cout << "io setup" << std::endl;
    const long long int modulus = (1L << 32) - (1 << 30) + 1;

    HE<MultiIOBase>* he = new HE<MultiIOBase>(num_party, io, &pool, party, modulus, 1, true, false, true);
    he->multiplication_keygen();

    std::cout << "p = " << he->cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2 << std::endl;
    std::cout << "log2 q = " << log2(he->cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    auto start                        = clock_start();
    ArithmeticCirc<MultiIOBase>* circ = new ArithmeticCirc<MultiIOBase>(num_party, party, io, he);
    double timeused                   = time_from(start);
    std::cout << "Arith triple gen\t" << timeused / (1000 * circ->num_triples_pool) << " ms" << std::endl;
    test_mul<MultiIOBase>(io, he, circ, circ->num_triples_pool);
    // delete he;
    delete io;
}
