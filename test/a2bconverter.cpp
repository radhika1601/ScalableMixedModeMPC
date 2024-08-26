#include "emp-aby/converter/a2bconverter.h"
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
void test_a2b(MPIOChannel<IO>* io, A2BConverter<IO>* converter, HE<IO>* he, size_t num = 1, double comm_offset = 0) {
    PRG prg;
    const long long int q = he->q;
    const int l           = ceil(log2(q));
    int64_t* a            = new int64_t[num];
    prg.random_data(a, num * sizeof(int64_t));
    for (int i = 0; i < num; ++i) {
        a[i] = ((a[i] % he->q) + he->q) % he->q;
    }

    bool* b    = new bool[num * l];
    auto start = clock_start();
    converter->convert(b, a, num);
    double timeused    = time_from(start);
    double online_comm = io->get_total_bytes_sent() - comm_offset;
    check<IO>(io, b, a, num, q, l, "A2B");
    std::cout << party << " A2B online time conversion per 32-bit: " << timeused / (num * 1000) << " ms\t" << std::endl;
    std::cout << party << " A2B online comm conversion per 32-bit: " << online_comm / (num) << " KB\t" << std::endl;
}

template <typename IO>
void test_rand_ab_shares(MPIOChannel<IO>* io, A2BConverter<IO>* converter, HE<IO>* he, size_t length = 1) {
    const int l = ceil(log2(he->q));
    bool* b     = new bool[length * l];
    int64_t* a  = new int64_t[length];

    int waste = converter->rand_ab_shares(a, b, length);
    check(io, b, a, length - waste, he->q, l, "Random shares gen", false);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Format: a2bconverter PartyID port num_parties" << std::endl;
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
    io->setup_ot_ios();
    std::cout << "io setup" << std::endl;

    size_t BUFFER_SZ                       = ((ferret_b13.n - ferret_b13.k - ferret_b13.t * ferret_b13.log_bin_sz_pre - 128) / 128) * 128;
    auto start                             = clock_start();
    MPSIMDCircExec<MultiIOBase>* simd_circ = new MPSIMDCircExec<MultiIOBase>(num_party, party, &pool, io);
    double timeused                        = time_from(start);
    double triple_comm                     = io->get_total_bytes_sent();
    // std::cout << party << "\tTriple Generation time\t" << timeused / (2 * BUFFER_SZ * 1000) << " ms\t" << std::endl;
    // std::cout << party << "\tTriple Generation comm\t" << triple_comm / (2 * BUFFER_SZ * 1000) << " KB\t" << std::endl;
    double per_triple_time      = timeused / (2 * BUFFER_SZ * 1000);
    double per_triple_comm      = triple_comm / (2 * BUFFER_SZ);
    const long long int modulus = (1L << 32) - (1L << 30) + 1;
    HE<MultiIOBase>* he         = new HE<MultiIOBase>(num_party, io, &pool, party, modulus);
    he->multiplication_keygen();
    he->rotation_keygen();
    std::cout << "p = " << he->cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2 << std::endl;
    std::cout << "log2 q = " << log2(he->cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;
    int pool_size = std::max(20 / num_party - 1, 0) * num_party + num_party;
    if (num_party > 16)
        pool_size = 16;
    std::cout << "pool size " << pool_size << std::endl;
    start = clock_start();
    A2BConverter<MultiIOBase>* converter =
        new A2BConverter<MultiIOBase>(num_party, party, io, &pool, he, simd_circ, pool_size);
    timeused            = time_from(start);
    double offline_comm = io->get_total_bytes_sent();
    int offline_pool    = (converter->ab_share_pool - converter->num_rejected);
    std::cout << party << "\tA2B offline\t" << timeused / (offline_pool * 1000) + (per_triple_time * 325) << " ms\t"
              << std::endl;
    std::cout << party << "\tA2B offline comm\t"
              << ((offline_comm - triple_comm) / (offline_pool)) + (per_triple_comm * 325) << " KB\t" << std::endl;

    int online_pool = min(converter->ab_share_pool - converter->num_rejected, BUFFER_SZ / 325);
    test_a2b<MultiIOBase>(io, converter, he, online_pool, offline_comm);
    delete he;
    delete io;
}
