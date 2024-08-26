#include "emp-aby/mp-circuit.hpp"
#include "emp-aby/io/multi-io.hpp"

#include <typeinfo>

int party, port, num_party;

const static int threads = 4;
int64_t mod              = ((1UL << 32) - (1UL << 30) + 1);

uint32_t mpc_main(uint32_t INPUT_A, uint32_t INPUT_B) {
    uint64_t x = (uint64_t)INPUT_A + (uint64_t)INPUT_B;

    if (x > mod)
        x = x - mod;
    return (uint32_t)x;
}

template <typename IO>
void test_circuit(const char* file, int party, MPIOChannel<IO>* io, const int num, MPSIMDCircExec<IO>* simd_circ,
                  Circuit<MPSIMDCircExec<IO>>* circuit) {
    PRG prg;

    bool* in = new bool[(circuit->n1 + circuit->n2) * num];
    prg.random_bool(in, (circuit->n1 + circuit->n2) * num);

    bool* out = new bool[(circuit->n3) * num];
    circuit->template compute<IO>(out, in, num);
    io->sync();
    if (party == ALICE) {
        bool* in_tmp  = new bool[(circuit->n1 + circuit->n2) * num];
        bool* out_tmp = new bool[(circuit->n3) * num];
        for (int i = 2; i <= num_party; ++i) {
            io->recv_bool(i, in_tmp, (circuit->n1 + circuit->n2) * num);
            io->recv_bool(i, out_tmp, (circuit->n3) * num);
            xorBools_arr(in, in, in_tmp, (circuit->n1 + circuit->n2) * num);
            xorBools_arr(out, out, out_tmp, circuit->n3 * num);
            io->flush(i);
        }
        for (int i = 0; i < num; ++i) {
            uint32_t A          = bool_to_int<uint32_t>(in + i * circuit->n1);
            uint32_t B          = bool_to_int<uint32_t>(in + num * circuit->n1 + i * circuit->n2);
            uint32_t output     = bool_to_int<uint32_t>(out + i * circuit->n3);
            uint32_t sim_output = mpc_main(A, B);
            if ((output + mod) % mod != (sim_output + mod) % mod) {
                std::cout << A << " " << B << " " << output << " " << sim_output << "\n";
                error("Test Failed!");
            }
        }
        std::cout << "Test passed" << std::endl;
    }
    else {
        io->send_bool(ALICE, in, (circuit->n1 + circuit->n2) * num);
        io->send_bool(ALICE, out, (circuit->n3) * num);
        io->flush(ALICE);
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cout << "Format: test_mp_bit_triple PartyID Port Circuit_file num_parties" << std::endl;
        exit(0);
    }
    parse_party_and_port(argv, &party, &port);
    const char* file = argv[3];
    num_party        = atoi(argv[4]);

    std::vector<std::pair<std::string, unsigned short>> net_config;

    if (argc == 6) {
        const char* file = argv[5];
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
    MPSIMDCircExec<MultiIOBase>* simd_circ = new MPSIMDCircExec<MultiIOBase>(num_party, party, &pool, io);
    std::cout << "SIMD CIRC setup \n";
    Circuit<MPSIMDCircExec<MultiIOBase>>* circuit = new Circuit<MPSIMDCircExec<MultiIOBase>>(file, party, simd_circ);
    std::cout << "Circuit setup done \n";

    test_circuit(file, party, io, 5000, simd_circ, circuit);
    test_circuit(file, party, io, 10000, simd_circ, circuit);
    delete io;
}
