#include "emp-aby/mp-circuit.hpp"
#include <typeinfo>

int party, port;

const static int threads = 1;

long long mpc_main(long long INPUT_A, long long INPUT_B, long long INPUT_B2) {
    long long value = (INPUT_A + INPUT_B) & ((1LL << 61) - 1);
    long long relu  = value < ((1LL << 60) - 1) ? value : 0;
    return relu + ((1LL << 61) - 1) - INPUT_B2;
}

void test_circuit(const char* file, int party, NetIO** ios, const int num, SIMDCircExec<NetIO>* simd_circ) {
    PRG prg;
    Circuit<SIMDCircExec<NetIO>>* circuit = new Circuit<SIMDCircExec<NetIO>>(file, party, simd_circ);

    bool* in = new bool[(circuit->n1 + circuit->n2) * num];

    bool* out = new bool[(circuit->n3) * num];

    for (int i = 0; i < 5; ++i) {
        prg.random_bool(in, (circuit->n1 + circuit->n2) * num);

        auto start = clock_start();
        circuit->compute<NetIO>(out, in, num);
        long long t = time_from(start);

        NetIO* io = ios[0];
        if (party == ALICE) {
            io->send_bool(in, (circuit->n1) * num);
            io->send_bool(out, (circuit->n3) * num);
        }
        else {
            bool* in1  = new bool[(circuit->n1) * num];
            bool* out1 = new bool[(circuit->n3) * num];
            io->recv_bool(in1, (circuit->n1) * num);
            io->recv_bool(out1, (circuit->n3) * num);
            xorBools_arr(out, out1, out, (circuit->n3) * num);
            for (int i = 0; i < num; ++i) {
                long long A          = bool_to_int<long long>(in1 + i * circuit->n1);
                long long B          = bool_to_int<long long>(in + i * circuit->n2);
                long long B2         = bool_to_int<long long>(in + i * circuit->n2 + 64);
                long long output     = bool_to_int<long long>(out + i * circuit->n3);
                long long sim_output = mpc_main(A, B, B2);
                if (output != sim_output) {
                    error("Test Failed!");
                }
            }
        }
        std::cout << "Compute Done: " << t << " us" << std::endl;
    }
    delete[] in;
    delete[] out;
}

int main(int argc, char** argv) {
    parse_party_and_port(argv, &party, &port);
    const char* file = argv[3];
    vector<NetIO*> ios;
    for (int i = 0; i < threads; ++i)
        ios.push_back(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port));

    SIMDCircExec<NetIO>* simd_circ = new SIMDCircExec<NetIO>(party, threads, ios.data());
    test_circuit(file, party, ios.data(), 10000, simd_circ);
    std::cout << "num and gates " << simd_circ->num_and_gates << " " << simd_circ->depth << "\n";
    delete simd_circ;
    for (auto io : ios)
        delete io;
}
