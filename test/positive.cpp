#include "emp-aby/mp-circuit.hpp"
#include <typeinfo>

int party, port;

const static int threads = 1;

unsigned char mpc_main(long long INPUT_A, long long INPUT_B) {
    long long value = (INPUT_A + INPUT_B) & ((1ULL << 61) - 1);
    return value > ((1ULL << 60) - 1) ? 1 : 0;
}

void test_circuit(const char* file, int party, NetIO** ios, const int num, SIMDCircExec<NetIO>* simd_circ) {
    PRG prg;
    Circuit<SIMDCircExec<NetIO>>* circuit = new Circuit<SIMDCircExec<NetIO>>(file, party, simd_circ);
    bool* in;
    if (party == ALICE) {
        in = new bool[(circuit->n1) * num];
        prg.random_bool(in, (circuit->n1) * num);
    }
    else {
        in = new bool[(circuit->n2) * num];
        prg.random_bool(in, (circuit->n2) * num);
    }

    bool* out  = new bool[(circuit->n3) * num];
    auto start = clock_start();
    for (int j = 0; j < 50; ++j) {
        circuit->compute<NetIO>(out, in, num);
    }
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
            long long A     = bool_to_int<long long>(in1 + i * circuit->n1);
            long long B     = bool_to_int<long long>(in + i * circuit->n2);
            bool output     = out[i];
            bool sim_output = mpc_main(A, B);
            if (output != sim_output) {
                error("Test Failed!");
            }
        }
    }

    std::cout << "Compute Done: " << t << " us" << std::endl;
}
int main(int argc, char** argv) {
    parse_party_and_port(argv, &party, &port);
    const char* file = argv[3];
    vector<NetIO*> ios;
    for (int i = 0; i < threads; ++i)
        ios.push_back(new NetIO(party == ALICE ? nullptr : "127.0.0.1", port));

    SIMDCircExec<NetIO>* simd_circ = new SIMDCircExec<NetIO>(party, threads, ios.data());
    test_circuit(file, party, ios.data(), 10000, simd_circ);
    delete simd_circ;
}
