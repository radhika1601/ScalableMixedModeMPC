#pragma once

#include "emp-aby/wire.h"
#include "emp-aby/simd_interface/mp-simd-exec.h"

#include <fstream>
#define _debug
namespace emp {

template <typename SIMDCirc>
class Circuit {
public:
    uint num_gates, num_wires, n1, n2, n3, party;
    vector<vector<int>> level_map;
    vector<Wire*> circuit;
    SIMDCirc* simd_circ;
    ~Circuit() {}
    Circuit(FILE* f) {
        this->from_file(f);
    }

    Circuit(const char* file, int party, SIMDCirc* simd_circ) {
        this->party = party;
        this->from_file(file);
        this->simd_circ = simd_circ;
    }

    void from_file(const char* file) {
        FILE* f = fopen(file, "r");
        if (f == nullptr) {
            std::cout << file << "\n";
            error("could not open file\n");
        }
        this->from_file(f);
        fclose(f);
    }

    void add_input_wires() {
        for (uint i = 0; i < (n1 + n2); ++i) {
            Wire* w    = new Wire(INPUT);
            circuit[i] = w;
            insert_level_map(w, i);
        }
    }

    void from_file(FILE* f) {
        circuit.clear();
        level_map.clear();
        int tmp = 0, in1 = 0, in2 = 0, out = 0, type = 0;
        fscanf(f, "%d%d\n", &num_gates, &num_wires);
        for (uint i = 0; i < num_wires; ++i) {
            Wire* w = nullptr;
            circuit.push_back(w);
        }
        fflush(f);
        fscanf(f, "%d%d%d\n", &n1, &n2, &n3);
        fflush(f);
        fscanf(f, "\n");
        fflush(f);
        this->add_input_wires();
        char str[10];
        Wire* new_wire;
        for (uint i = 0; i < num_gates; ++i) {
            fscanf(f, "%d", &tmp);
            if (tmp == 2) {
                fscanf(f, "%d%d%d%d%s", &tmp, &in1, &in2, &out, str);
                fflush(f);
                if (str[0] == 'A')
                    type = AND;
                else if (str[0] == 'X')
                    type = XOR;
                Wire* in1_wire = circuit[in1];
                Wire* in2_wire = circuit[in2];
                new_wire       = new Wire(type, in1_wire, in2_wire);
                circuit[out]   = new_wire;
                this->insert_level_map(new_wire, out);
            }
            else if (tmp == 1) {
                (void)fscanf(f, "%d%d%d%s", &tmp, &in1, &out, str);
                type           = INV;
                Wire* in1_wire = circuit[in1];
                new_wire       = new Wire(type, in1_wire);
                circuit[out]   = new_wire;
                this->insert_level_map(new_wire, out);
            }
        }
    }

    void insert_level_map(Wire* new_wire, int out) {
        if (level_map.size() < new_wire->level) {
            error("Missing level!");
        }
        else if (level_map.size() == new_wire->level) {
            level_map.push_back(std::vector<int>(1, out));
        }
        else {
            level_map[new_wire->level].push_back(out);
        }
    }

    template <typename IO>
    void compute(bool* out, bool* in, uint num, bool shared = false) {
#ifdef __debug
        auto t = clock_start();
#endif
        uint block_num = num / 128, bool_num = num % 128;
        Wire* w;
        if (!shared && std::is_same<SIMDCirc, SIMDCircExec<IO>>::value && (party == BOB)) {
            for (uint i = 0; i < n1; ++i) {
                w = circuit[i];
                w->initialise_value(block_num, bool_num);
                w->set = true;
            }
            for (uint i = 0; i < n2; ++i) {
                w = circuit[i + n1];
                w->initialise_value(block_num, bool_num);
                bool dummy[128] = {false};
                uint j;
                for (j = 0; j < num; ++j) {
                    dummy[j % 128] = in[i + j * n2];
                    if (j % 128 == 127) {
                        w->value[j / 128] = bool_to_block(dummy);
                    }
                }
                if (j % 128 != 0) {
                    for (uint k = 0; k < bool_num; ++k) {
                        w->rem_value[k] = dummy[k];
                    }
                }
                w->set = true;
            }
        }
        else {
            for (uint i = 0; i < n1; ++i) {
                w = circuit[i];
                w->initialise_value(block_num, bool_num);
                bool dummy[128] = {false};
                uint j;
                for (j = 0; j < num; ++j) {
                    dummy[j % 128] = in[i + j * n1];
                    if (j % 128 == 127) {
                        w->value[j / 128] = bool_to_block(dummy);
                    }
                }
                if (j % 128 != 0) {
                    for (uint k = 0; k < bool_num; ++k) {
                        w->rem_value[k] = dummy[k];
                    }
                }
                w->set = true;
            }
            if (!shared && std::is_same<SIMDCirc, SIMDCircExec<IO>>::value && (party == ALICE)) {
                for (uint i = 0; i < n2; ++i) {
                    w = circuit[i + n1];
                    w->initialise_value(block_num, bool_num);
                    w->set = true;
                }
            }
            else {
                for (uint i = 0; i < n2; ++i) {
                    w = circuit[i + n1];
                    w->initialise_value(block_num, bool_num);
                    bool dummy[128] = {false};
                    uint j;
                    for (j = 0; j < num; ++j) {
                        dummy[j % 128] = in[i + j * n2 + n1 * num];
                        if (j % 128 == 127) {
                            w->value[j / 128] = bool_to_block(dummy);
                        }
                    }
                    if (j % 128 != 0) {
                        for (uint k = 0; k < bool_num; ++k) {
                            w->rem_value[k] = dummy[k];
                        }
                    }
                    w->set = true;
                }
            }
        }
#ifdef __debug
        std::cout << "Inputs: " << time_from(t) << " us\n";
#endif
        for (uint i = 0; i < level_map.size(); ++i) {
#ifdef __debug
            t = clock_start();
#endif
            block *in1_and = nullptr, *in2_and = nullptr, *out_and = nullptr;
            bool *bool_in1_and = nullptr, *bool_in2_and = nullptr, *bool_out_and = nullptr;
            int num_ands = 0;
            for (auto it = level_map[i].begin(); it != level_map[i].end(); ++it) {
                Wire* w = circuit[*it];
                switch (w->type) {
                    case AND:
                        if (w->in1->set == false || w->in2->set == false) {
                            error("AND Unset wire being used");
                        }
                        in1_and      = (block*)realloc(in1_and, block_num * (num_ands + 1) * sizeof(block));
                        in2_and      = (block*)realloc(in2_and, block_num * (num_ands + 1) * sizeof(block));
                        bool_in1_and = (bool*)realloc(bool_in1_and, bool_num * (num_ands + 1) * sizeof(bool));
                        bool_in2_and = (bool*)realloc(bool_in2_and, bool_num * (num_ands + 1) * sizeof(bool));
                        memcpy(in1_and + num_ands * block_num, w->in1->value, block_num * sizeof(block));
                        memcpy(in2_and + num_ands * block_num, w->in2->value, block_num * sizeof(block));
                        memcpy(bool_in1_and + num_ands * bool_num, w->in1->rem_value, bool_num * sizeof(bool));
                        memcpy(bool_in2_and + num_ands * bool_num, w->in2->rem_value, bool_num * sizeof(bool));
                        num_ands++;
                        if ((w->set == true) && (w->num_required > 0)) {
                            std::cout << *it << " " << w->num_required << std::endl;
                            error("Set wire is being set again");
                        }
                        break;
                    default:
                        break;
                }
            }
#ifdef __debug
            std::cout << "level " << i << ": before AND " << time_from(t) << "\t";
            t = clock_start();
#endif
            if (num_ands != 0) {
                out_and      = (block*)malloc(block_num * num_ands * sizeof(block));
                bool_out_and = (bool*)malloc(bool_num * num_ands * sizeof(bool));
                simd_circ->and_gate(bool_out_and, bool_in1_and, bool_in2_and, num_ands * bool_num, out_and, in1_and,
                                    in2_and, num_ands * block_num);
            }
#ifdef __debug
            std::cout << "AND: " << time_from(t) << "\t";
            t = clock_start();
#endif
            int num_and_rev = 0;
            for (auto it = level_map[i].begin(); it != level_map[i].end(); ++it) {
                Wire* w = circuit[*it];
                switch (w->type) {
                    case AND:
                        w->initialise_value(block_num, bool_num);
                        w->set_value(out_and + (num_and_rev * block_num), block_num,
                                     bool_out_and + (num_and_rev * bool_num), bool_num);
                        num_and_rev++;
                        w->in1->reset_value();
                        w->in2->reset_value();
                        break;
                    case XOR:
                        w->initialise_value(block_num, bool_num);
                        if (w->in1->set == false || w->in2->set == false)
                            error("XOR Unset wire being used");

                        if ((w->set == true) && (w->num_required > 0)) {
                            std::cout << *it << " " << w->num_required << std::endl;
                            error("Set wire is being set again");
                        }
                        simd_circ->xor_gate((w->value), (w->in1->value), (w->in2->value), block_num);
                        simd_circ->xor_gate(w->rem_value, w->in1->rem_value, w->in2->rem_value, bool_num);
                        w->set = true;
                        w->in1->reset_value();
                        w->in2->reset_value();
                        break;
                    case INV:
                        w->initialise_value(block_num, bool_num);
                        if (w->in1->set == false)
                            error("INV Unset wire being used");

                        if ((w->set == true) && (w->num_required > 0)) {
                            std::cout << *it << " " << w->num_required << std::endl;
                            error("Set wire is being set again");
                        }
                        simd_circ->not_gate((w->value), (w->in1->value), block_num);
                        simd_circ->not_gate(w->rem_value, w->in1->rem_value, bool_num);
                        w->set = true;
                        w->in1->reset_value();
                        break;
                    default:
                        break;
                }
            }
#ifdef __debug
            std::cout << "after AND: " << time_from(t) << "\n";
#endif
        }
#ifdef __debug
        t = clock_start();
#endif
        for (uint i = 0; i < n3; ++i) {
            Wire* w = circuit[num_wires - n3 + i];
            for (uint circ_index = 0; circ_index < block_num; ++circ_index) {
                bool dummy[128] = {0};
                block_to_bool(dummy, w->value[circ_index]);
                for (uint j = 0; j < 128; ++j) {
                    if ((circ_index * 128 + j) < num) {
                        out[i + (circ_index * 128 + j) * n3] = dummy[j];
                    }
                }
            }
            for (uint circ_index = 0; circ_index < bool_num; ++circ_index)
                out[i + (circ_index + block_num * 128) * n3] = w->rem_value[circ_index];
        }
#ifdef __debug
        std::cout << "output:" << time_from(t) << " us\n";
#endif
    }
};
}  // namespace emp