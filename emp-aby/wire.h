#pragma once

#include "emp-aby/simd_interface/simd_exec.h"
#include <fstream>
namespace emp {

#define INPUT 0
#define AND   1
#define XOR   2
#define INV   3
//#define __debug__
class Wire {
public:
    int type;
    block* value;
    bool* rem_value;
    uint32_t level;
    Wire* in1        = nullptr;
    Wire* in2        = nullptr;
    bool set         = false;
    int num_required = 0;
    int num_used     = 0;

    Wire(int type) {
        if (type == INPUT) {
            this->type  = type;
            this->level = 0;
        }
        else {
            this->type = type;
        }
    }

    Wire(int type, Wire* in1) {
        if (type == INV) {
            this->type              = type;
            this->in1               = in1;
            this->level             = in1->level;
            this->in1->num_required = this->in1->num_required + 1;
        }
        else if (type == INPUT)
            error("Input wire has no imcoming!");
        else
            error("AND and XOR require 2 inputs!");
    }

    Wire(int type, Wire* in1, Wire* in2) {
        this->type              = type;
        this->in1               = in1;
        this->in2               = in2;
        this->in1->num_required = this->in1->num_required + 1;
        this->in2->num_required = this->in2->num_required + 1;
        switch (this->type) {
            case INPUT:
            case INV:
                error("Cannot set second incoming wire!");
            case AND:
                this->level = std::max(this->in1->level, this->in2->level) + 1;
                break;
            case XOR:
                this->level = std::max(this->in1->level, this->in2->level);
                break;
        }
    }
    void initialise_value(int n, int m = 0) {
        if (this->set == true && this->num_required > 0)
            error("Wire being initialised twice!");
        this->value = (block*)malloc(n * sizeof(block));
        memset(this->value, 0, n * sizeof(block));
        this->rem_value = (bool*)malloc(m * sizeof(bool));
        memset(this->rem_value, 0, m * sizeof(bool));
    }
    void set_value(block* value, int n, bool* rem_value, int m) {
        memcpy(this->value, value, n * sizeof(block));
        memcpy(this->rem_value, rem_value, m * sizeof(bool));
        if (this->num_required > 0)
            this->set = true;
    }
    void reset_value() {
        this->num_used = this->num_used + 1;
        if (this->num_used >= this->num_required) {
            free(this->value);
            this->set      = false;
            this->num_used = 0;
            free(this->rem_value);
        }
    }
};

}  // namespace emp