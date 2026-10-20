#include "emp-aby/simd_interface/simd.h"
#include <map>
using namespace emp;

// Checks that and_helper never hands out the same triple twice. The fake provider
// tags every triple it generates with a unique counter, so a reused triple shows up
// as a repeated tag. Runs locally, with no network.

block tag(int64_t i) {
    return makeBlock(0, i);
}
int64_t untag(block b) {
    return _mm_extract_epi64(b, 0);
}
int64_t untag(int64_t i) {
    return i;
}

struct FakeBtp {
    int64_t next = 0;
    size_t buffer_sz;

    void get_triple(int64_t* a, int64_t* b, int64_t* c) {
        for (size_t i = 0; i < buffer_sz; ++i, ++next)
            a[i] = b[i] = c[i] = next;
    }
    void get_triple(block* a, block* b, block* c) {
        for (size_t i = 0; i < buffer_sz; ++i, ++next)
            a[i] = b[i] = c[i] = tag(next);
    }
};

class FakeExec : public SIMDCircuitExecution<FakeBtp> {
public:
    FakeBtp btp;
    FakeBtp* getBtp() override {
        return &btp;
    }
    using SIMDCircuitExecution<FakeBtp>::and_helper;

    void and_gate(block*, block*, block*, size_t) override {}
    void xor_gate(block*, block*, block*, size_t) override {}
    void not_gate(block*, block*, size_t) override {}
    void and_gate(bool*, bool*, bool*, size_t, block*, block*, block*, size_t) override {}
    void and_gate(bool*, bool*, bool*, size_t) override {}
    void xor_gate(bool*, bool*, bool*, size_t) override {}
    void not_gate(bool*, bool*, size_t) override {}
};

template <typename T>
void test_no_reuse(const vector<size_t>& lengths, const char* name) {
    const size_t pool_sz = 100;
    FakeExec exec;
    exec.btp.buffer_sz = pool_sz;
    T *pool_a = new T[pool_sz], *pool_b = new T[pool_sz], *pool_c = new T[pool_sz];
    exec.btp.get_triple(pool_a, pool_b, pool_c);
    size_t num_triples = 0;

    std::map<int64_t, int> seen;
    for (size_t length : lengths) {
        T *a, *b, *c;
        bool delete_array = false;
        exec.and_helper<T>(a, b, c, length, delete_array, pool_a, pool_b, pool_c, pool_sz, num_triples);
        for (size_t i = 0; i < length; ++i) {
            int64_t id = untag(a[i]);
            if (id != untag(b[i]) || id != untag(c[i]) || id < 0 || id >= exec.btp.next)
                error("and_helper returned a triple the provider never generated");
            if (++seen[id] > 1) {
                std::cout << name << ": triple " << id << " reused\n";
                error("and_helper reused a triple");
            }
        }
        if (delete_array) {
            delete[] a;
            delete[] b;
            delete[] c;
        }
    }
    delete[] pool_a;
    delete[] pool_b;
    delete[] pool_c;
    std::cout << name << ": " << seen.size() << " triples, no reuse\n";
}

int main() {
    // Mixes batches below the pool, across a refill, and larger than the pool,
    // with the large batches arriving both on a fresh and a partly used pool.
    vector<vector<size_t>> cases = {
        {30, 250, 60, 50, 90},
        {250, 10, 250, 99, 1, 100},
        {99, 101, 99, 301, 7},
        {0, 100, 100, 200, 3, 97, 500},
    };
    for (auto& lengths : cases) {
        test_no_reuse<int64_t>(lengths, "int64_t");
        test_no_reuse<block>(lengths, "block");
    }
    return 0;
}
