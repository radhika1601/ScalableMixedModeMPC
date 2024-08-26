#pragma once

#include "emp-aby/he_interface.hpp"

namespace emp {

template <typename IO>
class ArithmeticCirc {
private:
    PRG prg;
    MPIOChannel<IO>* io;
    HE<IO>* he;
    int num_party, party;
    int64_t *triple_a, *triple_b, *triple_c;

public:
    size_t num_triples_pool, num_triples = 0;
    void get_triples(int64_t* triple_a, int64_t* triple_b, int64_t* triple_c);
    ArithmeticCirc(int num_party, int party, MPIOChannel<IO>* io, HE<IO>* he);
    ~ArithmeticCirc();
    void sum(int64_t* out, int64_t* in1, int64_t* in2, size_t length);
    void sub(int64_t* out, int64_t* in1, int64_t* in2, size_t length);
    void mult(int64_t* out, int64_t* in1, int64_t* in2, size_t length);
};

template <typename IO>
ArithmeticCirc<IO>::ArithmeticCirc(int num_party, int party, MPIOChannel<IO>* io, HE<IO>* he) {
    this->num_party        = num_party;
    this->party            = party;
    this->io               = io;
    this->he               = he;
    this->num_triples_pool = 20 * (he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2);

    triple_a = new int64_t[num_triples_pool];
    triple_b = new int64_t[num_triples_pool];
    triple_c = new int64_t[num_triples_pool];

    this->get_triples(triple_a, triple_b, triple_c);
}

template <typename IO>
ArithmeticCirc<IO>::~ArithmeticCirc() {
    delete[] triple_a;
    delete[] triple_b;
    delete[] triple_c;
}

template <typename IO>
void ArithmeticCirc<IO>::sum(int64_t* out, int64_t* in1, int64_t* in2, size_t length) {
    for (int i = 0; i < length; ++i)
        out[i] = (in1[i] + in2[i]) % he->q;
}

template <typename IO>
void ArithmeticCirc<IO>::sub(int64_t* out, int64_t* in1, int64_t* in2, size_t length) {
    for (int i = 0; i < length; ++i)
        out[i] = (in1[i] - in2[i]) % he->q;
}

template <typename IO>
void ArithmeticCirc<IO>::get_triples(int64_t* triple_a, int64_t* triple_b, int64_t* triple_c) {
    int batch_size = (he->cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2);
    prg.random_data(triple_a, num_triples_pool * sizeof(int64_t));
    prg.random_data(triple_b, num_triples_pool * sizeof(int64_t));
    for (size_t i = 0; i < num_triples_pool; ++i) {
        triple_a[i] %= he->q;
        triple_b[i] %= he->q;
        triple_a[i] = (he->q + triple_a[i]) % he->q;
        triple_b[i] = (he->q + triple_b[i]) % he->q;
    }
    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> a, b, c;
    for (int i = 0; i < num_triples_pool / batch_size; ++i) {
        lbcrypto::Plaintext p_a, p_b;
        vector<int64_t> tmp_a, tmp_b;
        tmp_a.resize(batch_size);
        tmp_b.resize(batch_size);
        memcpy(tmp_a.data(), triple_a + i * batch_size, batch_size * sizeof(int64_t));
        memcpy(tmp_b.data(), triple_b + i * batch_size, batch_size * sizeof(int64_t));

        p_a = he->cc->MakePackedPlaintext(tmp_a);
        p_b = he->cc->MakePackedPlaintext(tmp_b);

        a.push_back(he->cc->Encrypt(he->pk, p_a));
        b.push_back(he->cc->Encrypt(he->pk, p_b));
    }
    if (party == ALICE) {
        std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> tmp_a, tmp_b;
        for (int i = 2; i <= num_party; ++i) {
            he->deserialize_recv(tmp_a, i);
            he->deserialize_recv(tmp_b, i);
            for (int j = 0; j < num_triples_pool / batch_size; ++j) {
                he->cc->EvalAddInPlace(a[j], tmp_a[j]);
                he->cc->EvalAddInPlace(b[j], tmp_b[j]);
                he->cc->ModReduceInPlace(a[j]);
                he->cc->ModReduceInPlace(b[j]);
            }
        }

        for (int i = 0; i < num_triples_pool / batch_size; ++i) {
            auto tmp = he->cc->EvalMult(a[i], b[i]);
            c.push_back(he->cc->ModReduce(tmp));
        }
        he->serialize_sendall(c);
    }
    else {
        he->serialize_send(a, ALICE);
        he->serialize_send(b, ALICE);
        he->deserialize_recv(c, ALICE);
    }
    he->enc_to_share(c, triple_c, num_triples_pool);
}

template <typename IO>
void ArithmeticCirc<IO>::mult(int64_t* out, int64_t* in1, int64_t* in2, size_t length) {
    bool delete_array = false;
    int64_t *a, *b, *c;
    // std::cout << "In mult \n";
    if (length > num_triples_pool) {
        a = new int64_t[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
        b = new int64_t[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
        c = new int64_t[(length + num_triples_pool - 1) / num_triples_pool * num_triples_pool];
        for (uint i = 0; i < (length + num_triples_pool - 1) / num_triples_pool; ++i)
            this->get_triples(a + i * num_triples_pool, b + i * num_triples_pool, c + i * num_triples_pool);
        size_t tocp = min((length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length, num_triples);
        memcpy(triple_a, a + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
               tocp * sizeof(int64_t));
        memcpy(triple_b, b + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
               tocp * sizeof(int64_t));
        memcpy(triple_c, c + (length + num_triples_pool - 1) / num_triples_pool * num_triples_pool - length,
               tocp * sizeof(int64_t));
        num_triples  = 0;
        delete_array = true;
    }
    else if (length > num_triples_pool - num_triples) {
        a            = new int64_t[length];
        b            = new int64_t[length];
        c            = new int64_t[length];
        delete_array = true;
        memcpy(a, triple_a + num_triples, (num_triples_pool - num_triples) * sizeof(int64_t));
        memcpy(c, triple_c + num_triples, (num_triples_pool - num_triples) * sizeof(int64_t));
        memcpy(b, triple_b + num_triples, (num_triples_pool - num_triples) * sizeof(int64_t));
        get_triples(triple_a, triple_b, triple_c);
        memcpy(a + num_triples_pool - num_triples, triple_a,
               (length - (num_triples_pool - num_triples)) * sizeof(int64_t));
        memcpy(b + num_triples_pool - num_triples, triple_b,
               (length - (num_triples_pool - num_triples)) * sizeof(int64_t));
        memcpy(c + num_triples_pool - num_triples, triple_c,
               (length - (num_triples_pool - num_triples)) * sizeof(int64_t));
        num_triples = length - (num_triples_pool - num_triples);
    }
    else {
        a = triple_a + num_triples;
        b = triple_b + num_triples;
        c = triple_c + num_triples;
        num_triples += length;
    }

    int64_t *d = new int64_t[length], *e = new int64_t[length];

    for (int i = 0; i < length; ++i) {
        d[i] = (he->q + in1[i] - a[i]) % he->q;
        e[i] = (he->q + in2[i] - b[i]) % he->q;
        c[i] = (he->q + c[i]) % he->q;
    }

    // io->sync();
    if (party == ALICE) {
        int64_t *d0 = new int64_t[length], *e0 = new int64_t[length];

        for (int i = 2; i <= num_party; ++i) {
            io->recv_data(i, d0, length * sizeof(int64_t));
            io->recv_data(i, e0, length * sizeof(int64_t));
            for (int j = 0; j < length; ++j) {
                d[j] = (d[j] + d0[j]) % he->q;
                e[j] = (e[j] + e0[j]) % he->q;
            }
        }

        for (int i = 2; i <= num_party; ++i) {
            io->send_data(i, d, length * sizeof(int64_t));
            io->send_data(i, e, length * sizeof(int64_t));
            io->flush(i);
        }

        delete[] d0;
        delete[] e0;
    }
    else {
        io->send_data(ALICE, d, length * sizeof(int64_t));
        io->send_data(ALICE, e, length * sizeof(int64_t));
        io->flush(ALICE);
        io->recv_data(ALICE, d, length * sizeof(int64_t));
        io->recv_data(ALICE, e, length * sizeof(int64_t));
    }
    io->flush();

    for (uint i = 0; i < length; ++i) {
        long long int x = ((uint64_t)((uint64_t)e[i] * (uint64_t)a[i])) % he->q;
        long long int y = ((uint64_t)((uint64_t)d[i] * (uint64_t)b[i])) % he->q;

        out[i] = (c[i] + x + y) % he->q;
    }
    if (party == ALICE) {
        for (uint i = 0; i < length; ++i) {
            long long int x = ((uint64_t)((uint64_t)d[i] * (uint64_t)e[i])) % he->q;
            out[i]          = (out[i] + x + he->q) % he->q;
        }
    }
    delete[] d;
    delete[] e;
    if (delete_array) {
        delete[] a;
        delete[] b;
        delete[] c;
    }
}
}  // namespace emp
