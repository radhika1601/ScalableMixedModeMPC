#pragma once

#include "emp-aby/io/util.hpp"

namespace emp {
// Blocking TCP channel to one peer.
class MultiIOBase : public IOChannel {
public:
    int consocket = -1;

    MultiIOBase(int consocket, bool quiet = false) : consocket(consocket) {
        set_nodelay();
        if (!quiet)
            std::cout << "connected\n";
    }

    ~MultiIOBase() {
        close(consocket);
    }

    void set_nodelay() {
        const int one = 1;
        setsockopt(consocket, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
    }

    void set_delay() {
        const int zero = 0;
        setsockopt(consocket, IPPROTO_TCP, TCP_NODELAY, &zero, sizeof(zero));
    }

    void send_data_internal(const void* data, int64_t len) override {
        int64_t sent = 0;
        while (sent < len) {
            ssize_t res = write(consocket, (char*)data + sent, len - sent);
            if (res > 0)
                sent += res;
            else
                error("net_send_data\n");
        }
    }

    void recv_data_internal(void* data, int64_t len) override {
        int64_t recvd = 0;
        while (recvd < len) {
            ssize_t res = read(consocket, (char*)data + recvd, len - recvd);
            if (res > 0)
                recvd += res;
            else
                error("net_recv_data\n");
        }
    }
};

}  // namespace emp
