#pragma once

#include <shared_mutex>
#include <sys/select.h>
#include "emp-aby/io/util.hpp"

namespace emp {
class MultiIOBase : public IOChannel<MultiIOBase> {
private:
public:
    int consocket = -1;
    std::shared_mutex sock_mutex;
    bool continue_comm = false;
    std::deque<std::pair<int, void*>> recv_msg_queue[3];
    std::condition_variable recv_condition_vars[3];
    std::mutex recv_mutex[3];

    MultiIOBase(int consocket, bool quiet = false) : consocket(consocket), continue_comm(true) {
        set_nodelay();
        if (!quiet)
            std::cout << "connected\n";
    }

    void sync() {}
    void flush() {}

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

    void send_msg(const void* data, int len, MESSAGE_TYPE msg_type = NORM_MSG) {
        char* meta_buff = (char*)malloc(5);
        meta_buff[0]    = msg_type;
        memcpy(meta_buff + 1, &(len), 4);
        std::shared_lock lock(sock_mutex);
        this->send_data(meta_buff, 5);
        if(msg_type != TERMINATE_MSG)
            this->send_data(data, len);
        lock.unlock();
        free(meta_buff);
    }

    bool recv_msg() {
        char* meta_buff = (char*)malloc(5);
        std::shared_lock lock(sock_mutex);

        // add select here
        // select(consocket, );

        this->recv_data(meta_buff, 5);
        MESSAGE_TYPE recv_type = static_cast<MESSAGE_TYPE>(meta_buff[0]);
        if (recv_type == TERMINATE_MSG) {
            lock.unlock();
            this->continue_comm = false;
            return false;
        }

        int len;
        memcpy(&(len), meta_buff + 1, 4);
        // std::cout << "received length" << len << std::endl;
        void* data = malloc(len);
        this->recv_data(data, len);
        lock.unlock();

        std::unique_lock que_lock(recv_mutex[recv_type]);
        recv_msg_queue[recv_type].push_back(std::pair<int, void*>(len, data));
        que_lock.unlock();
        recv_condition_vars[recv_type].notify_one();
        return true;
    }

    void send_data_internal(const void* data, size_t len) {
        size_t sent = 0;
        while (sent < len) {
            size_t res = write(consocket, (char*)data + sent, len - sent);
            if (res > 0)
                sent += res;
            else
                error("net_send_data\n");
        }
    }

    void recv_data_internal(void* data, size_t len) {
        size_t recvd = 0;
        while (recvd < len) {
            size_t res = read(consocket, (char*)data + recvd, len - recvd);
            if (res > 0)
                recvd += res;
            else
                error("net_recv_data\n");
        }
    }
};

}  // namespace emp
