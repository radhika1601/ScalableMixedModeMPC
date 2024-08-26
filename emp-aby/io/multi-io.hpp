#pragma once

#include "emp-aby/io/mp_io_channel.h"
#include "emp-aby/io/multi-io-base.hpp"
#include <poll.h>

namespace emp {
class MultiIO : public MPIOChannel<MultiIOBase> {
private:
    /* data */
public:
    int party, num_party;
    int bind_port;
    string bind_address;
    std::vector<std::pair<std::string, unsigned short>> net_config;
    std::map<uint, MultiIOBase*> ios;
    std::map<uint, MultiIOBase*> ot_ios[2];
    bool continue_comm;
    std::future<void> background_recv_fut;
    MultiIO(int party, int num_party, std::vector<std::pair<std::string, unsigned short>>& net_config);
    ~MultiIO();

    void setup_ot_ios();
    void send_data(int dst, const void* data, int len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG);
    void recv_data(int src, void* data, int len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG);
    void* recv_data(int src, int& len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG);

    int get_total_bytes_sent();

    void flush(int idx = 0, int j = 0) {}
    void sync() {}

    // optimise later to send aligned data
    void send_bool(int dst, bool* data, int length, int j = 0) {
        send_data(dst, data, length, j, NORM_MSG);
    }

    // optimise later to recv aligned data
    void recv_bool(int src, bool* data, int length, int j = 0) {
        recv_data(src, data, length, j, NORM_MSG);
    }

    void send_block(int dst, const block* data, int length, int j = 0) {
        send_data(dst, data, length * sizeof(block), j, NORM_MSG);
    }

    void recv_block(int src, block* data, int length, int j = 0) {
        recv_data(src, data, length * sizeof(block), j, NORM_MSG);
    }

    MultiIOBase*& get(size_t idx, bool b = false) {
        if (b)
            return ot_ios[0][idx];
        else
            return ot_ios[1][idx];
    }

    void background_recv();
};

MultiIO::MultiIO(int party, int num_party, std::vector<std::pair<std::string, unsigned short>>& net_config)
    : party(party),
      num_party(num_party),
      bind_port(net_config[party - 1].second),
      bind_address(net_config[party - 1].first),
      net_config(net_config) {
    std::map<uint, int> socket_map;
    accept_base_connections(num_party - party, net_config[party - 1].second, net_config[party - 1].first, socket_map,
                            party, num_party);
    for (uint p = party - 1; p > 0; --p) {
        int consocket = request_base_connection(p, net_config, party, num_party);
        socket_map.emplace(std::pair(p, consocket));
    }
    for (auto& sock : socket_map) {
        MultiIOBase* io = new MultiIOBase(sock.second, true);
        ios.emplace(std::pair(sock.first, io));
    }
    socket_map.clear();
    continue_comm       = true;
    background_recv_fut = std::async([this]() {
        this->background_recv();
    });
}

void MultiIO::setup_ot_ios() {
    std::map<uint, int> socket_map;

    for (auto& conf : net_config) {
        conf.second++;
    }
    accept_base_connections(num_party - party, net_config[party - 1].second, net_config[party - 1].first, socket_map,
                            party, num_party);
    for (uint p = party - 1; p > 0; --p) {
        int consocket = request_base_connection(p, net_config, party, num_party);
        socket_map.emplace(std::pair(p, consocket));
    }
    for (auto& sock : socket_map) {
        MultiIOBase* io = new MultiIOBase(sock.second, true);
        ot_ios[0].emplace(std::pair(sock.first, io));
    }
    socket_map.clear();
    for (auto& conf : net_config) {
        conf.second++;
    }
    accept_base_connections(num_party - party, net_config[party - 1].second, net_config[party - 1].first, socket_map,
                            party, num_party);
    for (uint p = party - 1; p > 0; --p) {
        int consocket = request_base_connection(p, net_config, party, num_party);
        socket_map.emplace(std::pair(p, consocket));
    }
    for (auto& sock : socket_map) {
        MultiIOBase* io = new MultiIOBase(sock.second, true);
        ot_ios[1].emplace(std::pair(sock.first, io));
    }
    socket_map.clear();
}

MultiIO::~MultiIO() {
    for (auto& io : ios) {
        io.second->send_msg(nullptr, 0, TERMINATE_MSG);
    }
    bool c = false;
    for (auto& io : ios) {
        c |= io.second->continue_comm;
    }
    continue_comm = c;
    background_recv_fut.get();
    for (auto& io : ios) {
        io.second->~MultiIOBase();
    }
    for (auto& io : ot_ios[0]) {
        io.second->~MultiIOBase();
    }
    for (auto& io : ot_ios[1]) {
        io.second->~MultiIOBase();
    }
    net_config.clear();
    ios.clear();
}

void MultiIO::send_data(int dst, const void* data, int len, int j, MESSAGE_TYPE msg_type) {
    if (dst != 0 && dst != party) {
        MultiIOBase* io = ios[dst];
        io->send_msg(data, len, msg_type);
    }
    else {
        error("sending to invalid party");
    }
}

void MultiIO::recv_data(int src, void* data, int len, int j, MESSAGE_TYPE msg_type) {
    if (src != 0 && src != party) {
        MultiIOBase* io = ios[src];
        bool received   = false;
        while (!received) {
            std::unique_lock lock(io->recv_mutex[msg_type]);
            if (!io->recv_msg_queue[msg_type].empty()) {
                received     = true;
                int recv_len = io->recv_msg_queue[msg_type].front().first;
                if (len != recv_len) {
                    std::cout << "lengths" << len << " " << recv_len << "\n";
                    error("unequal length");
                }
                memcpy(data, io->recv_msg_queue[msg_type].front().second, len);
                io->recv_msg_queue[msg_type].pop_front();
                lock.unlock();
                return;
            }
            else {
                io->recv_condition_vars[msg_type].wait(
                    lock, [io, msg_type] { return !io->recv_msg_queue[msg_type].empty(); });
            }
        }
    }
    else {
        error("receive called for invalid party");
    }
}

void* MultiIO::recv_data(int src, int& len, int j, MESSAGE_TYPE msg_type) {
    void* data = nullptr;
    if (src != 0 and src != party) {
        MultiIOBase* io = ios[src];
        bool received   = false;
        while (!received) {
            std::unique_lock lock(io->recv_mutex[msg_type]);
            if (!io->recv_msg_queue[msg_type].empty()) {
                received     = true;
                int recv_len = io->recv_msg_queue[msg_type].front().first;
                len          = recv_len;
                data         = io->recv_msg_queue[msg_type].front().second;
                io->recv_msg_queue[msg_type].pop_front();
                lock.unlock();
                return data;
            }
            else {
                io->recv_condition_vars[msg_type].wait(
                    lock, [io, msg_type] { return !io->recv_msg_queue[msg_type].empty(); });
            }
        }
    }
    else {
        error("receive called for invalid party");
    }
    return data;
}

int MultiIO::get_total_bytes_sent() {
    double kb = 0;
    for (auto& io : this->ios) {
        kb += (double)(io.second->counter) / 1000;
    }
    for (auto& io : this->ot_ios[0]) {
        kb += (double)(io.second->counter) / 1000;
    }
    for (auto& io : this->ot_ios[1]) {
        kb += (double)(io.second->counter) / 1000;
    }
    return kb;
}

void MultiIO::background_recv() {
    struct pollfd* pfds;
    pfds = (struct pollfd*)calloc(num_party - 1, sizeof(struct pollfd));
    std::map<int, uint> socket_party;
    nfds_t i = 0;
    for (auto& io : ios) {
        pfds[i].fd     = io.second->consocket;
        pfds[i].events = POLLIN;
        socket_party.emplace(io.second->consocket, io.first);
        ++i;
    }

    while (continue_comm) {
        int ready = 0;
        while (ready < 1) {
            if (!continue_comm)
                return;
            ready = poll(pfds, num_party - 1, -1);
            if (ready == -1) {
                error("error: poll");
            }
        }
        bool c = true;
        for (nfds_t i = 0; i < num_party - 1; ++i) {
            if (pfds[i].revents != 0) {
                if (pfds[i].revents & POLLIN) {
                    int p = socket_party[pfds[i].fd];
                    c     = ios[p]->recv_msg();
                }
            }
        }
        if (!c) {
            for (auto& io : ios) {
                c |= io.second->continue_comm;
            }
            continue_comm = c;
        }
    }
}

}  // namespace emp
