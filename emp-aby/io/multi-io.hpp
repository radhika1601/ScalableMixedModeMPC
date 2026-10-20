#pragma once

#include "emp-aby/io/mp_io_channel.h"
#include "emp-aby/io/multi-io-base.hpp"

namespace emp {
// One blocking channel per peer for protocol messages, plus two per peer for OT.
// Sends and receives go straight to the socket. A send blocks once the socket
// buffer is full, so two parties must not both send large data before reading.
class MultiIO : public MPIOChannel<MultiIOBase> {
public:
    int party, num_party;
    int bind_port;
    string bind_address;
    std::vector<std::pair<std::string, unsigned short>> net_config;
    std::map<uint, MultiIOBase*> ios;
    std::map<uint, MultiIOBase*> ot_ios[2];
    MultiIO(int party, int num_party, std::vector<std::pair<std::string, unsigned short>>& net_config);
    ~MultiIO();

    void setup_ot_ios();

    void send_data(int dst, const void* data, int64_t len) {
        peer(dst)->send_data(data, len);
    }

    void recv_data(int src, void* data, int64_t len) {
        peer(src)->recv_data(data, len);
    }

    void send_bool(int dst, const bool* data, int64_t length) {
        peer(dst)->send_bool(data, length);
    }

    void recv_bool(int src, bool* data, int64_t length) {
        peer(src)->recv_bool(data, length);
    }

    void send_block(int dst, const block* data, int64_t length) {
        peer(dst)->send_block(data, length);
    }

    void recv_block(int src, block* data, int64_t length) {
        peer(src)->recv_block(data, length);
    }

    void flush(int idx = 0) {}
    void sync() {}

    int get_total_bytes_sent();

    MultiIOBase*& get(size_t idx, bool b = false) {
        if (b)
            return ot_ios[0][idx];
        else
            return ot_ios[1][idx];
    }

private:
    MultiIOBase* peer(int p) {
        if (p == 0 || p == party || ios.find(p) == ios.end())
            error("invalid party");
        return ios[p];
    }

    void connect_all(std::map<uint, MultiIOBase*>& out);
};

void MultiIO::connect_all(std::map<uint, MultiIOBase*>& out) {
    std::map<uint, int> socket_map;
    accept_base_connections(num_party - party, net_config[party - 1].second, net_config[party - 1].first, socket_map,
                            party, num_party);
    for (uint p = party - 1; p > 0; --p) {
        int consocket = request_base_connection(p, net_config, party, num_party);
        socket_map.emplace(std::pair(p, consocket));
    }
    for (auto& sock : socket_map) {
        out.emplace(std::pair(sock.first, new MultiIOBase(sock.second, true)));
    }
}

MultiIO::MultiIO(int party, int num_party, std::vector<std::pair<std::string, unsigned short>>& net_config)
    : party(party),
      num_party(num_party),
      bind_port(net_config[party - 1].second),
      bind_address(net_config[party - 1].first),
      net_config(net_config) {
    connect_all(ios);
}

void MultiIO::setup_ot_ios() {
    for (int i = 0; i < 2; ++i) {
        for (auto& conf : net_config) {
            conf.second++;
        }
        connect_all(ot_ios[i]);
    }
}

MultiIO::~MultiIO() {
    for (auto& io : ios) {
        delete io.second;
    }
    for (auto& io : ot_ios[0]) {
        delete io.second;
    }
    for (auto& io : ot_ios[1]) {
        delete io.second;
    }
}

int MultiIO::get_total_bytes_sent() {
    double kb = 0;
    for (auto& io : this->ios) {
        kb += (double)(io.second->send_counter) / 1000;
    }
    for (auto& io : this->ot_ios[0]) {
        kb += (double)(io.second->send_counter) / 1000;
    }
    for (auto& io : this->ot_ios[1]) {
        kb += (double)(io.second->send_counter) / 1000;
    }
    return kb;
}

}  // namespace emp
