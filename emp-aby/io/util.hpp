#pragma once

#include <map>
#include "emp-tool/emp-tool.h"

namespace emp {
void accept_base_connections(int expected_connections, int bind_port, string bind_address,
                             std::map<uint, int>& socket_map, int party, int num_party) {
    int reuse = 1;
    struct sockaddr_in dest;
    struct sockaddr_in serv;

    memset(&serv, 0, sizeof(serv));
    serv.sin_family      = AF_INET;
    serv.sin_addr.s_addr = inet_addr(bind_address.c_str()); /* set our address to any interface */
    serv.sin_port        = htons(bind_port);                /* set the server port number */
    socklen_t socksize   = sizeof(struct sockaddr_in);
    int mysocket         = socket(AF_INET, SOCK_STREAM, 0);

    setsockopt(mysocket, SOL_SOCKET, SO_REUSEADDR, (const char*)&reuse, sizeof(reuse));
    if (bind(mysocket, (struct sockaddr*)&serv, sizeof(struct sockaddr)) < 0) {
        error("error: bind");
    }

    if (listen(mysocket, expected_connections) < 0) {
        error("error: listen");
    }
    for (int i = 0; i < expected_connections; ++i) {
        int consocket = accept(mysocket, (struct sockaddr*)&dest, &socksize);

        const int one = 1;
        setsockopt(consocket, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
        int p;
        read(consocket, &p, sizeof(p));
        if (p > num_party) {
            error("connected to wrong party");
        }
        if (socket_map.find(p) != socket_map.end()) {
            error("already connected to this party");
        }
        if (send(consocket, &party, sizeof(party), 0) < 0) {
            error("error: send");
        }
        socket_map.emplace(std::make_pair(p, consocket));
    }
}

int request_base_connection(int p, std::vector<std::pair<std::string, unsigned short>>& net_config, int party,
                            int num_party) {
    struct sockaddr_in dest;
    memset(&dest, 0, sizeof(dest));
    dest.sin_family      = AF_INET;
    dest.sin_addr.s_addr = inet_addr(net_config[p - 1].first.c_str());
    dest.sin_port        = htons(net_config[p - 1].second);

    int consocket;

    while (1) {
        consocket = socket(AF_INET, SOCK_STREAM, 0);

        if (connect(consocket, (struct sockaddr*)&dest, sizeof(struct sockaddr)) == 0) {
            break;
        }

        close(consocket);
        usleep(1000);
    }
    const int one = 1;
    setsockopt(consocket, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
    send(consocket, &party, sizeof(party), 0);
    int new_p;
    read(consocket, &new_p, sizeof(new_p));
    if (new_p != p) {
        std::cout << party << ":  " << new_p << " != " << p << std::endl;
        error("connected to wrong party");
    }

    return consocket;
}

}  // namespace emp