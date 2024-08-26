#pragma once

#include <emp-tool/emp-tool.h>

namespace emp {
enum MESSAGE_TYPE : uint8_t { NORM_MSG = 0, BOOT_REQ_MSG = 1, BOOT_RSP_MSG = 2, TERMINATE_MSG = 3 };
template <typename T>
class MPIOChannel {
public:
    virtual void send_data(int dst, const void* data, int len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG) = 0;
    virtual void recv_data(int src, void* data, int len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG)       = 0;
    virtual void* recv_data(int src, int& len, int j = 0, MESSAGE_TYPE msg_type = NORM_MSG)                 = 0;
    virtual void send_bool(int dst, bool* data, int length, int j = 0)                                      = 0;
    virtual void recv_bool(int src, bool* data, int length, int j = 0)                                      = 0;
    virtual void send_block(int dst, const block* data, int length, int j = 0)                              = 0;
    virtual void recv_block(int src, block* data, int length, int j = 0)                                    = 0;
    virtual void sync()                                                                                     = 0;
    virtual void flush(int idx = 0, int j = 0)                                                              = 0;
    virtual T*& get(size_t idx, bool b = false)                                                             = 0;
    virtual ~MPIOChannel()                                                                                  = 0;
    virtual int get_total_bytes_sent()                                                                      = 0;
};

template <typename IO>
MPIOChannel<IO>::~MPIOChannel() {}
}  // namespace emp
