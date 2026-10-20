#pragma once

#include "emp-aby/utils.h"

namespace emp {
template <typename T>
class MPIOChannel {
public:
    virtual void send_data(int dst, const void* data, int64_t len)       = 0;
    virtual void recv_data(int src, void* data, int64_t len)             = 0;
    virtual void send_bool(int dst, const bool* data, int64_t length)    = 0;
    virtual void recv_bool(int src, bool* data, int64_t length)          = 0;
    virtual void send_block(int dst, const block* data, int64_t length)  = 0;
    virtual void recv_block(int src, block* data, int64_t length)        = 0;
    virtual void sync()                                                  = 0;
    virtual void flush(int idx = 0)                                      = 0;
    virtual T*& get(size_t idx, bool b = false)                          = 0;
    virtual ~MPIOChannel()                                               = 0;
    virtual int get_total_bytes_sent()                                   = 0;
};

template <typename IO>
MPIOChannel<IO>::~MPIOChannel() {}
}  // namespace emp
