#ifndef IO_H
#define IO_H

#include <errno.h>
#include <sys/types.h>
#include <unistd.h>

typedef struct {
    int fd = 0;
    char* name;
    char* buf;
    char* crt_ptr;
    uint32_t in_len = 0;
    uint32_t buf_size = 0;
} file_io_t;

static inline void init_file_io(file_io_t* f_io, char* name, int fd, uint32_t sz) {
    f_io->fd = fd;
    f_io->name = name;
    f_io->buf = new char[sz];
    f_io->crt_ptr = f_io->buf;
    f_io->buf_size = sz;
}

static inline void reset_io(file_io_t* f_io) {
    f_io->in_len = 0;
    f_io->crt_ptr = f_io->buf;
    close(f_io->fd);
    f_io->fd = open(f_io->name, O_RDONLY);
}

static inline void close_io(file_io_t* f_io) {
    delete[] f_io->buf;
}

static inline ssize_t xread(int fd, char* buf, size_t count) {
    ssize_t len, tlen;
    //    printf("Read %d \n", read_num++);
    tlen = 0;
    do {
        len = read(fd, buf, count);
        if (len == -1) {
            if (errno == EINTR) continue;
            return -1;
        }

        if (len == 0) return tlen;

        buf += len;
        count -= len;
        tlen += len;
    } while (count);

    return tlen;
}

/* Move remainder_length bytes starting at in_end pointer to the begining of f_io buffer*/
static inline void move_remainder_io(file_io_t* f_io, char* in_end) {
    uint32_t r = (in_end == f_io->buf + f_io->in_len) ? 0 : f_io->in_len - (in_end - f_io->buf) + 1;
    assert(r < f_io->buf_size);
    if (r > 0) {
        memmove(f_io->buf, in_end, r);
        f_io->in_len = r - 1;
    } else
        f_io->in_len = 0;
    f_io->crt_ptr = f_io->buf;
}

static inline ssize_t fill_buf_io(file_io_t* f_io) {
    ssize_t len, tlen;
    //    printf("Read %d \n", read_num++);
    char* buf = &f_io->buf[f_io->in_len];
    size_t count = f_io->buf_size - f_io->in_len;
    tlen = 0;
    do {
        len = read(f_io->fd, buf, count);
        if (len == -1) {
            if (errno == EINTR) continue;
            f_io->in_len = f_io->in_len - 1;
            return -1;
        }

        if (len == 0) break;

        buf += len;
        count -= len;
        tlen += len;
    } while (count);

    f_io->in_len += tlen;
    f_io->crt_ptr = f_io->buf;

    return tlen;
}

static inline ssize_t fill_buf_io(file_io_t* f_io, char* in_end) {
    move_remainder_io(f_io, in_end);
    ssize_t len, tlen;
    //    printf("Read %d \n", read_num++);
    char* buf = &f_io->buf[f_io->in_len];
    size_t count = f_io->buf_size - f_io->in_len;
    tlen = 0;
    do {
        len = read(f_io->fd, buf, count);
        if (len == -1) {
            if (errno == EINTR) continue;
            f_io->in_len = f_io->in_len - 1;
            return -1;
        }

        if (len == 0) break;

        buf += len;
        count -= len;
        tlen += len;
    } while (count);

    f_io->in_len += tlen;
    f_io->crt_ptr = f_io->buf;

    return tlen;
}

#endif