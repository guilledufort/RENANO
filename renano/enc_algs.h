// MIT License

// Copyright (c) 2021 Guillermo Dufort y Álvarez

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#ifndef ENC_ALGS_H
#define ENC_ALGS_H

#include <omp.h>

#include <iostream>

#include "compressor.h"

typedef struct
{
    int aligned = 0;
    int decompress = 0;

    file_io_t in_fq, in_paf;

    uint AVG_CANT, B_CTX, B_CTX_LEN, NS_MODEL_SIZE;

    uint16_t *ctx_avgs_sums;
    uint16_t *ctx_avgs_err_sums;
    uint32_t *ctx_err_avgs_total;

    uint32_t total_reads = 0;
    uint32_t enc_reads = 0;
} global_structs_t;

void init_global_structs(enano_params *p, context_models *cm, global_structs_t *gs);
void init_global_structs(enano_params *p, global_structs_t *gs);
void init_global_stats(context_models *cm, global_structs_t *gs);
void delete_global_stats(context_models *cm, global_structs_t *gs);
void copy_average_stats(Compressor *c, global_structs_t *gs);
void update_stats(context_models *cm, global_structs_t *gs, Compressor **comps, uc blocks_loaded);

bool load_data(global_structs_t *gs, enano_params *p, Compressor **comps, int update_load,
               uint &blocks_loaded);
bool load_data_decode(int in_fd, Compressor **comps, int update_load,
                      uint &blocks_loaded);

/* Parallelized Encode/Decode algorithms*/
int encode(enano_params &p);
int decode(enano_params &p);

/* Single thread Encode/Decode algorithms*/
int encode_st(enano_params &p);
int decode_st(enano_params &p);

#endif