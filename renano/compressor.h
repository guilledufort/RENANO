// MIT License

// Copyright (c) 2020 Guillermo Dufort y Álvarez

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <fcntl.h>
#include <inttypes.h>
#include <string.h>
#include <sys/types.h>
#include <omp.h>

#include "alignments.h"
#include "paf_process.h"
#include "params.h"
#include "io.h"

#ifndef WIN32

#include <unistd.h>

#else
#include "unistd.h"
#endif

#define NDEBUG
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>

typedef struct
{
    //Read lengths
    SIMPLE_MODEL<256> model_len1;
    SIMPLE_MODEL<256> model_len2;
    SIMPLE_MODEL<256> model_len3;
    SIMPLE_MODEL<2> model_same_len;

    // Names
    SIMPLE_MODEL<256> model_name_prefix[256];
    SIMPLE_MODEL<256> model_name_suffix[256];
    SIMPLE_MODEL<256> model_name_len[256];
    SIMPLE_MODEL<128> model_name_middle[8192];

    // Basecalls models;
    BASE_MODEL<uint8_t> *model_seq;
    ali_models *am;

    // Qualities
    SIMPLE_MODEL<QUANT_D_CANT> model_qual_quant[CTX_CNT];
    SIMPLE_MODEL<QMAX - QUANT_D_CANT> quant_top;
} context_models;

// Takes a context model and updates its accumulated frequencies of the models to prepare for fast encoding/decoding
void update_AccFreqs(enano_params *p, context_models *ctx_m, bool decode);
/*
 * The enano class itself
 */
class Compressor {
   public:
    Compressor();
    Compressor(enano_params *p);
    ~Compressor();

    enano_params *p;
    context_models *cm;
    ali_comp_t *ac;

    uint AVG_CANT, B_CTX, B_MASK, B_CTX_LEN, NS_MODEL_SIZE;
    // Quality
    uint16_t *ctx_avgs_sums;
    uint16_t *ctx_avgs_err_sums;
    uint32_t *ctx_err_avgs_total;

    /* Compression metrics */
    uint64_t base_in, base_out;
    uint64_t qual_in, qual_out;
    uint64_t name_in, name_out;
    uint64_t total_in, total_out;

    char *decode_buf;
    char not_nl[256];
    char not_end_name[256];
    unsigned char QDif[8][14] = {
        {0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7},
        {1, 0, 2, 3, 4, 5, 5, 6, 6, 6, 7},
        {2, 1, 0, 3, 4, 5, 5, 6, 6, 6, 7},
        {3, 2, 1, 0, 3, 4, 5, 5, 6, 6, 7},
        {3, 3, 2, 1, 0, 4, 5, 5, 6, 6, 7},
        {3, 3, 3, 2, 1, 0, 4, 5, 6, 6, 7},
        {3, 3, 3, 2, 2, 1, 0, 4, 5, 6, 7},
        {0, 1, 2, 3, 4, 3, 4, 3, 4, 5, 6}};
    unsigned char QBin[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14,
                               14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                               15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                               15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                               15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
                               15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15};
    /* --- Buffers */
    // Input and output buffers; Output buffer needs to be more than BLK_SIZE
    char out_buf[BLK_SIZE + BLK_SIZE / 2];
    char name_buf[BLK_SIZE / 4];
    char seq_buf[BLK_SIZE];
    char qual_buf[BLK_SIZE];
    int name_len_a[BLK_SIZE / 9];
    int seq_len_a[BLK_SIZE / 9];
    char out_lens[BLK_SIZE];   // seq_len
    char out_names[BLK_SIZE];  // name
    char out_seqs[BLK_SIZE];   // seq
    char out_quals[BLK_SIZE];  // qual
    int sz_lens, sz_names, sz_seqs, sz_quals;
    char *in_buf_lens, *in_buf_names, *in_buf_seqs, *in_buf_quals;

    int out_ind;  // index into out_buf.
    int comp_len;
    int uncomp_len;
    int ns;
    int seq_len;
    int last_len;
    char last_name[1024];  // Last name
    int last_name_len;     // Length of last name
    int last_p_len;        // Length of last common prefix
    int last_s_len;

    ali_list_t *store_reference(file_io_t *in_fq, file_io_t *in_paf, int out_fd, uint32_t &ref_len, uint32_t &ref_sz);
    char *decode_reference(ali_comp_t *ac, int in_fd, uint32_t &ref_len);

    void soft_reset();
    //Copies the statistical stats from cmp to this
    void copy_stats(context_models *ctx_m);

    void encode_len(RangeCoder *rc, int len);
    int decode_len(RangeCoder *rc);

    void encode_name(RangeCoder *rc, char *name, int len);
    int decode_name(RangeCoder *rc, char *name);

    void encode_seq(RangeCoder *rc, char *seq, int len);
    void decode_seq(RangeCoder *rc, char *seq, int len);

    inline uint get_context(unsigned char s, unsigned char q1, unsigned char q2, uint &s_prev_ctx, uint &Q_prev_ctx);
    void encode_qual(RangeCoder *rc, char *seq, char *qual, int len);
    void decode_qual(RangeCoder *rc, char *seq, char *qual, int len);

    int fq_parse_reads(file_io_t *in_fq, file_io_t *in_paf);

    void compress_lens_names_seqs(bool aligned);
    void compress_quals();

    void decompress_name_and_bc_seq();
    void decompress_name_and_bc_seq_aligned();
    void decompress_qual();

    int fq_compress();
    void fq_decompress();

    bool output_block(int out_fd);
};

#endif