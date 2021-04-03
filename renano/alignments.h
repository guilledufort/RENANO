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
#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

#include <assert.h>
#include <math.h>
#include <vector>
#include <list>
#include <ctype.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <omp.h>

/* Range Coder:
 * This is using Eugene Shelwien's code from coders6c2.zip.
 */
#include "clr.h"

/*
 * Order 0 models, optimsed for various sizes of alphabet.
 * order0_coder is the original Dmitry Shkarin code, tweaked a bit to
 * make RangeCoder a parameter and to work as a template for adjusting
 * the alphabet size.
 *
 * Simple_model is James Bonfield implementation for FQZComp, adhering to the same API as
 * Dmitry Shkarin's original ORDER_0_CODER model.
 *
 * Base model N is James Bonfield's model specialising in symbols 0, 1, 2, 3, with
 * dedicated code for handling N (it encodes whichever of 0, 1, 2, 3
 * is most common and returns that value).
 */

#include "base_model.h"  // BASE_MODEL
#include "params.h"
#include "io.h"
#include "simple_model.h"  // SIMPLE_MODEL

#define DECODE_INT(a) ((a)[0] + ((a)[1] << 8) + ((a)[2] << 16) + ((a)[3] << 24))
#define UPDATE_CONTEXT(ctx, b) ((ctx * 5 + b) % NS_MODEL_SIZE)

typedef struct
{
    char *l = NULL;
    uint32_t top = 0;
    uint32_t idx = 0;
    uint32_t size = 0;
    uint32_t i_pos = 0;
} char_l_t;

typedef char_l_t cs_string_t;

//Turns len chars starting on *bc_str to its reverse complement
static inline void bc_rev_comp(char *bc_str, int len) {
    std::reverse(bc_str, bc_str + len);
    for (int i = 0; i < len; i++)
        bc_str[i] = REV_BC[int(bc_str[i])];
}

//Turns to uppercase len chars starting on bc_str
static inline void make_upper_case(char *bc_str, int len) {
    for (int i = 0; i < len; i++)
        bc_str[i] = toupper(bc_str[i]);
}

typedef struct {
    uint32_t srt, end, r_srt;
} ali_reg;

//Struct that represents an alignment
typedef struct
{
    bool aligned = 1;
    std::string q_name;
    char *q_seq;
    // uint32_t q_len;
    uint32_t q_start;
    uint32_t q_end;
    uint8_t strand;
    std::string t_name;
    // uint32_t t_len;
    uint32_t t_start, t_start_new;
    uint32_t t_end, t_end_new;
    // uint32_t num_matches;
    uint32_t ali_len;
    // uint8_t map_qual;
    cs_string_t cs;
    uint32_t num_blk;
    int32_t t_idx;
    std::list<ali_reg*>::iterator it_reg;
} alignment_t;

typedef struct
{
    SIMPLE_MODEL<MAX_ALI_R + 1> ali_qttys_m;

    SIMPLE_MODEL<256> q_starts_m1;
    SIMPLE_MODEL<256> q_starts_m2;
    SIMPLE_MODEL<256> q_starts_m3;

    SIMPLE_MODEL<256> q_ends_m1;
    SIMPLE_MODEL<256> q_ends_m2;
    SIMPLE_MODEL<256> q_ends_m3;

    SIMPLE_MODEL<256> t_starts_m1;
    SIMPLE_MODEL<256> t_starts_m2;
    SIMPLE_MODEL<256> t_starts_m3;
    SIMPLE_MODEL<256> t_starts_m4;

    SIMPLE_MODEL<256> t_ends_m1;
    SIMPLE_MODEL<256> t_ends_m2;
    SIMPLE_MODEL<256> t_ends_m3;

    SIMPLE_MODEL<2> strands_m;

    SIMPLE_MODEL<256> t_idxs_m1;
    SIMPLE_MODEL<256> t_idxs_m2;
    SIMPLE_MODEL<256> t_idxs_m3;
    SIMPLE_MODEL<256> t_idxs_m4;

    SIMPLE_MODEL<256> cs_matches_m1;
    SIMPLE_MODEL<256> cs_matches_m2;
    SIMPLE_MODEL<256> cs_skips_m1;
    SIMPLE_MODEL<256> cs_skips_m2;
    SIMPLE_MODEL<256> cs_insertions_m1;
    SIMPLE_MODEL<256> cs_insertions_m2;

    BASE_MODEL<uint8_t> *model_seq;

} ali_models;

typedef struct {
    uint8_t qtty;
    uint32_t alis_pos;
} ali_info;

typedef struct {
    alignment_t *alis;
    uint32_t ali_top = 0;
    uint32_t ali_idx = 0;
    uint32_t max_alis = 0;
    uint32_t num_blk = 0;
    std::list<uint32_t> blk_lims;
    std::list<ali_reg*> ali_regs;
    std::unordered_map<std::string, ali_info> alis_info_hm;
} ali_list_t;

typedef struct
{
    RangeCoder rc_ali_qtty, rc_seq, rc_q_start, rc_q_end, rc_t_start, rc_t_end, rc_strand, rc_t_idx, rc_ins_mtch, rc_match, rc_skip, rc_ins;
    int32_t total_sz = 0, sz_ali_qtty, sz_seq, sz_q_start, sz_q_end, sz_t_start, sz_t_end, sz_strand, sz_t_idx, sz_ins_mtch, sz_match, sz_skip, sz_ins, sz_lens;
    int32_t ali_qtty_cnt = 0, q_start_cnt = 0, q_end_cnt = 0, t_start_cnt = 0, t_end_cnt = 0, strand_cnt = 0, t_idx_cnt = 0, ins_mtch_cnt = 0, match_cnt = 0, skip_cnt = 0, ins_cnt = 0, seq_cnt = 0, lens_cnt = 0;
    int32_t tot_alis = 0, tot_matches = 0, tot_strands = 0;
    /* Ali output buffers */
    char out_ali_qttys[SMALL_BLK_SIZE];
    char out_q_starts[SMALL_BLK_SIZE];
    char out_q_ends[SMALL_BLK_SIZE];
    char out_strands[SMALL_BLK_SIZE];
    char out_t_starts[SMALL_BLK_SIZE];
    char out_t_ends[SMALL_BLK_SIZE];
    char out_t_idxs[SMALL_BLK_SIZE];

    char out_cs_matches[MEDIUM_BLK_SIZE];
    char out_cs_skips[4*SMALL_BLK_SIZE];
    char out_cs_insertions[4*SMALL_BLK_SIZE];
    char out_ali_bcs[4*MEDIUM_BLK_SIZE];

    ali_models *am;
    enano_params *p;
    int prev_ns = 0;

    ali_list_t *al;

    int *L;
    uint32_t NS_MODEL_SIZE;

    uint8_t last_ins_match = INSER;
    uint8_t min_match_len = MIN_MATCH_LEN;
    uint last_bc_ctx = 0;

} ali_comp_t;

void init_ali_list_t(ali_list_t *al, uint32_t buf_size);
void init_ali_comp_t(ali_comp_t *ac, uint32_t buf_size);
void init_encode(ali_comp_t *ac);
void finish_encode(ali_comp_t *ac);
void set_inputs(ali_comp_t *ac, char *&in_buf);
void init_decode(ali_comp_t *ac);
void finish_decode(ali_comp_t *ac);
void reset_ali_comp_t(ali_comp_t *ac);
void delete_ali_list_t(ali_list_t *al);
void delete_ali_comp_t(ali_comp_t *ac);
void make_ali_stats_global(ali_comp_t *ac);
uint32_t search_index(std::string r_name, global_index_t *g_idx);

uint32_t encode_min_reference(ali_list_t *al, ali_comp_t *ac, global_index_t *g_idx, uint32_t &num_seg);

void encode_ali_seq(std::string q_name, ali_comp_t *ac, char *seq, uint32_t q_len, global_index_t *g_idx, int aligned);
int decode_ali_seq(char *seq_p, ali_comp_t *ac, global_index_t *g_idx, uint32_t seq_len);
void decode_seq_a(ali_comp_t *ac, RangeCoder *rc, char *seq, uint32_t len);
uint32_t get_enc_size(ali_comp_t *ac);
uint32_t decode_seg_n(ali_comp_t *ac, RangeCoder *rc);
#endif