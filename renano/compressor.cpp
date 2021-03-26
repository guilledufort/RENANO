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

#include "compressor.h"

// #define __DEBUG_LOG__

#define DEBUG_PRINTS

/* -------------------------------------------------------------------------
 * Constructors and destructors
 */

Compressor::Compressor() {
    Compressor(NULL);
}

Compressor::Compressor(enano_params *g_p) {
    p = g_p;
    cm = new context_models;

    B_CTX_LEN = p->llevel;
    B_CTX = (1 << (B_CTX_LEN * A_LOG));
    AVG_CANT = (B_CTX * Q_CTX);
    B_MASK = (B_CTX - 1);
    NS_MODEL_SIZE = pow5[p->klevel];

    for (uint i = 0; i < 256; i++) {
        not_nl[i] = 1;
        not_end_name[i] = 1;
    }
    not_nl['\r'] = not_nl['\n'] = 0;
    not_end_name['\r'] = not_end_name['\n'] = not_end_name[' '] = 0;

    cm->model_seq = new BASE_MODEL<uint8_t>[NS_MODEL_SIZE];
    ctx_avgs_sums = new uint16_t[AVG_CANT];
    ctx_avgs_err_sums = new uint16_t[AVG_CANT];
    ctx_err_avgs_total = new uint32_t[Q_CTX];

    if (p->aligned) {
        cm->am = new ali_models;
        ac = new ali_comp_t;
        ac->am = cm->am;
        ac->am->model_seq = cm->model_seq;
        ac->p = p;
        ac->NS_MODEL_SIZE = pow5[p->klevel];
        init_ali_comp_t(ac, MAX_ALIS_L);
    }

    memset(ctx_avgs_sums, 0, AVG_CANT * sizeof(uint16_t));
    memset(ctx_avgs_err_sums, 0, AVG_CANT * sizeof(uint16_t));

    for (uint s_ctx = 0; s_ctx < B_CTX; s_ctx++) {
        for (uint dif = 0; dif < DIF_CANT; dif++) {
            for (uint q_quant = 0; q_quant < Q_LOG_CANT; q_quant++) {
                uint avg_ctx = (s_ctx << TOTAL_Q_LOG) + (dif << Q_LOG) + q_quant;
                ctx_avgs_sums[avg_ctx] = q_quant << AVG_SHIFT;
            }
        }
    }

    memset(ctx_err_avgs_total, 0, Q_CTX * sizeof(uint32_t));

    /* Name settings */
    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;

    name_in = name_out = 0;
    base_in = base_out = 0;
    qual_in = qual_out = 0;
    total_in = total_out = 0;
}

void update_AccFreqs(enano_params *p, context_models *g_ctx_m, bool decode) {
    uint i;
    for (i = 0; i < 256; i++) {
        g_ctx_m->model_name_prefix[i].updateModelAccFrecs(decode);
        g_ctx_m->model_name_suffix[i].updateModelAccFrecs(decode);
        g_ctx_m->model_name_len[i].updateModelAccFrecs(decode);
    }

    for (i = 0; i < 8192; i++)
        g_ctx_m->model_name_middle[i].updateModelAccFrecs(decode);

    for (i = 0; i < CTX_CNT; i++)
        g_ctx_m->model_qual_quant[i].updateModelAccFrecs(decode);

    g_ctx_m->quant_top.updateModelAccFrecs(decode);
    g_ctx_m->model_len1.updateModelAccFrecs(decode);
    g_ctx_m->model_len2.updateModelAccFrecs(decode);
    g_ctx_m->model_len3.updateModelAccFrecs(decode);
    g_ctx_m->model_same_len.updateModelAccFrecs(decode);

    if (p->aligned) {
        g_ctx_m->am->ali_qttys_m.updateModelAccFrecs(decode);

        g_ctx_m->am->q_starts_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->q_starts_m2.updateModelAccFrecs(decode);
        g_ctx_m->am->q_starts_m3.updateModelAccFrecs(decode);

        g_ctx_m->am->q_ends_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->q_ends_m2.updateModelAccFrecs(decode);
        g_ctx_m->am->q_ends_m3.updateModelAccFrecs(decode);

        g_ctx_m->am->t_starts_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->t_starts_m2.updateModelAccFrecs(decode);
        g_ctx_m->am->t_starts_m3.updateModelAccFrecs(decode);
        g_ctx_m->am->t_starts_m4.updateModelAccFrecs(decode);

        // g_ctx_m->am->t_ends_m1.updateModelAccFrecs(decode);
        // g_ctx_m->am->t_ends_m2.updateModelAccFrecs(decode);
        // g_ctx_m->am->t_ends_m3.updateModelAccFrecs(decode);

        // g_ctx_m->am->strands_m.updateModelAccFrecs(decode);
        
        if (p->aligned == REF_ALI){
            g_ctx_m->am->t_idxs_m1.updateModelAccFrecs(decode);
            g_ctx_m->am->t_idxs_m2.updateModelAccFrecs(decode);
            g_ctx_m->am->t_idxs_m3.updateModelAccFrecs(decode);
            g_ctx_m->am->t_idxs_m4.updateModelAccFrecs(decode);
        }

        g_ctx_m->am->cs_matches_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->cs_matches_m2.updateModelAccFrecs(decode);
        g_ctx_m->am->cs_skips_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->cs_skips_m2.updateModelAccFrecs(decode);
        g_ctx_m->am->cs_insertions_m1.updateModelAccFrecs(decode);
        g_ctx_m->am->cs_insertions_m2.updateModelAccFrecs(decode);
    }
}

void output_int32(int32_t i, int out_fd) {
    char out[4];

    out[0] = (i >> 0) & 0xff;
    out[1] = (i >> 8) & 0xff;
    out[2] = (i >> 16) & 0xff;
    out[3] = (i >> 24) & 0xff;

    int sz = write(out_fd, out, 4);
}

uint32_t read_int32(int in_fd) {
    char in[4];
    read(in_fd, in, 4);
    return DECODE_INT((unsigned char *)(in));
}

ali_list_t *Compressor::store_reference(file_io_t *in_fq, file_io_t *in_paf, int out_fd, uint32_t &ref_len, uint32_t &ref_sz) {
    printf("Storing reference... \n");
    double start_time = omp_get_wtime();
    ali_list_t *al = new ali_list_t;
    init_ali_list_t(al, p->num_reads);
    uint64_t sum_targets = paf_parse_all_alignments(in_paf, al, p, &p->g_idx);

    uint32_t max_len = MIN(p->g_idx.total_sz, sum_targets)/4;
    char* ref_output = new char[max_len];

    ac->rc_seq.output(ref_output);
    ac->rc_seq.StartEncode();
    uint32_t num_seg = 0;
    ref_len = encode_min_reference(al, ac, &p->g_idx, num_seg);
    ac->rc_seq.FinishEncode();

    ref_sz = ac->rc_seq.size_out();
    assert((double)ref_sz*100/max_len < 100);
    printf("ref_sz over max: %0.3f%% (must be less than 100%%)\n", (double)ref_sz*100/max_len);
    printf("Ref entropy: %0.3f \n", (double)ref_sz*8/ref_len);
    printf("ref_len: %ld \n", ref_len);
    printf("ref_sz: %ld \n", ref_sz);

    p->g_idx.ref_len = ref_len;

    output_int32(ref_len, out_fd);
    output_int32(ref_sz, out_fd);
    write(out_fd, ref_output, ref_sz);
    delete []ref_output;

    reset_io(in_paf);
    printf("Finished in %.2fs \n", omp_get_wtime() - start_time);
    return al;
}

char *Compressor::decode_reference(ali_comp_t *ac, int in_fd, uint32_t &ref_len) {
    ref_len = read_int32(in_fd);
    uint32_t ref_sz = read_int32(in_fd);

    printf("ref_len: %ld \n", ref_len);
    printf("ref_sz: %ld \n", ref_sz);

    char *ref_seq = new char[ref_len];
    memset(ref_seq, 0, ref_len * sizeof(char));

    char* ref_buf = new char[ref_sz];
    read(in_fd, ref_buf, ref_sz);

    ac->rc_seq.input(ref_buf);
    ac->rc_seq.StartDecode();
    decode_seq_a(ac, &ac->rc_seq, ref_seq, ref_len);

    ac->rc_seq.FinishDecode();
    delete [] ref_buf;

    return ref_seq;
}

void Compressor::soft_reset() {
    /* Name settings */
    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;

    // reset_ali_comp_t(ac);
}

// Copies global g_ctx_m stats to the compressor.
void Compressor::copy_stats(context_models *g_ctx_m) {
    //Save pointerS cause they are going to get modified by the next lines
    BASE_MODEL<uint8_t> *model_seq_ptr = cm->model_seq;
    ali_models *am_ptr = cm->am;

    *cm = *(g_ctx_m);
    memcpy(model_seq_ptr, g_ctx_m->model_seq, sizeof(BASE_MODEL<uint8_t>) * NS_MODEL_SIZE);

    if (p->aligned) {
        cm->am = am_ptr;
        *(cm->am) = *(g_ctx_m->am);
    }
    //Set pointers again
    cm->model_seq = model_seq_ptr;
    if (p->aligned)
        cm->am->model_seq = model_seq_ptr;
}

Compressor::~Compressor() {
    delete[] ctx_avgs_sums;
    delete[] ctx_avgs_err_sums;
    delete[] ctx_err_avgs_total;
    if (p->aligned) {
        delete_ali_comp_t(ac);
    }
}

/* -------------------------------------------------------------------------
 * Name model
 */
void Compressor::encode_name(RangeCoder *rc, char *name, int len) {
    int p_len, s_len;  // prefix and suffix length
    int i, j, k, last_char;

    // Prefix
    for (i = 0; i < len && i < last_name_len; i++) {
        if (name[i] != last_name[i])
            break;
    }
    p_len = i;

    // Suffix
    for (i = len - 1, j = last_name_len - 1; i >= 0 && j >= 0; i--, j--) {
        if (name[i] != last_name[j])
            break;
    }
    s_len = len - 1 - i;
    if (len - s_len - p_len < 0)
        s_len = len - p_len;

    cm->model_name_prefix[last_p_len].encodeSymbol(rc, p_len, p->upd_m, p->max_comp);
    cm->model_name_suffix[last_s_len].encodeSymbol(rc, s_len, p->upd_m, p->max_comp);
    cm->model_name_len[last_name_len].encodeSymbol(rc, len, p->upd_m, p->max_comp);

    last_p_len = p_len;
    last_s_len = s_len;

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;

        cm->model_name_middle[last_char].encodeSymbol(rc, name[i] & 0x7f, p->upd_m, p->max_comp);

        if (name[i] == ' ' && last_name[j] != ' ')
            j++;
        if (name[i] != ' ' && last_name[j] == ' ')
            j--;
        if (name[i] == ':' && last_name[j] != ':')
            j++;
        if (name[i] != ':' && last_name[j] == ':')
            j--;

        if (name[i] == ':' || name[i] == ' ')
            k = (k + 3) >> 2 << 2;

        lc2 = name[i] == last_name[j];
    }

    memcpy(last_name, name, len);
    last_name_len = len;
}

int Compressor::decode_name(RangeCoder *rc, char *name) {
    int p_len, s_len, len;  // prefix and suffix length
    int i, j, k;
    int last_char;

    p_len = cm->model_name_prefix[last_p_len].decodeSymbol(rc, p->upd_m, p->max_comp);
    s_len = cm->model_name_suffix[last_s_len].decodeSymbol(rc, p->upd_m, p->max_comp);
    len = cm->model_name_len[last_name_len].decodeSymbol(rc, p->upd_m, p->max_comp);

    last_p_len = p_len;
    last_s_len = s_len;

    for (i = 0; i < p_len; i++)
        name[i] = last_name[i];

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        unsigned char c;

        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;
        c = cm->model_name_middle[last_char].decodeSymbol(rc, p->upd_m, p->max_comp);
        //c = 'x';
        name[i] = c;

        if (c == ' ' && last_name[j] != ' ')
            j++;
        if (c != ' ' && last_name[j] == ' ')
            j--;
        if (c == ':' && last_name[j] != ':')
            j++;
        if (c != ':' && last_name[j] == ':')
            j--;

        if (name[i] == ':' || name[i] == ' ')
            k = (k + 3) >> 2 << 2;

        lc2 = c == last_name[j];
    }

    for (j = last_name_len - s_len; i < len; i++, j++)
        name[i] = last_name[j];

    memcpy(last_name, name, len);
    last_name_len = len;

    return len;
}

/* -------------------------------------------------------------------------
 * Sequence length model
 */
void Compressor::encode_len(RangeCoder *rc, int len) {
    if (len != last_len) {
        cm->model_same_len.encodeSymbol(rc, 0, p->upd_m, p->max_comp);
        cm->model_len1.encodeSymbol(rc, len & 0xff, p->upd_m, p->max_comp);
        cm->model_len2.encodeSymbol(rc, (len >> 8) & 0xff, p->upd_m, p->max_comp);
        cm->model_len3.encodeSymbol(rc, (len >> 16) & 0xff, p->upd_m, p->max_comp);
    } else {
        cm->model_same_len.encodeSymbol(rc, 1, p->upd_m, p->max_comp);
    }
}

int Compressor::decode_len(RangeCoder *rc) {
    if (cm->model_same_len.decodeSymbol(rc, p->upd_m, p->max_comp)) {
        return last_len;
    } else {
        int l1 = cm->model_len1.decodeSymbol(rc, p->upd_m, p->max_comp);
        int l2 = cm->model_len2.decodeSymbol(rc, p->upd_m, p->max_comp);
        int l3 = cm->model_len3.decodeSymbol(rc, p->upd_m, p->max_comp);
        last_len = l1 + (l2 << 8) + (l3 << 16);
        return last_len;
    }
}

/* -------------------------------------------------------------------------
 * Sequence model
 */

void Compressor::encode_seq(RangeCoder *rc, char *seq, int len) {
    int last;
    // Corresponds to a sequence of NS consecutive 'N'
    last = 0;

    for (int i = 0; i < len; i++) {
        unsigned char b = L[(unsigned char)seq[i]];
        cm->model_seq[last].encodeSymbol(rc, b, p->upd_m);
        last = UPDATE_CONTEXT(last, b);
    }
}

void Compressor::decode_seq(RangeCoder *rc, char *seq, int len) {
    int last;
    const char *dec = "ACGTN";

    // Corresponds to a sequence of NS consecutive 'N'
    last = 0;

    for (int i = 0; i < len; i++) {
        unsigned char b;
        b = cm->model_seq[last].decodeSymbol(rc, p->upd_m);
        *seq++ = dec[b];
        last = UPDATE_CONTEXT(last, b);
    }
}

/* -------------------------------------------------------------------------
 * Quality model
 */

inline uint Compressor::get_context(unsigned char b, unsigned char q1, unsigned char q2, uint &B_prev_ctx, uint &Q_prev_ctx) {
    int err;
    uint nctx;
    uint avg_err;

    //Update averages
    nctx = (B_prev_ctx << TOTAL_Q_LOG) + Q_prev_ctx;
    err = q1 - (DIV_ROUND(ctx_avgs_sums[nctx], AVG_SHIFT));
    ctx_avgs_sums[nctx] += err;
    uint abs_err = ABS(err);
    ctx_avgs_err_sums[nctx] += abs_err - (DIV_ROUND(ctx_avgs_err_sums[nctx], AVG_SHIFT));
    ctx_err_avgs_total[Q_prev_ctx] += abs_err - (DIV_ROUND(ctx_err_avgs_total[Q_prev_ctx], TOTAL_ERR_SHIFT));

    //We don't consider N for context. That saves memory.
    if (b == 'N')
        b = 'A';
    B_prev_ctx = ((B_prev_ctx << A_LOG) + L[b]) & B_MASK;

    //Current q_ctx
    uint q1_quant = QBin[q1];

    int dif_ctx = 0;
    if (q1 < 7) {
        dif_ctx = QDif[q1][MIN(q2, 10)];
    } else {
        int dif = q2 - q1;

        if (dif == 0) {
            dif_ctx = 0;
        } else {
            if (dif < 0)
                dif_ctx = QDif[7][MIN(-2 * dif - 1, 9)];
            else
                dif_ctx = QDif[7][MIN(2 * dif, 10)];
        }
    }

    Q_prev_ctx = ((dif_ctx << Q_LOG) + q1_quant);

    nctx = (B_prev_ctx << TOTAL_Q_LOG) + Q_prev_ctx;

    uint avg = QBin[DIV_ROUND(ctx_avgs_sums[nctx], AVG_SHIFT)];
    uint total_err_avg = (ctx_err_avgs_total[Q_prev_ctx] >> (TOTAL_ERR_SHIFT - AVG_SHIFT));
    avg_err = ctx_avgs_err_sums[nctx];
    uint err_c = 0;

    if (avg_err < (total_err_avg >> 1))
        err_c = 0;
    else if (avg_err < total_err_avg)
        err_c = 1;
    else if (avg_err < (total_err_avg << 1))
        err_c = 2;
    else
        err_c = 3;

    return (err_c << (TOTAL_Q_LOG + LOG_AVGS)) + (avg << (TOTAL_Q_LOG)) + Q_prev_ctx;
}

void Compressor::encode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    int i, next_b;
    next_b = 1 + B_CTX_LEN / 2;
    uint B_prev_ctx = 0, Q_prev_ctx = 0, ctx = 0;
    uint q1 = 0, q2 = 0;

    // Get first context
    for (i = 0; i < next_b; i++) {
        if (i < len)
            ctx = get_context(seq[i], q1, q2, B_prev_ctx, Q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, B_prev_ctx, Q_prev_ctx);
    }

    for (i = 0; i < len; i++, next_b++) {
        q1 = (qual[i] - '!') & (QMAX - 1);

        if (q1 < QUANT_D_MAX) {
            cm->model_qual_quant[ctx].encodeSymbol(rc, q1, p->upd_m, p->max_comp);
        } else {
            if (p->upd_m) {
                if (p->max_comp) {
                    cm->model_qual_quant[ctx].encodeSymbolOrder(rc, QUANT_D_MAX);
                    cm->quant_top.encodeSymbolOrder(rc, q1 - QUANT_D_MAX);
                } else {   
                    cm->model_qual_quant[ctx].encodeSymbolRegular(rc, QUANT_D_MAX);
                    cm->quant_top.encodeSymbolRegular(rc, q1 - QUANT_D_MAX);
                }
            } else {    
                cm->model_qual_quant[ctx].encodeSymbolNoUpdate(rc, QUANT_D_MAX);
                cm->quant_top.encodeSymbolNoUpdate(rc, q1 - QUANT_D_MAX);
            }
        }

        if (next_b < len)
            ctx = get_context(seq[next_b], q1, q2, B_prev_ctx, Q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, B_prev_ctx, Q_prev_ctx);

        q2 = q1;
    }
}

void Compressor::decode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    int i;
    int q1 = 0, q2 = 0;
    int next_b = 1 + B_CTX_LEN / 2;
    uint b_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;

    for (i = 0; i < next_b; i++) {
        if (i < len)
            ctx = get_context(seq[i], q1, q2, b_prev_ctx, q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, b_prev_ctx, q_prev_ctx);
    }
    assert(ctx < CTX_CNT);
    for (i = 0; i < len; i++, next_b++) {
        unsigned q1;
        q1 = (unsigned char)cm->model_qual_quant[ctx].decodeSymbol(rc, p->upd_m, p->max_comp);
        if (q1 == QUANT_D_MAX) {
            q1 += (unsigned char)cm->quant_top.decodeSymbol(rc, p->upd_m, p->max_comp);
        }

        if (next_b < len)
            ctx = get_context(seq[next_b], q1, q2, b_prev_ctx, q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, b_prev_ctx, q_prev_ctx);

        qual[i] = q1 + '!';

        q2 = q1;
    }
}
/* --------------------------------------------------------------------------
 * Compression functions.
 */
void Compressor::compress_lens_names_seqs(bool aligned) {
    RangeCoder rc_lens, rc_names, rc_seqs;

    rc_lens.output(out_lens);
    rc_names.output(out_names);
    if (!aligned)
        rc_seqs.output(out_seqs);

    rc_lens.StartEncode();
    rc_names.StartEncode();
    if (aligned)
        init_encode(ac);
    else
        rc_seqs.StartEncode();

    char *name_p = name_buf;
    char *seq_p = seq_buf;

    for (int i = 0; i < ns; i++) {
        encode_len(&rc_lens, seq_len_a[i]);
        encode_name(&rc_names, name_p, name_len_a[i]);
        if (aligned) {
            uint j = 0;
            while (not_end_name[(uc)name_p[j]])
                j++;

            std::string r_name = std::string(name_p, j);
            encode_ali_seq(r_name, ac, i, seq_p, seq_len_a[i], &p->g_idx, p->aligned);
        } else
            encode_seq(&rc_seqs, seq_p, seq_len_a[i]);
        name_p += name_len_a[i];
        seq_p += seq_len_a[i];
    }

    rc_lens.FinishEncode();
    rc_names.FinishEncode();
    if (aligned)
        finish_encode(ac);
    else
        rc_seqs.FinishEncode();

    sz_lens = rc_lens.size_out();
    sz_names = rc_names.size_out();
    if (aligned)
        sz_seqs = get_enc_size(ac);
    else
        sz_seqs = rc_seqs.size_out();

    name_in += name_p - name_buf;
    name_out += sz_names;

    base_in += seq_p - seq_buf;
    base_out += sz_seqs;
}

/* Quality values */
void Compressor::compress_quals() {
    char *qual_p = qual_buf;
    char *seq_p = seq_buf;
    RangeCoder rc;

    rc.output(out_quals);
    rc.StartEncode();

    for (int i = 0; i < ns; i++) {
        encode_qual(&rc, seq_p, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
        seq_p += seq_len_a[i];
    }

    rc.FinishEncode();

    sz_quals = rc.size_out();
    qual_in += qual_p - qual_buf;
    qual_out += sz_quals;
}

static inline void write_int32(char *&out_p, int32_t i) {
    *out_p++ = (i >> 0) & 0xff; /* Number of sequences */
    *out_p++ = (i >> 8) & 0xff;
    *out_p++ = (i >> 16) & 0xff;
    *out_p++ = (i >> 24) & 0xff;
}
// #define __DEBUG_BLOCKS__

int Compressor::fq_compress() {
    /* Encode seq len, we have a dependency on this for seq/qual */
    char *out = out_buf + 4;

    compress_lens_names_seqs(p->aligned);

    compress_quals();

    /* Concatenate compressed output into a single block */
    char *out_p = out;

    write_int32(out_p, ns);
    write_int32(out_p, sz_lens);
    write_int32(out_p, sz_names);
    if (p->aligned) {
        write_int32(out_p, ac->sz_ali_qtty);
        write_int32(out_p, ac->sz_seq);
        write_int32(out_p, ac->sz_q_start);
        write_int32(out_p, ac->sz_q_end);
        write_int32(out_p, ac->sz_t_start);
        // write_int32(out_p, ac->sz_t_end);
        write_int32(out_p, ac->sz_strand);
        if (ac->p->aligned == REF_ALI)
            write_int32(out_p, ac->sz_t_idx);
        write_int32(out_p, ac->sz_match);
        write_int32(out_p, ac->sz_skip);
        write_int32(out_p, ac->sz_ins);
    } else
        write_int32(out_p, sz_seqs);
    write_int32(out_p, sz_quals);

#ifdef __DEBUG_BLOCKS__
    printf("ns : %ld\n", ns);
    printf("lens : %ld\n", sz_lens);
    printf("names : %ld\n", sz_names);
    if (p->aligned)
        printf("bases : %ld\n", ac->total_sz);
    else
        printf("bases : %ld\n", sz_seqs);
    printf("quals : %ld\n", sz_quals);
#endif

    memcpy(out_p, out_lens, sz_lens);
    out_p += sz_lens;
    memcpy(out_p, out_names, sz_names);
    out_p += sz_names;
    if (p->aligned) {
        memcpy(out_p, ac->out_ali_qttys, ac->sz_ali_qtty);
        out_p += ac->sz_ali_qtty;
        memcpy(out_p, ac->out_ali_bcs, ac->sz_seq);
        out_p += ac->sz_seq;
        memcpy(out_p, ac->out_q_starts, ac->sz_q_start);
        out_p += ac->sz_q_start;
        memcpy(out_p, ac->out_q_ends, ac->sz_q_end);
        out_p += ac->sz_q_end;
        memcpy(out_p, ac->out_t_starts, ac->sz_t_start);
        out_p += ac->sz_t_start;
        // memcpy(out_p, ac->out_t_ends, ac->sz_t_end);
        // out_p += ac->sz_t_end;
        memcpy(out_p, ac->out_strands, ac->sz_strand);
        out_p += ac->sz_strand;
        if (ac->p->aligned == REF_ALI) {
            memcpy(out_p, ac->out_t_idxs, ac->sz_t_idx);
            out_p += ac->sz_t_idx;
        }
        memcpy(out_p, ac->out_cs_matches, ac->sz_match);
        out_p += ac->sz_match;
        memcpy(out_p, ac->out_cs_skips, ac->sz_skip);
        out_p += ac->sz_skip;
        memcpy(out_p, ac->out_cs_insertions, ac->sz_ins);
        out_p += ac->sz_ins;
    } else {
        memcpy(out_p, out_seqs, sz_seqs);
        out_p += sz_seqs;
    }
    memcpy(out_p, out_quals, sz_quals);
    out_p += sz_quals;

    comp_len = out_p - out;

    assert(comp_len < BLK_SIZE + BLK_SIZE / 2);

    return 0;
}
/* --------------------------------------------------------------------------
 * Decompression functions.
 */

void Compressor::decompress_name_and_bc_seq(void) {
    RangeCoder rc_name, rc_seq;
    rc_seq.input(in_buf_seqs);
    rc_seq.StartDecode();
    rc_name.input(in_buf_names);
    rc_name.StartDecode();

    char *seq_p = seq_buf;
    char *name_p = name_buf;

    int32_t j = 0, k = 0;

    for (int i = 0; i < ns; i++) {
        *name_p++ = '@';
        uint32_t n_len = decode_name(&rc_name, name_p);
        decode_seq(&rc_seq, seq_p, seq_len_a[i]);
        long s_len = seq_len_a[i];

        name_p += n_len;
        *name_p++ = '\n';
        j += n_len + 2;
        seq_p += s_len;
        assert(seq_p <= seq_buf + BLK_SIZE);
    }
    rc_seq.FinishDecode();
    rc_name.FinishDecode();
}

void Compressor::decompress_name_and_bc_seq_aligned(void) {
    RangeCoder rc_name;
    rc_name.input(in_buf_names);
    rc_name.StartDecode();
    init_decode(ac);

    char *seq_p = seq_buf;
    char *name_p = name_buf;

    int32_t j = 0, k = 0;

    // Clear alignment lists
    reset_ali_comp_t(ac);

    for (int i = 0; i < ns; i++) {
        *name_p++ = '@';
        uint32_t n_len = decode_name(&rc_name, name_p);
        decode_ali_seq(seq_p, ac, &p->g_idx, seq_len_a[i]);
        long s_len = seq_len_a[i];

        name_p += n_len;
        *name_p++ = '\n';
        j += n_len + 2;
        seq_p += s_len;
    }
    rc_name.FinishDecode();
    finish_decode(ac);
}

void Compressor::decompress_qual(void) {
    RangeCoder rc_qual;
    rc_qual.input(in_buf_quals);
    rc_qual.StartDecode();

    char *seq_p = seq_buf;
    char *qual_p = qual_buf;
    for (int i = 0; i < ns; i++) {
        decode_qual(&rc_qual, seq_p, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
        seq_p += seq_len_a[i];
    }
    rc_qual.FinishDecode();
}

/* Decompress a single block */
void Compressor::fq_decompress() {
    char *name_p, *seq_p, *qual_p;
    char *in = decode_buf;
    uint32_t sz_seqs;
    uint32_t offset = 0;
    uint32_t nseqs = DECODE_INT((unsigned char *)(in + offset));
    offset += 4;
    uint32_t sz_lens = DECODE_INT((unsigned char *)(in + offset));
    offset += 4;
    uint32_t sz_names = DECODE_INT((unsigned char *)(in + offset));
    offset += 4;
    if (p->aligned) {
        ac->total_sz = 0;
        ac->total_sz += ac->sz_ali_qtty = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_seq = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_q_start = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_q_end = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_t_start = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        // ac->total_sz += ac->sz_t_end = DECODE_INT((unsigned char *)(in + offset));
        // offset += 4;
        ac->total_sz += ac->sz_strand = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        if (ac->p->aligned == REF_ALI) {
            ac->total_sz += ac->sz_t_idx = DECODE_INT((unsigned char *)(in + offset));
            offset += 4;
        }
        ac->total_sz += ac->sz_match = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_skip = DECODE_INT((unsigned char *)(in + offset));
        offset += 4;
        ac->total_sz += ac->sz_ins = DECODE_INT((unsigned char *)(in + offset));
    } else
        sz_seqs = DECODE_INT((unsigned char *)(in + offset));

    offset += 4;
    uint32_t sz_quals = DECODE_INT((unsigned char *)(in + offset));
    offset += 4;
    in += offset;

#ifdef __DEBUG_BLOCKS__
    printf("ns : %ld\n", nseqs);
    printf("lens : %ld\n", sz_lens);
    printf("names : %ld\n", sz_names);
    if (p->aligned)
        printf("bases : %ld\n", ac->total_sz);
    else
        printf("bases : %ld\n", sz_seqs);
    printf("quals : %ld\n", sz_quals);
#endif
    ns = nseqs;

    in_buf_lens = in;
    in += sz_lens;
    in_buf_names = in;
    in += sz_names;
    if (p->aligned) {
        set_inputs(ac, in);
    } else {
        in_buf_seqs = in;
        in += sz_seqs;
    }
    in_buf_quals = in;
    in += sz_quals;

    RangeCoder rc0;
    rc0.input(in_buf_lens);
    rc0.StartDecode();

    for (int i = 0; i < ns; i++)
        seq_len_a[i] = decode_len(&rc0);
    rc0.FinishDecode();

    if (p->aligned)
        decompress_name_and_bc_seq_aligned();
    else
        decompress_name_and_bc_seq();

    decompress_qual();

    /* Stick together the arrays into out_buf */
    out_ind = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    for (int i = 0; i < ns; i++) {
        /* name */
        char *aux_name_p = name_p;
        while ((out_buf[out_ind++] = *name_p++) != '\n')
            ;

        /* seq */
        for (int j = 0; j < seq_len_a[i]; j++)
            out_buf[out_ind++] = *seq_p++;

        out_buf[out_ind++] = '\n';
        out_buf[out_ind++] = '+';
        if (p->duplicate_names) {
            while ((out_buf[out_ind++] = *(++aux_name_p)) != '\n')
                ;
        } else {
            out_buf[out_ind++] = '\n';
        }
        /* qual */
        for (int j = 0; j < seq_len_a[i]; j++) {
            out_buf[out_ind++] = *qual_p++;
        }
        out_buf[out_ind++] = '\n';
    }

    uncomp_len = out_ind;
}

int Compressor::fq_parse_reads(file_io_t *in_fq, file_io_t *in_paf) {
    char *in = in_fq->buf;
    int end_hash = 0;
    int i, j;

    char *name_p = name_buf;
    char *seq_p = seq_buf;
    char *qual_p = qual_buf;

    char *in_end;
    int remainder_length;

    ns = 0;

    /* Parse and separate into name, seq, qual buffers */
    seq_len = 0;

    uint len;

    uint in_len = in_fq->in_len;

    std::string r_name;

    if (p->aligned) 
        reset_ali_comp_t(ac);
    

    for (i = 0; i < in_len;) {
        /* Name */
        if (in[i] != '@')
            return -1;

        j = i;
        i++;

        while (i < in_len && not_end_name[(uc)in[i]])
            i++;

        if (p->aligned)
            r_name = std::string(&in[j + 1], i - j - 1);

        while (i < in_len && not_nl[(uc)in[i]])
            i++;

        len = i - j - 1;
        name_len_a[ns] = len;
        memcpy(name_p, &in[j + 1], len * sizeof(char));
        name_p += len;

        if (++i >= in_len)
            break;

        /* Sequence */
        for (j = i; i < in_len && not_nl[(uc)in[i]]; i++)
            ;

        len = i - j;
        seq_len_a[ns] = len;
        uint seq_pos = j;
        memcpy(seq_p, &in[j], len * sizeof(char));

        seq_p += len;

        if (++i >= in_len)
            break;

        /* +name, assume to be identical to @name */
        if (in[i] != '+')
            return -2;

        for (; i < in_len && not_nl[(uc)in[i]]; i++)
            ;
        if (++i >= in_len)
            break;

        /* Quality */
        if (i + seq_len_a[ns] > in_len) {
            i = in_len + 1;
            break;
        }

        memcpy(qual_p, &in[i], seq_len_a[ns] * sizeof(char));
        qual_p += seq_len_a[ns];
        i += seq_len_a[ns];

        /* In this case i can be equal to in_len*/
        if (++i > in_len)
            break;

        end_hash = i;

        if (seq_len == 0)
            seq_len = seq_len_a[ns];
        else if (seq_len != seq_len_a[ns])
            seq_len = -1;
        
        if (p->aligned == REF_ALI) {
            paf_parse_read_alignments(r_name, in_paf, ac->al, 0, &p->g_idx);
        } else if (p->aligned == STORE_REF_ALI)
            paf_parse_read_cs_strings(r_name, in_paf, ac->al, 0, &p->g_idx);

        ns++;
    }

    /* rl = first_pos - last_pos + 1
     (we dont add 1 because last pos is already one more than needed) */
    remainder_length = i - end_hash;
    in_end = in + end_hash;

    if (i == BLK_SIZE)
        remainder_length += 1;

    /* We maybe ended on a partial fastq entry, so start from there */
    move_remainder_io(in_fq, in_end);//, remainder_length);

    /* Parse the .paf file alignments for the ns number of reads.*/
    if (p->aligned) 
        p->enc_reads += ns;

    return 0;
}

bool Compressor::output_block(int out_fd) {
    if (p->aligned)
        make_ali_stats_global(ac);

    out_buf[0] = (comp_len >> 0) & 0xff;
    out_buf[1] = (comp_len >> 8) & 0xff;
    out_buf[2] = (comp_len >> 16) & 0xff;
    out_buf[3] = (comp_len >> 24) & 0xff;
    comp_len += 4;

    int sz = write(out_fd, out_buf, comp_len);
    total_out += sz;
    return (sz == comp_len);
}
