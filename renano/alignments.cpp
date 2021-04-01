#include "alignments.h"

void encode_qtty(ali_comp_t *ac, RangeCoder *rc, int8_t qtty) {
    ac->ali_qtty_cnt++;
    ac->am->ali_qttys_m.encodeSymbol(rc, qtty, ac->p->upd_m, ac->p->max_comp);
}
int8_t decode_qtty(ali_comp_t *ac, RangeCoder *rc) {
    return ac->am->ali_qttys_m.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
}
void encode_q_start_end(ali_comp_t *ac, RangeCoder *rc_q_start, RangeCoder *rc_q_end, uint32_t q_start, uint32_t q_end, uint32_t prev_q_start) {
    ac->q_start_cnt++;
    ac->q_end_cnt++;
    uint32_t ali_len = q_end - q_start;
    q_start -= prev_q_start;

    ac->am->q_starts_m1.encodeSymbol(rc_q_start, q_start & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->q_starts_m2.encodeSymbol(rc_q_start, (q_start >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->q_starts_m3.encodeSymbol(rc_q_start, (q_start >> 16) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert(((q_start >> 24) & 0xff) == 0);

    ac->am->q_ends_m1.encodeSymbol(rc_q_end, ali_len & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->q_ends_m2.encodeSymbol(rc_q_end, (ali_len >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->q_ends_m3.encodeSymbol(rc_q_end, (ali_len >> 16) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert(((ali_len >> 24) & 0xff) == 0);
}
void decode_q_start_end(ali_comp_t *ac, RangeCoder *rc_q_start, RangeCoder *rc_q_end, uint32_t &q_start, uint32_t &q_end, uint32_t prev_q_start) {
    uint32_t l1 = ac->am->q_starts_m1.decodeSymbol(rc_q_start, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->q_starts_m2.decodeSymbol(rc_q_start, ac->p->upd_m, ac->p->max_comp);
    uint32_t l3 = ac->am->q_starts_m3.decodeSymbol(rc_q_start, ac->p->upd_m, ac->p->max_comp);
    q_start = l1 + (l2 << 8) + (l3 << 16);
    q_start += prev_q_start;

    l1 = ac->am->q_ends_m1.decodeSymbol(rc_q_end, ac->p->upd_m, ac->p->max_comp);
    l2 = ac->am->q_ends_m2.decodeSymbol(rc_q_end, ac->p->upd_m, ac->p->max_comp);
    l3 = ac->am->q_ends_m3.decodeSymbol(rc_q_end, ac->p->upd_m, ac->p->max_comp);
    q_end = l1 + (l2 << 8) + (l3 << 16) + q_start;
}
void encode_t_start(ali_comp_t *ac, RangeCoder *rc_t_start, uint32_t t_start) {
    ac->t_start_cnt++;
    ac->am->t_starts_m1.encodeSymbol(rc_t_start, t_start & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_starts_m2.encodeSymbol(rc_t_start, (t_start >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_starts_m3.encodeSymbol(rc_t_start, (t_start >> 16) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_starts_m4.encodeSymbol(rc_t_start, (t_start >> 24) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert((((uint64_t)t_start >> 32) & 0xff) == 0);
}
int32_t decode_t_start(ali_comp_t *ac, RangeCoder *rc_t_start) {
    uint32_t l1 = ac->am->t_starts_m1.decodeSymbol(rc_t_start, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->t_starts_m2.decodeSymbol(rc_t_start, ac->p->upd_m, ac->p->max_comp);
    uint32_t l3 = ac->am->t_starts_m3.decodeSymbol(rc_t_start, ac->p->upd_m, ac->p->max_comp);
    uint32_t l4 = ac->am->t_starts_m4.decodeSymbol(rc_t_start, ac->p->upd_m, ac->p->max_comp);
    return l1 + (l2 << 8) + (l3 << 16) + (l4 << 24);
}
void encode_strand(ali_comp_t *ac, RangeCoder *rc, int8_t strand) {
    ac->strand_cnt++;
    ac->am->strands_m.encodeSymbolBinary(rc, strand);
}
int8_t decode_strand(ali_comp_t *ac, RangeCoder *rc) {
    return ac->am->strands_m.decodeSymbolBinary(rc);
}
void encode_t_idx(ali_comp_t *ac, RangeCoder *rc, uint32_t t_idx) {
    ac->t_idx_cnt++;
    ac->am->t_idxs_m1.encodeSymbol(rc, t_idx & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_idxs_m2.encodeSymbol(rc, (t_idx >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_idxs_m3.encodeSymbol(rc, (t_idx >> 16) & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->t_idxs_m4.encodeSymbol(rc, (t_idx >> 24) & 0xff, ac->p->upd_m, ac->p->max_comp);
}
uint32_t decode_t_idx(ali_comp_t *ac, RangeCoder *rc) {
    uint32_t t_idx;
    uint32_t l1 = ac->am->t_idxs_m1.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->t_idxs_m2.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l3 = ac->am->t_idxs_m3.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l4 = ac->am->t_idxs_m4.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    t_idx = l1 + (l2 << 8) + (l3 << 16) + (l4 << 24);
    return t_idx;
}
void encode_seg_n(ali_comp_t *ac, RangeCoder *rc, uint32_t seg_n) {
    encode_t_start(ac, rc, seg_n);
}
uint32_t decode_seg_n(ali_comp_t *ac, RangeCoder *rc) {
    return decode_t_start(ac, rc);
}
void encode_cs_match(ali_comp_t *ac, RangeCoder *rc, uint32_t match) {
    ac->match_cnt++;
    ac->tot_matches += match;
    ac->am->cs_matches_m1.encodeSymbol(rc, match & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->cs_matches_m2.encodeSymbol(rc, (match >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert(((match >> 16) & 0xff) == 0);
}
uint32_t decode_cs_match(ali_comp_t *ac, RangeCoder *rc) {
    uint32_t l1 = ac->am->cs_matches_m1.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->cs_matches_m2.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    return l1 + (l2 << 8);
}
void encode_cs_skip(ali_comp_t *ac, RangeCoder *rc, uint32_t skip) {
    ac->skip_cnt++;
    ac->am->cs_skips_m1.encodeSymbol(rc, skip & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->cs_skips_m2.encodeSymbol(rc, (skip >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert(((skip >> 16) & 0xff) == 0);
}
uint32_t decode_cs_skip(ali_comp_t *ac, RangeCoder *rc) {
    uint32_t l1 = ac->am->cs_skips_m1.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->cs_skips_m2.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    return l1 + (l2 << 8);
}
void encode_cs_ins(ali_comp_t *ac, RangeCoder *rc, uint32_t ins) {
    ac->ins_cnt++;
    ac->am->cs_insertions_m1.encodeSymbol(rc, ins & 0xff, ac->p->upd_m, ac->p->max_comp);
    ac->am->cs_insertions_m2.encodeSymbol(rc, (ins >> 8) & 0xff, ac->p->upd_m, ac->p->max_comp);
    assert(((ins >> 16) & 0xff) == 0);
}
uint32_t decode_cs_ins(ali_comp_t *ac, RangeCoder *rc) {
    uint32_t l1 = ac->am->cs_insertions_m1.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    uint32_t l2 = ac->am->cs_insertions_m2.decodeSymbol(rc, ac->p->upd_m, ac->p->max_comp);
    return l1 + (l2 << 8);
}
void encode_seq(ali_comp_t *ac, RangeCoder *rc, char *seq, uint32_t len) {
    uint NS_MODEL_SIZE = ac->NS_MODEL_SIZE;
    // Corresponds to a sequence of NS consecutive 'N'
    for (uint32_t i = 0; i < len; i++) {
        ac->seq_cnt++;
        unsigned char b = L[(unsigned char)seq[i]];
        assert(ac->last_bc_ctx < NS_MODEL_SIZE);
        assert(&ac->am->model_seq[ac->last_bc_ctx] != NULL);
        ac->am->model_seq[ac->last_bc_ctx].encodeSymbol(rc, b, ac->p->upd_m);
        ac->last_bc_ctx = UPDATE_CONTEXT(ac->last_bc_ctx, b);
    }
}
void encode_not_b(ali_comp_t *ac, RangeCoder *rc, char b_c, char not_b_c) {
    uint NS_MODEL_SIZE = ac->NS_MODEL_SIZE;
    // Corresponds to a sequence of NS consecutive 'N'
    ac->seq_cnt++;
    unsigned char b = L[(unsigned char)b_c];
    unsigned char not_b = L[(unsigned char)not_b_c];
    assert(ac->last_bc_ctx < NS_MODEL_SIZE);
    assert(&ac->am->model_seq[ac->last_bc_ctx] != NULL);
    ac->am->model_seq[ac->last_bc_ctx].encodeSymbolNotB(rc, b, not_b, ac->p->upd_m);
    ac->last_bc_ctx = UPDATE_CONTEXT(ac->last_bc_ctx, b);
}
void decode_seq_a(ali_comp_t *ac, RangeCoder *rc, char *seq, uint32_t len) {
    uint NS_MODEL_SIZE = ac->NS_MODEL_SIZE;
    const char *dec = "ACGTN";

    for (uint32_t i = 0; i < len; i++) {
        unsigned char b;
        b = ac->am->model_seq[ac->last_bc_ctx].decodeSymbol(rc, ac->p->upd_m);
        *seq++ = dec[b];
        ac->last_bc_ctx = UPDATE_CONTEXT(ac->last_bc_ctx, b);
    }
}
void decode_not_b(ali_comp_t *ac, RangeCoder *rc, char *seq, char not_b_c) {
    uint NS_MODEL_SIZE = ac->NS_MODEL_SIZE;
    const char *dec = "ACGTN";
    unsigned char not_b = L[(unsigned char)not_b_c];
    unsigned char b = ac->am->model_seq[ac->last_bc_ctx].decodeSymbolNotB(rc, not_b, ac->p->upd_m);
    *seq++ = dec[b];
    ac->last_bc_ctx = UPDATE_CONTEXT(ac->last_bc_ctx, b);
    
}
void upd_seq_model(ali_comp_t *ac, char *seq, int len) {
    uint NS_MODEL_SIZE = ac->NS_MODEL_SIZE;
    for (int i = 0; i < len; i++) {
        unsigned char b = L[(unsigned char)seq[i]];
        ac->last_bc_ctx = UPDATE_CONTEXT(ac->last_bc_ctx, b);
    }
}
void init_ali_list_t(ali_list_t *al, uint32_t max_reads) {
    al->alis = new alignment_t[max_reads * AVG_ALI_PER_R];
    al->max_alis = max_reads * AVG_ALI_PER_R;
    al->ali_top = 0;
    al->ali_idx = 0;
    al->alis_info_hm.reserve(max_reads);
}
void init_ali_comp_t(ali_comp_t *ac, uint32_t max_reads) {
    if (ac->p->aligned == REF_ALI || ac->p->decompress) {
        ac->al = new ali_list_t;
        init_ali_list_t(ac->al, max_reads);
    }
    ac->last_ins_match = INSER;
    ac->last_bc_ctx = 0;
}
void init_encode(ali_comp_t *ac) {
    if (ac->p->aligned == REF_ALI) {
        ac->al->ali_idx = 0;
    }

    ac->rc_ali_qtty.output(ac->out_ali_qttys);
    ac->rc_seq.output(ac->out_ali_bcs);
    ac->rc_q_start.output(ac->out_q_starts);
    ac->rc_q_end.output(ac->out_q_ends);
    ac->rc_t_start.output(ac->out_t_starts);
    ac->rc_strand.output(ac->out_strands);
    if (ac->p->aligned == REF_ALI)
        ac->rc_t_idx.output(ac->out_t_idxs);
    ac->rc_match.output(ac->out_cs_matches);
    ac->rc_skip.output(ac->out_cs_skips);
    ac->rc_ins.output(ac->out_cs_insertions);

    ac->rc_ali_qtty.StartEncode();
    ac->rc_seq.StartEncode();
    ac->rc_q_start.StartEncode();
    ac->rc_q_end.StartEncode();
    ac->rc_t_start.StartEncode();
    ac->rc_strand.StartEncodeBinary();
    if (ac->p->aligned == REF_ALI)
        ac->rc_t_idx.StartEncode();
    ac->rc_match.StartEncode();
    ac->rc_skip.StartEncode();
    ac->rc_ins.StartEncode();
}
void finish_encode(ali_comp_t *ac) {
    ac->rc_ali_qtty.FinishEncode();
    ac->rc_seq.FinishEncode();
    ac->rc_q_start.FinishEncode();
    ac->rc_q_end.FinishEncode();
    ac->rc_t_start.FinishEncode();
    ac->rc_strand.FinishEncodeBinary();
    if (ac->p->aligned == REF_ALI)
        ac->rc_t_idx.FinishEncode();
    ac->rc_match.FinishEncode();
    ac->rc_skip.FinishEncode();
    ac->rc_ins.FinishEncode();
}
uint32_t get_enc_size(ali_comp_t *ac) {
    ac->total_sz = 0;
    ac->total_sz += ac->sz_ali_qtty = ac->rc_ali_qtty.size_out();
    ac->total_sz += ac->sz_seq = ac->rc_seq.size_out();
    ac->total_sz += ac->sz_q_start = ac->rc_q_start.size_out();
    ac->total_sz += ac->sz_q_end = ac->rc_q_end.size_out();
    ac->total_sz += ac->sz_t_start = ac->rc_t_start.size_out();
    ac->total_sz += ac->sz_strand = ac->rc_strand.size_out();
    if (ac->p->aligned == REF_ALI)
        ac->total_sz += ac->sz_t_idx = ac->rc_t_idx.size_out();
    ac->total_sz += ac->sz_match = ac->rc_match.size_out();
    ac->total_sz += ac->sz_skip = ac->rc_skip.size_out();
    ac->total_sz += ac->sz_ins = ac->rc_ins.size_out();

#ifdef __PRINT_ALI_STATS__
    printf("Total sz: %d \n", ac->total_sz);
    printf("sz_ali_qtty: %0.3f\n", (double)ac->sz_ali_qtty / ac->total_sz);
    printf("sz_seq: %0.3f\n", (double)ac->sz_seq / ac->total_sz);
    printf("sz_q_start: %0.3f\n", (double)ac->sz_q_start / ac->total_sz);
    printf("sz_q_end: %0.3f\n", (double)ac->sz_q_end / ac->total_sz);
    printf("sz_t_start: %0.3f\n", (double)ac->sz_t_start / ac->total_sz);
    printf("sz_strand: %0.3f\n", (double)ac->sz_strand / ac->total_sz);
    if (ac->p->aligned == REF_ALI)
        printf("sz_t_idx: %0.3f\n", (double)ac->sz_t_idx / ac->total_sz);
    printf("sz_match: %0.3f\n", (double)ac->sz_match / ac->total_sz);
    printf("sz_skip: %0.3f\n", (double)ac->sz_skip / ac->total_sz);
    printf("sz_ins: %0.3f\n", (double)ac->sz_ins / ac->total_sz);
#endif

    // printf("sz_ali_qtty: %0.3f\n", (double)ac->sz_ali_qtty * 100 / SMALL_BLK_SIZE);
    assert(ac->sz_ali_qtty <= SMALL_BLK_SIZE);
    // printf("sz_seq: %0.3f\n", (double)ac->sz_seq * 100 / (4*MEDIUM_BLK_SIZE));
    assert(ac->sz_seq <= 4*MEDIUM_BLK_SIZE);
    // printf("sz_q_start: %0.3f\n", (double)ac->sz_q_start * 100 / SMALL_BLK_SIZE);
    assert(ac->sz_q_start <= SMALL_BLK_SIZE);
    // printf("sz_q_end: %0.3f\n", (double)ac->sz_q_end * 100 / SMALL_BLK_SIZE);
    assert(ac->sz_q_end <= SMALL_BLK_SIZE);
    // printf("sz_t_start: %0.3f\n", (double)ac->sz_t_start * 100 / SMALL_BLK_SIZE);
    assert(ac->sz_t_start <= SMALL_BLK_SIZE);
    // printf("sz_t_end: %0.3f\n", (double)ac->sz_t_end * 100 / SMALL_BLK_SIZE);
    // assert(ac->sz_t_end <= SMALL_BLK_SIZE);
    // printf("sz_strand: %0.3f\n", (double)ac->sz_strand * 100 / SMALL_BLK_SIZE);
    assert(ac->sz_strand <= SMALL_BLK_SIZE);
    // printf("sz_t_idx: %0.3f\n", (double)ac->sz_t_idx * 100 / SMALL_BLK_SIZE);
    if (ac->p->aligned == REF_ALI)
        assert(ac->sz_t_idx <= SMALL_BLK_SIZE);
    // printf("sz_match: %0.3f\n", (double)ac->sz_match * 100 / (MEDIUM_BLK_SIZE));
    assert(ac->sz_match <= MEDIUM_BLK_SIZE);
    // printf("sz_skip: %0.3f\n", (double)ac->sz_skip * 100 / (4*SMALL_BLK_SIZE));
    assert(ac->sz_skip <= 4*SMALL_BLK_SIZE);
    // printf("sz_ins: %0.3f\n", (double)ac->sz_ins * 100 / (4*SMALL_BLK_SIZE));
    assert(ac->sz_ins <= 4*SMALL_BLK_SIZE);
    return ac->total_sz;
}

void make_ali_stats_global(ali_comp_t *ac) {
    ac->p->g_stats.ali_stats.sz_ali_qtty += ac->sz_ali_qtty;
    ac->p->g_stats.ali_stats.sz_seq += ac->sz_seq;
    ac->p->g_stats.ali_stats.sz_lens += ac->sz_lens;
    ac->p->g_stats.ali_stats.sz_q_start += ac->sz_q_start;
    ac->p->g_stats.ali_stats.sz_q_end += ac->sz_q_end;
    ac->p->g_stats.ali_stats.sz_t_start += ac->sz_t_start;
    ac->p->g_stats.ali_stats.sz_t_end += ac->sz_t_end;
    ac->p->g_stats.ali_stats.sz_strand += ac->sz_strand;
    if (ac->p->aligned == REF_ALI)
        ac->p->g_stats.ali_stats.sz_t_idx += ac->sz_t_idx;
    ac->p->g_stats.ali_stats.sz_match += ac->sz_match;
    ac->p->g_stats.ali_stats.sz_skip += ac->sz_skip;
    ac->p->g_stats.ali_stats.sz_ins += ac->sz_ins;
    ac->p->g_stats.ali_stats.tot_alis += ac->tot_alis;
    ac->p->g_stats.ali_stats.tot_matches += ac->tot_matches;
    ac->p->g_stats.ali_stats.tot_strands += ac->tot_strands;
    ac->p->g_stats.ali_stats.ali_qtty_cnt += ac->ali_qtty_cnt;
    ac->p->g_stats.ali_stats.seq_cnt += ac->seq_cnt;
    ac->p->g_stats.ali_stats.lens_cnt += ac->lens_cnt;
    ac->p->g_stats.ali_stats.q_start_cnt += ac->q_start_cnt;
    ac->p->g_stats.ali_stats.q_end_cnt += ac->q_end_cnt;
    ac->p->g_stats.ali_stats.t_start_cnt += ac->t_start_cnt;
    ac->p->g_stats.ali_stats.t_end_cnt += ac->t_end_cnt;
    ac->p->g_stats.ali_stats.strand_cnt += ac->strand_cnt;
    if (ac->p->aligned == REF_ALI)
        ac->p->g_stats.ali_stats.t_idx_cnt += ac->t_idx_cnt;
    ac->p->g_stats.ali_stats.match_cnt += ac->match_cnt;
    ac->p->g_stats.ali_stats.skip_cnt += ac->skip_cnt;
    ac->p->g_stats.ali_stats.ins_cnt += ac->ins_cnt;
}

void set_inputs(ali_comp_t *ac, char *&in_buf) {
    ac->rc_ali_qtty.input(in_buf);
    in_buf += ac->sz_ali_qtty;
    ac->rc_seq.input(in_buf);
    in_buf += ac->sz_seq;
    ac->rc_q_start.input(in_buf);
    in_buf += ac->sz_q_start;
    ac->rc_q_end.input(in_buf);
    in_buf += ac->sz_q_end;
    ac->rc_t_start.input(in_buf);
    in_buf += ac->sz_t_start;
    ac->rc_t_end.input(in_buf);
    in_buf += ac->sz_t_end;
    ac->rc_strand.input(in_buf);
    in_buf += ac->sz_strand;
    if (ac->p->aligned == REF_ALI) {
        ac->rc_t_idx.input(in_buf);
        in_buf += ac->sz_t_idx;
    }
    ac->rc_match.input(in_buf);
    in_buf += ac->sz_match;
    ac->rc_skip.input(in_buf);
    in_buf += ac->sz_skip;
    ac->rc_ins.input(in_buf);
    in_buf += ac->sz_ins;
}

void init_decode(ali_comp_t *ac) {
    ac->rc_ali_qtty.StartDecode();
    ac->rc_seq.StartDecode();
    ac->rc_q_start.StartDecode();
    ac->rc_q_end.StartDecode();
    ac->rc_t_start.StartDecode();
    ac->rc_t_end.StartDecode();
    ac->rc_strand.StartDecodeBinary();
    if (ac->p->aligned == REF_ALI)
        ac->rc_t_idx.StartDecode();
    ac->rc_match.StartDecode();
    ac->rc_skip.StartDecode();
    ac->rc_ins.StartDecode();
}
void finish_decode(ali_comp_t *ac) {
    ac->rc_ali_qtty.FinishEncode();
    ac->rc_seq.FinishDecode();
    ac->rc_q_start.FinishDecode();
    ac->rc_q_end.FinishDecode();
    ac->rc_t_start.FinishDecode();
    ac->rc_t_end.FinishDecode();
    ac->rc_strand.FinishDecodeBinary();
    if (ac->p->aligned == REF_ALI)
        ac->rc_t_idx.FinishDecode();
    ac->rc_match.FinishDecode();
    ac->rc_skip.FinishDecode();
    ac->rc_ins.FinishDecode();
}

void reset_ali_comp_t(ali_comp_t *ac) {
    ac->ali_qtty_cnt = 0;
    ac->q_start_cnt = 0;
    ac->q_end_cnt = 0;
    ac->t_start_cnt = 0;
    ac->t_end_cnt = 0;
    ac->strand_cnt = 0;
    ac->t_idx_cnt = 0;
    ac->match_cnt = 0;
    ac->skip_cnt = 0;
    ac->ins_cnt = 0;
    ac->seq_cnt = 0;
    ac->tot_alis = 0;
    ac->tot_matches = 0;
    ac->tot_strands = 0;

    if (ac->p->aligned == REF_ALI) {
        ac->al->ali_top = 0;
        ac->al->ali_idx = 0;
        ac->al->alis_info_hm.clear();
    }

    ac->last_ins_match = INSER;
    ac->last_bc_ctx = 0;
}
void delete_ali_list_t(ali_list_t *al) {
    delete[] al->alis;
    al->ali_regs.clear();
    std::unordered_map<std::string, ali_info>().swap(al->alis_info_hm);
    delete al;
}
void delete_ali_comp_t(ali_comp_t *ac) {
    if (ac->p->aligned == REF_ALI)
        delete_ali_list_t(ac->al);
    delete ac;
}
uint32_t search_index(std::string r_name, global_index_t *g_idx) {
    return g_idx->reads_name_hm[r_name];
}
static inline char *search_index_for_tseq(int32_t idx, global_index_t *g_idx) {
    return g_idx->bcs_array[idx];
}

static inline void encode_triplet(ali_comp_t *ac, char* t_seq, uint32_t i_t, char* q_seq, uint32_t i_q,
                                uint32_t skip_len, uint32_t ins_len, uint32_t m_len, bool &first) {
    encode_cs_skip(ac, &ac->rc_skip, skip_len);
    encode_cs_ins(ac, &ac->rc_ins, ins_len);
    encode_cs_match(ac, &ac->rc_match, m_len);
    //CHANGE
    if (ac->p->aligned == REF_ALI && (skip_len > 0 && ins_len > 0)) {
        encode_not_b(ac, &ac->rc_seq, q_seq[i_q], t_seq[i_t]);
        if (ins_len > 1)
            encode_seq(ac, &ac->rc_seq, &q_seq[i_q + 1], ins_len - 1);
    } else {
        first = false;
        encode_seq(ac, &ac->rc_seq, &q_seq[i_q], ins_len);
    }
    upd_seq_model(ac, &q_seq[i_q + ins_len], m_len);
}

static inline void encode_parts(ali_comp_t *ac, char* t_seq, uint32_t i_t, char* q_seq, uint32_t i_q,
                                uint32_t skip_len, uint32_t ins_len, uint32_t m_len, bool &first) {

    i_q -= ins_len;
    while (skip_len > MAX_INS_LEN || ins_len > MAX_SKIP_LEN ) {
         if (ins_len > MAX_INS_LEN) {
            // printf("MAX_INS: %u", ins_len);
            encode_triplet(ac, t_seq, i_t, q_seq, i_q, 0, MAX_INS_LEN, 0, first);
            i_q += MAX_INS_LEN;
            ins_len -= MAX_INS_LEN;
        } else if (skip_len > MAX_SKIP_LEN) {
            // printf("MAX_SKIP: %u", skip_len);
            bool first_aux = first;
            encode_triplet(ac, t_seq, i_t, q_seq, i_q, MAX_SKIP_LEN, 0, 0, first_aux);
            i_t += skip_len;
            skip_len -= MAX_SKIP_LEN;
        }
    }
    
    encode_triplet(ac, t_seq, i_t, q_seq, i_q, skip_len, ins_len, MIN(m_len, MAX_MATCH_LEN), first);

    if (m_len > MAX_MATCH_LEN) {
        // printf("MAX_MATCH: %u", m_len);
        i_q += ins_len + MAX_MATCH_LEN;
        m_len -= MAX_MATCH_LEN;
        while (m_len > MAX_MATCH_LEN) {
            encode_triplet(ac, t_seq, i_t, q_seq, i_q, 0, 0, MAX_MATCH_LEN, first);
            i_q += MAX_MATCH_LEN;
            m_len -= MAX_MATCH_LEN;
        }
        encode_triplet(ac, t_seq, i_t, q_seq, i_q, 0, 0, m_len, first);
    }
}

static inline void decode_parts(ali_comp_t *ac, uint32_t i_q, uint32_t i_t, uint32_t &skip_len,
                                 uint32_t &ins_len, uint32_t &match_len, char* seq_buf, char* t_seq, bool &first) {
    skip_len = decode_cs_skip(ac, &ac->rc_skip);
    ins_len = decode_cs_ins(ac, &ac->rc_ins);
    match_len = decode_cs_match(ac, &ac->rc_match);
    //CHANGE
    if (ac->p->aligned == REF_ALI && (skip_len > 0 && ins_len > 0)) {
        decode_not_b(ac, &ac->rc_seq, &seq_buf[i_q], t_seq[i_t]);
        if (ins_len > 1)
            decode_seq_a(ac, &ac->rc_seq, &seq_buf[i_q + 1], ins_len - 1);
    } else {
        first = false;
        decode_seq_a(ac, &ac->rc_seq, &seq_buf[i_q], ins_len);
    }
    memcpy(&seq_buf[i_q + ins_len], &t_seq[i_t + skip_len], match_len);
    upd_seq_model(ac, &seq_buf[i_q + ins_len], match_len);
}
/* It decodes the alignment ali with the values present in ali_list to position seq_buf*/
static inline void decode_ali(alignment_t *ali, ali_comp_t *ac, global_index_t *g_idx, char *seq_buf) {
    uint32_t ali_len = ali->q_end - ali->q_start;
    uint32_t match_len, skip_len, ins_len;
    char *t_seq = search_index_for_tseq(ali->t_idx, g_idx);
    uint32_t i_t = ali->t_start, i_q = 0;
    bool first = true;
    while (i_q < ali_len) {
            decode_parts(ac, i_q, i_t, skip_len,
                        ins_len, match_len, seq_buf, t_seq, first);
            i_q += ins_len + match_len;
            i_t += skip_len + match_len;
    }
    assert(i_q == ali_len);

    make_upper_case(seq_buf, ali_len);
    if (ali->strand == REV) {
        bc_rev_comp(seq_buf, ali_len);
    }
}

static inline uint get_int_cs(char *cs, uint cs_len, uint &i) {
    uint res = 0;
    while (i < cs_len && std::isdigit(cs[i]))
        res = res * 10 + (int)(cs[i++] - '0');
    return res;
}
static inline uint get_slen_cs(char *cs, uint cs_len, uint &i) {
    uint res = 0;
    while (i < cs_len && IS_BC[(int)cs[i]]) {
        res++;
        i++;
    }
    return res;
}
static inline void enc_cs_string(alignment_t *ali, ali_comp_t *ac, global_index_t *g_idx, char *seq) {
    uint32_t ali_len = ali->q_end - ali->q_start;
    uint32_t ins_len = 0, i = 0, op_len = 0, skip_len = 0;
    char *cs = ali->cs.l;
    uint32_t cs_len = ali->cs.top;
    char *t_seq = search_index_for_tseq(ali->t_idx, g_idx);
    uint32_t i_t = ali->t_start, i_q = 0;

    char *q_seq;
    if (ali->strand == FOR) {
        q_seq = &seq[ali->q_start];
    } else {
        q_seq = new char[ali_len];
        memcpy(q_seq, &seq[ali->q_start], ali_len);
        bc_rev_comp(q_seq, ali_len);
    }

    bool first = true;

    while (i_q < ali_len) {
        char op = cs[i++];
        //MATCH MODE
        if (op == ':') {
            op_len = get_int_cs(cs, cs_len, i);

            if (i_q + op_len > ali_len)
                op_len = ali_len - i_q;

            if (op_len < ac->min_match_len) {
                skip_len += op_len;
                ins_len += op_len;
                i_q += op_len;
            } else {
                encode_parts(ac, t_seq, i_t, q_seq, i_q,
                                skip_len, ins_len, op_len, first);
                i_q += op_len;
                i_t += skip_len + op_len;
                skip_len = 0;
                ins_len = 0;
            }
            // INSERTION
        } else if (op == '+') {
            op_len = get_slen_cs(cs, cs_len, i);
            if (i_q + op_len > ali_len)
                op_len = ali_len - i_q;
            ins_len += op_len;
            i_q += op_len;
            // SNP
        } else if (op == '*') {
            i += 2;
            i_q += 1;
            ins_len += 1;
            skip_len += 1;
            // DELETION
        } else if (op == '-') {
            skip_len += get_slen_cs(cs, cs_len, i);
        }
    }
    assert(i_q == ali_len);
    // There can be a remaining insertion, as cs_strings ar segmented so they do not overlap
    if (ins_len > 0)
        encode_parts(ac, t_seq, i_t, q_seq, i_q,
                                skip_len, ins_len, 0, first);
    
    if (ali->strand == REV)
        delete[] q_seq;
}
static inline void enc_cs_string_sr(alignment_t *ali, ali_comp_t *ac, global_index_t *g_idx, char *seq) {
    uint32_t ali_len = ali->q_end - ali->q_start;
    uint32_t ins_len = 0, i = 0, op_len = 0, skip_len = 0, s_gap = 0, del_len = 0;
    char *cs = ali->cs.l;
    uint32_t cs_len = ali->cs.top;
    char *t_seq_old = search_index_for_tseq(ali->t_idx, g_idx);
    // i_t_new is the start alignment position on the new artificial contiguous reference
    uint32_t i_t_new = ali->t_start_new, i_q = 0;
    // i_t_old is the original start of the alignment. We use this position to recover the bases associated to each
    // insertion opertaion.
    uint32_t i_t_old = ali->t_start;

    char *q_seq;
    if (ali->strand == FOR) {
        q_seq = &seq[ali->q_start];
    } else {
        q_seq = new char[ali_len];
        memcpy(q_seq, &seq[ali->q_start], ali_len);
        bc_rev_comp(q_seq, ali_len);
    }
    bool first = true;
    
    // it holds a pointer to a list of valid regions.
    // it is always initialized to the first region relevant to the alignemnt.
    // The valid region can either overlap with the alignemnt, or be the next valid region.
    std::list<ali_reg*>::iterator it = ali->it_reg;
    // s_gap denotes the part of the skip that resides in a discarded section of the reference
    s_gap = 0;
    if ((*it)->srt > i_t_new)
        s_gap = (*it)->srt - i_t_new;

    while (i_q < ali_len) {
        char op = cs[i++];
        //MATCH MODE
        if (op == ':') {
            op_len = get_int_cs(cs, cs_len, i);

            if (i_q + op_len > ali_len)
                op_len = ali_len - i_q;

            uint z = 0;
            while (z < op_len) {
                // In case the current position i_t_new + skip_len has surpassed the last valid region,
                // we get the next one.
                while (i_t_new + skip_len >= (*it)->end) {
                    uint aux = (*it)->end;
                    it++;
                    s_gap += (*it)->srt - aux;
                }
                // s_len represents the next section we have to skip to reach the first position available in the new reference.
                // s_len can be 0 if we are inside a valid region of the reference, or it can be larger than op_len,
                // in which case we are skipping the match operator entirely.
                // We dont have to update the real skip len, as it is not present in te encoded reference.
                uint s_len = MAX(0, MIN((int)((*it)->srt - (i_t_new + skip_len)), (int)(op_len - z)));
                if (s_len > 0) {
                    z += s_len;
                    ins_len += s_len;
                    i_q += s_len;
                    skip_len += s_len;
                }
                // If there is at least a remaining section of the match inside of a valid region.
                if (z < op_len) {
                    // We calculate how many bases of the match are inside the valid region.
                    uint m_len = MIN((int)(op_len - z), (int)((*it)->end - (i_t_new + skip_len)));
                    z += m_len;
                    // If the match is too small we add it to the insertion and skip operations
                    if (m_len < ac->min_match_len) {
                        skip_len += m_len;
                        ins_len += m_len;
                        i_q += m_len;
                    } else {
                        // If it is the first triplet encoding we need to update the initial position of the alignment
                        encode_parts(ac, t_seq_old, i_t_old + s_gap, q_seq, i_q,
                                skip_len - s_gap, ins_len, m_len, first);
                        
                        i_q += m_len;
                        i_t_new += skip_len + m_len;
                        i_t_old += skip_len + m_len;
                        s_gap = 0;
                        skip_len = 0;
                        ins_len = 0;
                    }
                }
            }
            // INSERTION
        } else if (op == '+') {
            op_len = get_slen_cs(cs, cs_len, i);
            if (i_q + op_len > ali_len)
                op_len = ali_len - i_q;
            ins_len += op_len;
            i_q += op_len;
            // SNP
        } else if (op == '*') {
            i += 2;
            i_q += 1;
            ins_len += 1;
            skip_len += 1;
            // DELETION
        } else if (op == '-') {
            del_len = get_slen_cs(cs, cs_len, i);
            skip_len += del_len;
        }
    }

    assert(i_q == ali_len);
    // There can be a remaining insertion, as cs_strings are segmented so they do not overlap
     if (ins_len > 0)
         encode_parts(ac, t_seq_old, i_t_old, q_seq, i_q,
                                0, ins_len, 0, first);
    
    if (ali->strand == REV)
        delete[] q_seq;
}
static inline void encode_ali(alignment_t *ali, ali_comp_t *ac, global_index_t *g_idx, char *seq, uint32_t prev_q_start) {
    ac->tot_alis += ali->q_end - ali->q_start;
    ac->tot_strands += ali->strand;
    encode_q_start_end(ac, &ac->rc_q_start, &ac->rc_q_end, ali->q_start, ali->q_end, prev_q_start);
    encode_strand(ac, &ac->rc_strand, ali->strand);
    if (ac->p->aligned == STORE_REF_ALI) {
        // We encode the position of the first valid region as the alignment start.
        uint32_t reg_offset = 0;
        if (ali->t_start_new > (*ali->it_reg)->srt)
            reg_offset = ali->t_start_new - (*ali->it_reg)->srt;
        
        encode_t_start(ac, &ac->rc_t_start, (*ali->it_reg)->r_srt + reg_offset);
        // We encode the encoding transformation.
        enc_cs_string_sr(ali, ac, g_idx, seq);
    } else {
        encode_t_idx(ac, &ac->rc_t_idx, ali->t_idx);
        encode_t_start(ac, &ac->rc_t_start, ali->t_start);
        enc_cs_string(ali, ac, g_idx, seq);
    }
}

void encode_ali_seq(std::string q_name, ali_comp_t *ac, char *seq, uint32_t q_len, global_index_t *g_idx, int aligned) {
    alignment_t *ali_ptr;
    uint8_t qtty = 0;
    ali_info ali_i;

    std::unordered_map<std::string, ali_info>::iterator got = ac->al->alis_info_hm.find(q_name);
    if (got != ac->al->alis_info_hm.end()) {
        ali_i = got->second;
        qtty = ali_i.qtty;
        ali_ptr = &ac->al->alis[ali_i.alis_pos];
    }

    uint32_t next_enc_pos = 0;
    uint32_t prev_q_start = 0;
    uint8_t real_qtty = 0;
    for (uint32_t j = 0; j < qtty; j++) {
        alignment_t ali = *ali_ptr++;

        if (ali.aligned == true) {
            real_qtty++;
            if (ali.q_start > next_enc_pos)
                encode_seq(ac, &ac->rc_seq, seq + next_enc_pos, ali.q_start - next_enc_pos);

            encode_ali(&ali, ac, g_idx, seq, prev_q_start);
            prev_q_start = ali.q_start;
            next_enc_pos = ali.q_end;
        }
        if (ali.cs.l != NULL) {
            delete[] ali.cs.l;
            ali.cs.l = NULL;
        }
    }
    encode_qtty(ac, &ac->rc_ali_qtty, real_qtty);
    if (aligned != STORE_REF_ALI)
        ac->al->ali_idx += qtty;
    // Encode remaining basecall sequence
    if (next_enc_pos < q_len)
        encode_seq(ac, &ac->rc_seq, seq + next_enc_pos, q_len - next_enc_pos);
}

void get_dec_ali(alignment_t &ali, ali_comp_t *ac, uint32_t prev_q_start) {
    decode_q_start_end(ac, &ac->rc_q_start, &ac->rc_q_end, ali.q_start, ali.q_end, prev_q_start);
    ali.t_start = decode_t_start(ac, &ac->rc_t_start);
    ali.ali_len = ali.q_end - ali.q_start;
    if (ac->p->aligned == STORE_REF_ALI)
        ali.t_idx = 0;
    else
        ali.t_idx = decode_t_idx(ac, &ac->rc_t_idx);
    ali.strand = decode_strand(ac, &ac->rc_strand);
}

#define PRINT_ALIS

int decode_ali_seq(char *seq_p, ali_comp_t *ac, global_index_t *g_idx, uint32_t seq_len) {
    ac->al->ali_idx = 0;
    alignment_t ali;

    uint8_t qtty = decode_qtty(ac, &ac->rc_ali_qtty);
    assert(qtty <= MAX_ALI_R);

    uint32_t next_enc_pos = 0;
    uint32_t prev_q_start = 0;
    for (uint32_t j = 0; j < qtty; j++) {
        get_dec_ali(ali, ac, prev_q_start);
        if (ali.q_start > next_enc_pos) {
            decode_seq_a(ac, &ac->rc_seq, seq_p + next_enc_pos, ali.q_start - next_enc_pos);
        }

        decode_ali(&ali, ac, g_idx, seq_p + ali.q_start);
        next_enc_pos = ali.q_end;
        prev_q_start = ali.q_start;
    }
    // Encode remaining basecall sequence
    if (next_enc_pos < seq_len) {
        decode_seq_a(ac, &ac->rc_seq, seq_p + next_enc_pos, seq_len - next_enc_pos);
    }
    return 0;
}


typedef struct {
    alignment_t* ali;
    uint32_t pos;
    bool end;
} ali_mark;

// If end1 is larger than end2, then in case of a tie the ends are processed before the beginnings
bool ali_marks_sorter(ali_mark am1, ali_mark am2) {
    return (am1.ali->t_idx < am2.ali->t_idx) ||
           ((am1.ali->t_idx  == am2.ali->t_idx) && (am1.pos < am2.pos))|| 
           ((am1.ali->t_idx  == am2.ali->t_idx) && (am1.pos == am2.pos) && (am1.end > am2.end));
}

char* fake_ref;

uint32_t encode_min_reference(ali_list_t *al, ali_comp_t *ac, global_index_t *g_idx, uint32_t &num_seg) {
    //Create a marks array that stores starting and ending positions of every alignment
    ali_mark * ali_marks = new ali_mark[al->ali_top*2];
    for (uint i = 0; i < al->ali_top; i++) {
        ali_marks[2*i].ali = &al->alis[i];
        ali_marks[2*i].pos = al->alis[i].t_start;
        ali_marks[2*i].end = false;
        ali_marks[2*i+1].ali = &al->alis[i];
        ali_marks[2*i+1].pos = al->alis[i].t_end;
        ali_marks[2*i+1].end = true;
    }
    //Sort the marks by position
    std::sort(ali_marks, ali_marks + 2 * al->ali_top, ali_marks_sorter);
    uint32_t reg_st = 0, reg_end = 0, r_start = 0, ovlp = 0;

    uint i = 0;
    uint32_t last_enc = 0, next_r_pos = 0, ref_len = 0;
    ali_reg* ar = NULL;
    uint32_t num_reg = 0, t_idx_ovlp = 0;
    
    for (i = 0; i < al->ali_top*2; i++) {
        alignment_t* ali = ali_marks[i].ali;
        if (ali->aligned) {
            //BEGIN MARK
            if (ali_marks[i].end == false) {
                ovlp++;
                num_reg++;
                if (ar == NULL) {
                    ar = new ali_reg;
                    al->ali_regs.push_back(ar);
                }
                //if its the start of a new region
                if (num_reg == 1) {
                    r_start = ali->t_start_new;
                }
                //if its the start of a new overlapping
                if (ovlp == 2) {
                    reg_st = ali->t_start_new;
                    t_idx_ovlp = ali->t_idx;
                }
                ali->t_start_new += next_r_pos - r_start;
                ali->it_reg = --al->ali_regs.end();;
            //END MARK
            } else {
                ovlp--;
                //if its the end of an overlapping region
                if (ovlp == 1) {
                    reg_end = ali->t_end_new;
                    uint32_t reg_len = reg_end - reg_st;
                    ar->r_srt = ref_len;
                    ar->srt = reg_st + next_r_pos - r_start;
                    ar->end = reg_end + next_r_pos - r_start;
                    assert(ar->srt >= last_enc);

                    encode_seq(ac, &ac->rc_seq, &g_idx->bcs_array[t_idx_ovlp][reg_st], reg_len);
                    ref_len += reg_len;

                    last_enc = ar->end;
                    ar = NULL;
                    num_seg++;
                }
                //if its the end of a region
                if (ovlp == 0) {
                    //if there are no overlappings in the region
                    if (num_reg == 1) {
                        ali->aligned = false;
                        al->ali_regs.pop_back();
                        delete ar;
                        ar = NULL;
                    } else {
                        uint32_t r_len = ali->t_end_new - r_start;
                        ali->t_end_new += next_r_pos - r_start;
                        next_r_pos += r_len;
                    }
                    num_reg = 0;
                } else {
                    ali->t_end_new += next_r_pos - r_start;
                }    
            }
        }
    }
    //We add one dummy region at the end to simplify border cases
    ar = new ali_reg;
    ar->srt = next_r_pos;
    ar->end = ar->srt = next_r_pos+1;
    al->ali_regs.push_back(ar);
    printf("Total segments: %d \n", num_seg);

    delete[] ali_marks;
    assert(next_r_pos <= g_idx->total_sz);
    return ref_len;
}