#include "enc_algs.h"

extern uint32_t tot_a0, equal_a0, tot_0, equal_0;

void print_enc_results(enano_params *p) {
    if (p->aligned) {
        printf("Total %% of aligned base call strings: %0.1f%%, %% of matches: %0.1f%%\n", (double)p->g_stats.ali_stats.tot_alis * 100 / p->g_stats.base_in, (double)p->g_stats.ali_stats.tot_matches * 100 / p->g_stats.ali_stats.tot_alis);
        printf("Stream \t %% of compressed base calls \t Avg. symbol len. (bits)\n");
        double b_ent = (double)p->g_stats.ali_stats.sz_seq * 8 / p->g_stats.ali_stats.seq_cnt;
        printf("B \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_seq * 100 / p->g_stats.base_out,
               b_ent);
        printf("L \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_lens * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_lens * 8 / p->g_stats.ali_stats.lens_cnt);
        printf("Q \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_ali_qtty * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_ali_qtty * 8 / p->g_stats.ali_stats.ali_qtty_cnt);
        printf("S \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_q_start * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_q_start * 8 / p->g_stats.ali_stats.q_start_cnt);
        printf("E \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_q_end * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_q_end * 8 / p->g_stats.ali_stats.q_end_cnt);
        if (p->aligned == REF_ALI)
            printf("U \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_t_idx * 100 / p->g_stats.base_out,
                (double)p->g_stats.ali_stats.sz_t_idx * 8 / p->g_stats.ali_stats.t_idx_cnt);
        printf("S\' \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_t_start * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_t_start * 8 / p->g_stats.ali_stats.t_start_cnt);
        printf("D \t %0.1f%% \t %0.3f, reversed: %0.3f\n", (double)p->g_stats.ali_stats.sz_strand * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_strand * 8 / p->g_stats.ali_stats.strand_cnt, (double)p->g_stats.ali_stats.tot_strands * 100 / p->g_stats.ali_stats.strand_cnt);
        printf("N \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_ins * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_ins * 8 / p->g_stats.ali_stats.ins_cnt);
        printf("K \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_skip * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_skip * 8 / p->g_stats.ali_stats.skip_cnt);
        printf("M \t %0.1f%% \t %0.3f\n", (double)p->g_stats.ali_stats.sz_match * 100 / p->g_stats.base_out,
               (double)p->g_stats.ali_stats.sz_match * 8 / p->g_stats.ali_stats.match_cnt);
    }

    printf(
        "Stream <original size in bytes> -> <compressed size in bytes> "
        "(<compression ratio>)\n");

    printf("IDs   %llu -> %llu (%0.3f)\n", (long long unsigned int)p->g_stats.name_in, (long long unsigned int)p->g_stats.name_out,
           (double)p->g_stats.name_out / p->g_stats.name_in);
    if (p->aligned == STORE_REF_ALI)
        printf("Bases %llu -> %llu (%0.3f) - without ref: (%0.3f) - ref: (%0.3f%%) \n", (long long unsigned int)p->g_stats.base_in, (long long unsigned int)p->g_stats.base_out,
           (double)p->g_stats.base_out / p->g_stats.base_in, (double)(p->g_stats.base_out - p->g_stats.ref_sz) / p->g_stats.base_in,
            (double)p->g_stats.ref_sz*100/p->g_stats.base_out);
    else 
        printf("Bases %llu -> %llu (%0.3f)\n", (long long unsigned int)p->g_stats.base_in, (long long unsigned int)p->g_stats.base_out,
           (double)p->g_stats.base_out / p->g_stats.base_in);
    printf("Quals %llu -> %llu (%0.3f)\n", (long long unsigned int)p->g_stats.qual_in, (long long unsigned int)p->g_stats.qual_out,
           (double)p->g_stats.qual_out / p->g_stats.qual_in);
    printf("Total %llu -> %llu (%0.3f)\n", (long long unsigned int)p->g_stats.total_in, (long long unsigned int)p->g_stats.total_out,
           (double)p->g_stats.total_out / p->g_stats.total_in);
}

void init_global_stats(context_models *cm, global_structs_t *gs) {
    cm->model_seq = new BASE_MODEL<uint8_t>[gs->NS_MODEL_SIZE];
    if (gs->aligned) {
        cm->am = new ali_models;
        cm->am->model_seq = cm->model_seq;
    }
    gs->ctx_avgs_sums = new uint16_t[gs->AVG_CANT];
    gs->ctx_avgs_err_sums = new uint16_t[gs->AVG_CANT];
    gs->ctx_err_avgs_total = new uint32_t[Q_CTX];

    memset(gs->ctx_avgs_sums, 0, gs->AVG_CANT * sizeof(uint16_t));
    memset(gs->ctx_avgs_err_sums, 0, gs->AVG_CANT * sizeof(uint16_t));

    for (uint s_ctx = 0; s_ctx < gs->B_CTX; s_ctx++) {
        for (uint dif = 0; dif < DIF_CANT; dif++) {
            for (uint q_quant = 0; q_quant < Q_LOG_CANT; q_quant++) {
                uint avg_ctx = (s_ctx << TOTAL_Q_LOG) + (dif << Q_LOG) + q_quant;
                gs->ctx_avgs_sums[avg_ctx] = q_quant << AVG_SHIFT;
            }
        }
    }

    memset(gs->ctx_err_avgs_total, 0, Q_CTX * sizeof(uint32_t));
}

void init_global_files_and_variables(enano_params *p, global_structs_t *gs) {
    init_file_io(&gs->in_fq, p->fastq_name, p->in_fd, BLK_SIZE);
    if (p->aligned && !p->decompress)
        init_file_io(&gs->in_paf, p->paf_name, p->in_paf, BLK_SIZE);
    gs->aligned = p->aligned;
    gs->decompress = p->decompress;
    gs->B_CTX_LEN = p->llevel;
    gs->B_CTX = (1 << (gs->B_CTX_LEN * A_LOG));
    gs->AVG_CANT = (gs->B_CTX * Q_CTX);
    gs->NS_MODEL_SIZE = pow5[p->klevel];
}

void init_global_structs(enano_params *p, context_models *cm, global_structs_t *gs) {
    init_global_files_and_variables(p, gs);
    init_global_stats(cm, gs);
}

void delete_global_stats(context_models *cm, global_structs_t *gs) {
    delete[] cm->model_seq;
    if (gs->aligned)
        delete cm->am;
    delete cm;
    delete[] gs->ctx_avgs_sums;
    delete[] gs->ctx_avgs_err_sums;
    delete[] gs->ctx_err_avgs_total;
    close_io(&gs->in_fq);
    if (gs->aligned && !gs->decompress)
        close_io(&gs->in_paf);
}

void copy_average_stats(Compressor *c, global_structs_t *gs) {
    memcpy(c->ctx_avgs_sums, gs->ctx_avgs_sums, gs->AVG_CANT * sizeof(uint16_t));
    memcpy(c->ctx_avgs_err_sums, gs->ctx_avgs_err_sums, gs->AVG_CANT * sizeof(uint16_t));
    memcpy(c->ctx_err_avgs_total, gs->ctx_err_avgs_total, Q_CTX * sizeof(uint32_t));
}

void update_stats(context_models *cm, global_structs_t *gs, Compressor **comps, uc blocks_loaded) {
    uint i;

    for (i = 0; i < gs->AVG_CANT; i++) {
        uint sum_ctx_avgs_sums = 0;
        uint sum_ctx_avgs_err_sums = 0;
        for (uc c = 0; c < blocks_loaded; c++) {
            sum_ctx_avgs_sums += comps[c]->ctx_avgs_sums[i];
            sum_ctx_avgs_err_sums += comps[c]->ctx_avgs_err_sums[i];
        }
        gs->ctx_avgs_sums[i] = round((double)sum_ctx_avgs_sums / blocks_loaded);
        gs->ctx_avgs_err_sums[i] = round((double)sum_ctx_avgs_err_sums / blocks_loaded);
    }

    for (i = 0; i < Q_CTX; i++) {
        uint sum_ctx_err_avgs_total = 0;
        for (uc c = 0; c < blocks_loaded; c++) {
            sum_ctx_err_avgs_total += comps[c]->ctx_err_avgs_total[i];
        }
        gs->ctx_err_avgs_total[i] =
            round((double)sum_ctx_err_avgs_total / blocks_loaded);
    }

    void **models = new void *[blocks_loaded];

    for (i = 0; i < gs->NS_MODEL_SIZE; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_seq[i]);
        }
        cm->model_seq[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < 256; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_name_prefix[i]);
        }
        cm->model_name_prefix[i].mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_name_suffix[i]);
        }
        cm->model_name_suffix[i].mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_name_len[i]);
        }
        cm->model_name_len[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < 8192; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_name_middle[i]);
        }
        cm->model_name_middle[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < CTX_CNT; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->model_qual_quant[i]);
        }
        cm->model_qual_quant[i].mix_array(models, blocks_loaded);
    }

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void *)(&comps[c]->cm->quant_top);
    }
    cm->quant_top.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void *)(&comps[c]->cm->model_len1);
    }
    cm->model_len1.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void *)(&comps[c]->cm->model_len2);
    }
    cm->model_len2.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void *)(&comps[c]->cm->model_len3);
    }
    cm->model_len3.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void *)(&comps[c]->cm->model_same_len);
    }
    cm->model_same_len.mix_array(models, blocks_loaded);

    if (gs->aligned) {
        // Ali quantities
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->ali_qttys_m);
        }
        cm->am->ali_qttys_m.mix_array(models, blocks_loaded);

        // Ali query starts
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_starts_m1);
        }
        cm->am->q_starts_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_starts_m2);
        }
        cm->am->q_starts_m2.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_starts_m3);
        }
        cm->am->q_starts_m3.mix_array(models, blocks_loaded);

        // Ali query ends
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_ends_m1);
        }
        cm->am->q_ends_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_ends_m2);
        }
        cm->am->q_ends_m2.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->q_ends_m3);
        }
        cm->am->q_ends_m3.mix_array(models, blocks_loaded);

        //Ali target starts
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_starts_m1);
        }
        cm->am->t_starts_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_starts_m2);
        }
        cm->am->t_starts_m2.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_starts_m3);
        }
        cm->am->t_starts_m3.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_starts_m4);
        }
        cm->am->t_starts_m4.mix_array(models, blocks_loaded);

        // Ali target ends
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_ends_m1);
        }
        cm->am->t_ends_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_ends_m2);
        }
        cm->am->t_ends_m2.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_ends_m3);
        }
        cm->am->t_ends_m3.mix_array(models, blocks_loaded);

        // Ali strands
        // for (uint c = 0; c < blocks_loaded; c++) {
        //     models[c] = (void *)(&comps[c]->cm->am->strands_m);
        // }
        // cm->am->strands_m.mix_array(models, blocks_loaded);

        // Ali target indexes
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_idxs_m1);
        }
        cm->am->t_idxs_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_idxs_m2);
        }
        cm->am->t_idxs_m2.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_idxs_m3);
        }
        cm->am->t_idxs_m3.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->t_idxs_m4);
        }
        cm->am->t_idxs_m4.mix_array(models, blocks_loaded);

        // cs matches
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_matches_m1);
        }
        cm->am->cs_matches_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_matches_m2);
        }
        cm->am->cs_matches_m2.mix_array(models, blocks_loaded);

        // cs target skips
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_skips_m1);
        }
        cm->am->cs_skips_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_skips_m2);
        }
        cm->am->cs_skips_m2.mix_array(models, blocks_loaded);

        // cs insertions
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_insertions_m1);
        }
        cm->am->cs_insertions_m1.mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void *)(&comps[c]->cm->am->cs_insertions_m2);
        }
        cm->am->cs_insertions_m2.mix_array(models, blocks_loaded);
    }

    delete[] models;
}

/*
 * A blocking read that refuses to return truncated reads.
 */
bool load_data(global_structs_t *gs, Compressor **comps, int update_load,
               uint &blocks_loaded) {
    int sz;

    u_char comp_id = 0;
    blocks_loaded = 0;

    while (comp_id < update_load &&
           (sz = fill_buf_io(&gs->in_fq)) > 0) {
        Compressor *c = comps[comp_id];
        int error = 0;
        if (0 > (error = c->fq_parse_reads(&gs->in_fq, &gs->in_paf))) {
            printf("Failure to parse and/or compress. Error %d \n", error);
            return true;
        }
        c->total_in += sz;
        blocks_loaded++;
        comp_id++;
    }

    return blocks_loaded <= 0;
}

bool load_data_decode(int in_fd, Compressor **comps, int update_load,
                      uint &blocks_loaded) {
    int res = 0;
    unsigned char len_buf[4];
    int sz;

    blocks_loaded = 0;

    u_char comp_id = 0;

    while (comp_id < update_load && (4 == (sz = read(in_fd, len_buf, 4)))) {
        int32_t comp_len = (len_buf[0] << 0) + (len_buf[1] << 8) +
                           (len_buf[2] << 16) + (len_buf[3] << 24);

        int rem_len = comp_len, in_off = 0;

        do {
            errno = 0;
            int tmp_len = read(in_fd, comps[comp_id]->decode_buf + in_off, rem_len);
            if (errno == EINTR && tmp_len == -1)
                continue;

            if (tmp_len == -1) {
                printf("Abort: read failed, %d.\n", errno);
                perror("foo");
                res = -1;
                goto error;
            }
            if (tmp_len == 0) {
                printf("Abort: truncated read, %d.\n", errno);
                res = -1;
                goto error;
            }
            rem_len -= tmp_len;
            in_off += tmp_len;
        } while (rem_len);

        blocks_loaded++;
        comp_id++;
    }
error:
    return (blocks_loaded <= 0) || (res == -1);
}

/*
 * Encode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */
/*
 * Parse one block at a time. Blocks may not terminate on exact fastq
 * boundaries, so we need to know where we ended processing and move
 * that partial fastq entry back to the start for the next block.
 *
 * We write out the block size too so we can decompress block at a time.
 */
int encode(enano_params &p) {
    int res = 0;
    bool decode = false;

    uint block_num = 0;

    uint BLK_UPD_FREQ = p.blk_upd_freq;
    uint BLK_UPD_THRESH = p.blk_upd_thresh + 1;

    u_char cant_compressors = MAX(BLK_UPD_FREQ, p.num_threads);

    Compressor **comps = new Compressor *[cant_compressors];
    for (uint i = 0; i < cant_compressors; i++) {
        comps[i] = new Compressor(&p);
    }
    context_models* cm = new context_models;
    global_structs_t gs;
    init_global_structs(&p, cm, &gs);

    bool finished = false;

    // Update stats
    uint blocks_loaded;
    uint update_blocks = 0;

    uint num_batches = 1 + BLK_UPD_THRESH / BLK_UPD_FREQ;
    uint *update_batches = new uint[1 + BLK_UPD_THRESH / BLK_UPD_FREQ];
    update_batches[0] = 1;
    for (uint i = 1; i < num_batches; i++)
        update_batches[i] = BLK_UPD_FREQ;

    uint batch = 0;
    uint update_load =
        update_batches[batch];  // MIN(BLK_UPD_FREQ, BLK_UPD_THRESH);

    p.enc_reads = 0;
    ali_list_t *g_al = NULL;
    uint32_t ref_len, ref_sz;

    if (p.aligned == STORE_REF_ALI) {
        // We need to store the Reference that is aligned
        Compressor *comp = comps[0];
        g_al = comp->store_reference(&gs.in_paf, p.out_fd, ref_len, ref_sz);
        
        for (uint i = 0; i < cant_compressors; i++)
            comps[i]->ac->al = g_al;
        
        memcpy(cm->model_seq, comp->cm->model_seq, sizeof(BASE_MODEL<uint8_t>) * comp->NS_MODEL_SIZE);
            
        printf("Original ref sz: %u -> Actually used: %u (%0.3f) \n", p.g_idx.total_sz, ref_len, double(ref_len) / p.g_idx.total_sz);
    }

    printf(
        "Starting adaptative encoding with %d threads, for %d (%d + 1) blocks, and update every "
        "%d blocks... \n", p.num_threads,
        BLK_UPD_THRESH, BLK_UPD_THRESH - 1, BLK_UPD_FREQ);

    while (update_blocks < BLK_UPD_THRESH &&
           !(finished = load_data(&gs, comps, update_load, blocks_loaded))) {
#pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            copy_average_stats(comps[i], &gs);
            comps[i]->copy_stats(cm);
            comps[i]->soft_reset();
            comps[i]->fq_compress();
        }

        for (uint i = 0; i < blocks_loaded; i++) {
            if (!comps[i]->output_block(p.out_fd)) {
                printf("Enc1 - Abort: truncated write.\n");
                finished = true;
                res = -1;
                break;
            }
        }
        update_stats(cm, &gs, comps, blocks_loaded);

        update_blocks += blocks_loaded;
        block_num += blocks_loaded;
        batch += 1;
        update_load = MIN((BLK_UPD_THRESH - update_blocks), update_batches[batch]);
    }

    delete[] update_batches;

    printf("Starting parallelized fast encoding...\n");
    // Update context models accumulated probabilities.
    update_AccFreqs(&p, cm, decode);
    // No updates from now on
    p.upd_m = false;
    for (uint i = 0; i < cant_compressors; i++) {
        delete[] comps[i]->cm->model_seq;
        if (p.aligned)
            delete comps[i]->cm->am;
        delete comps[i]->cm;
        comps[i]->cm = cm;
        if (p.aligned) {
            comps[i]->ac->am = cm->am;
            comps[i]->ac->am->model_seq = cm->model_seq;
        }
    }

    while (!finished) {
        uint blocks_loaded = 0;
        finished = load_data(&gs, comps, p.num_threads, blocks_loaded);
#pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            copy_average_stats(comps[i], &gs);
            comps[i]->soft_reset();
            comps[i]->fq_compress();
        }
        for (uint i = 0; i < blocks_loaded; i++) {
            if (!comps[i]->output_block(p.out_fd)) {
                printf("Enc2 - Abort: truncated write.\n");
                finished = true;
                res = -1;
                break;
            }
        }
        block_num += blocks_loaded;
    }

    printf("Total encoded blocks: %d \n", block_num);

    // We initialize total_out in H_SZ for the H_SZ bytes of the header
    p.g_stats.total_in = 0, p.g_stats.total_out = H_SZ;
    for (uint i = 0; i < cant_compressors; i++) {
        p.g_stats.name_in += comps[i]->name_in;
        p.g_stats.name_out += comps[i]->name_out;
        p.g_stats.base_in += comps[i]->base_in;
        p.g_stats.base_out += comps[i]->base_out;
        p.g_stats.qual_in += comps[i]->qual_in;
        p.g_stats.qual_out += comps[i]->qual_out;
        p.g_stats.total_in += comps[i]->total_in;
        p.g_stats.total_out += comps[i]->total_out;
    }

    for (uint i = 0; i < cant_compressors; i++)
        delete comps[i];

    delete[] comps;
    delete_global_stats(cm, &gs);

    if (p.aligned == STORE_REF_ALI) {
        // We add the compressed reference sz to the base out and total out
        p.g_stats.ref_sz = ref_sz;
        p.g_stats.base_out += ref_sz;
        p.g_stats.total_out += ref_sz;
        delete_ali_list_t(g_al);
    }

    print_enc_results(&p);

    return res;
}

/*
 * Decode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */

/*
 * Parse one block at a time. Blocks may not terminate on exact fastq
 * boundaries, so we need to know where we ended processing and move
 * that partial fastq entry back to the start for the next block.
 *
 * We write out the block size too so we can decompress block at a time.
 */
int decode(enano_params &p) {
    int res = 0;
    bool decode = true;

    printf("Starting decoding with %d threads... \n", p.num_threads);

    uint block_num = 0;

    uint BLK_UPD_FREQ = p.blk_upd_freq;
    uint BLK_UPD_THRESH = p.blk_upd_thresh + 1;

    u_char cant_compressors = MAX(BLK_UPD_FREQ, p.num_threads);

    Compressor **comps = new Compressor *[cant_compressors];
    for (uint i = 0; i < cant_compressors; i++) {
        comps[i] = new Compressor(&p);
        comps[i]->decode_buf = new char[BLK_SIZE];
    }
    context_models *cm = new context_models;
    if (p.aligned)
        cm->am = new ali_models;

    global_structs_t gs;
    init_global_structs(&p, cm, &gs);

    bool finished = false;

    // Update stats
    uint blocks_loaded;
    uint update_blocks = 0;
    uint num_batches = 1 + BLK_UPD_THRESH / BLK_UPD_FREQ;
    uint *update_batches = new uint[1 + BLK_UPD_THRESH / BLK_UPD_FREQ];
    update_batches[0] = 1;
    for (uint i = 1; i < num_batches; i++)
        update_batches[i] = BLK_UPD_FREQ;

    uint batch = 0;
    uint update_load = update_batches[batch];

    if (p.aligned == STORE_REF_ALI) {
        // We need to decode the Reference that is aligned
        printf("Decoding stored reference sequence... \n");
        Compressor *comp = comps[0];
        uint32_t ref_len;
        char *ref_seq = comp->decode_reference(comp->ac, p.in_fd, ref_len);
        memcpy(cm->model_seq, comp->cm->model_seq, sizeof(BASE_MODEL<uint8_t>) * comp->NS_MODEL_SIZE);
        add_g_idx(p.g_idx, "", ref_seq);
    }

    p.enc_reads = 0;

    printf("Starting decoding with context model update... \n");

    while (update_blocks < BLK_UPD_THRESH &&
           !(finished =
                 load_data_decode(p.in_fd, comps, update_load, blocks_loaded))) {
#pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            copy_average_stats(comps[i], &gs);
            comps[i]->copy_stats(cm);
            comps[i]->soft_reset();
            comps[i]->fq_decompress();
        }
        // Write output
        for (uint i = 0; i < blocks_loaded; i++) {
            if (comps[i]->uncomp_len !=
                write(p.out_fd, comps[i]->out_buf, comps[i]->uncomp_len)) {
                printf("Dec1 - Abort: truncated write.\n");
                res = -1;
                goto finishdecode;
            }
        }
        update_stats(cm, &gs, comps, blocks_loaded);

        update_blocks += blocks_loaded;
        block_num += blocks_loaded;
        batch += 1;
        update_load = MIN((BLK_UPD_THRESH - update_blocks), update_batches[batch]);
    }

    delete[] update_batches;

    if (!finished) {
        printf("Starting parallelized fast decoding... \n");
        // Update context models accumulated probabilities.
        // This is for fast arithmetic encoding
        update_AccFreqs(&p, cm, decode);
        // No updates from now on
        p.upd_m = false;
        for (uint i = 0; i < cant_compressors; i++) {
            if (p.aligned)
                delete comps[i]->cm->am;
            delete comps[i]->cm;
            comps[i]->cm = cm;
            if (p.aligned) {
                comps[i]->ac->am = cm->am;
                comps[i]->ac->am->model_seq = cm->model_seq;
            }
        }
        // Parallelized decompression with fixed stats
        while (!finished) {
            uint blocks_loaded = 0;
            finished = load_data_decode(p.in_fd, comps, p.num_threads, blocks_loaded);
#pragma omp parallel for
            for (uint i = 0; i < blocks_loaded; i++) {
                copy_average_stats(comps[i], &gs);
                comps[i]->soft_reset();
                comps[i]->fq_decompress();
            }
            for (uint i = 0; i < blocks_loaded; i++) {
                if (comps[i]->uncomp_len !=
                    write(p.out_fd, comps[i]->out_buf, comps[i]->uncomp_len)) {
                    printf("Dec2 - Abort: truncated write.\n");
                    res = -1;
                    goto finishdecode;
                }
            }
            block_num += blocks_loaded;
        }
    }
// We use this goto flag to break the double loop
finishdecode:

    printf("Total decoded blocks: %d\n", block_num);

    for (uint i = 0; i < cant_compressors; i++) {
        delete comps[i]->decode_buf;
        if (comps[i]->cm != cm)
            delete comps[i]->cm;
        delete comps[i];
    }
    delete[] comps;
    delete_global_stats(cm, &gs);
    return res;
}

int encode_st(enano_params &p) {
    int res = 0;
    printf("Starting encoding in Max Compression mode... \n");
    uint block_num = 0;
    uint cant_compressors = 1;

    Compressor **comps = new Compressor *[cant_compressors];
    comps[0] = new Compressor(&p);

    global_structs_t gs;
    init_global_files_and_variables(&p, &gs);

    p.enc_reads = 0;
    ali_list_t *g_al = NULL;
    uint32_t ref_len, ref_sz;

    if (p.aligned == STORE_REF_ALI) {
        // We need to store the Reference that is aligned
        Compressor *comp = comps[0];
        g_al = comp->store_reference(&gs.in_paf, p.out_fd, ref_len, ref_sz);
        comp->ac->al = g_al;
        printf("Original ref sz: %u -> Actually used: %u (%0.3f) \n", p.g_idx.total_sz, ref_len, double(ref_len) / p.g_idx.total_sz);
    }

    // Update stats
    uint blocks_loaded;
    while (!load_data(&gs, comps, 1, blocks_loaded)) {
        // printf("Block: %d.\n", block_num);
        comps[0]->fq_compress();
        if (!comps[0]->output_block(p.out_fd)) {
            printf("Abort: truncated write.\n");
            res = -1;
            break;
        }
        block_num += blocks_loaded;
    }

    printf("Total encoded blocks: %d \n", block_num);

    p.g_stats.name_in = 0, p.g_stats.name_out = 0, p.g_stats.base_in = 0, p.g_stats.base_out = 0, p.g_stats.qual_in = 0,
    p.g_stats.qual_out = 0;
    // We initialize p.g_stats.total_out in 9 for the 9 bytes of the header
    p.g_stats.total_in = 0, p.g_stats.total_out = H_SZ;

    p.g_stats.name_in += comps[0]->name_in;
    p.g_stats.name_out += comps[0]->name_out;
    p.g_stats.base_in += comps[0]->base_in;
    p.g_stats.base_out += comps[0]->base_out;
    p.g_stats.qual_in += comps[0]->qual_in;
    p.g_stats.qual_out += comps[0]->qual_out;
    p.g_stats.total_in += comps[0]->total_in;
    p.g_stats.total_out += comps[0]->total_out;

    if (p.aligned == STORE_REF_ALI) {
        // We add the compressed reference sz to the base out and total out
        p.g_stats.ref_sz = ref_sz;
        p.g_stats.base_out += ref_sz;
        p.g_stats.total_out += ref_sz;
        delete_ali_list_t(g_al);
    }
    delete comps[0];
    delete[] comps;

    print_enc_results(&p);

    return res;
}

int decode_st(enano_params &p) {
    int res = 0;
    printf("Starting decoding Max Compression mode... \n");
    uint block_num = 0;
    u_char cant_compressors = 1;

    Compressor **comps = new Compressor *[cant_compressors];
    comps[0] = new Compressor(&p);
    comps[0]->decode_buf = new char[BLK_SIZE];

    global_structs_t gs;
    init_global_files_and_variables(&p, &gs);
    p.enc_reads = 0;

    if (p.aligned == STORE_REF_ALI) {
        // We need to decode the Reference that is aligned
        printf("Decoding stored reference sequence... \n");
        Compressor *comp = comps[0];
        uint32_t ref_len;
        char *ref_seq = comp->decode_reference(comp->ac, p.in_fd, ref_len);
        add_g_idx(p.g_idx, "", ref_seq);
    }
    // Update stats
    uint blocks_loaded;
    while (!load_data_decode(p.in_fd, comps, 1, blocks_loaded)) {
        comps[0]->fq_decompress();
        // Write output
        if (comps[0]->uncomp_len !=
            write(p.out_fd, comps[0]->out_buf, comps[0]->uncomp_len)) {
            printf("Abort: truncated write.\n");
            res = -1;
            break;
        }
        block_num += blocks_loaded;
    }

    delete comps[0]->decode_buf;
    delete comps[0];
    delete[] comps;
    printf("Total decoded blocks: %d\n", block_num);
    return res;
}
