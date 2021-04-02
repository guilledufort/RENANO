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
#ifndef PARAMS_H
#define PARAMS_H

#include <math.h>
#include <zlib.h>
#include <omp.h>

#include <unordered_map>

#include "constants.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
    int64_t ali_qtty_cnt = 0, sz_ali_qtty = 0;
    int64_t q_start_cnt = 0, sz_q_start = 0;
    int64_t q_end_cnt = 0, sz_q_end = 0;
    int64_t t_start_cnt = 0, sz_t_start = 0;
    int64_t t_end_cnt = 0, sz_t_end = 0;
    int64_t strand_cnt = 0, sz_strand = 0;
    int64_t t_idx_cnt = 0, sz_t_idx = 0;
    // int64_t ins_mtch_cnt = 0, sz_ins_mtch = 0;
    int64_t match_cnt = 0, sz_match = 0;
    int64_t skip_cnt = 0, sz_skip = 0;
    int64_t ins_cnt = 0, sz_ins = 0;
    int64_t seq_cnt = 0, sz_seq = 0;
    int64_t lens_cnt = 0, sz_lens = 0;
    int64_t tot_alis = 0, tot_matches = 0, tot_strands = 0;
} global_ali_stats_t;

typedef struct {
    /* Compression metrics */
    global_ali_stats_t ali_stats;
    uint64_t ref_sz = 0;
    uint64_t base_in = 0, base_out = 0;
    uint64_t qual_in = 0, qual_out = 0;
    uint64_t name_in = 0, name_out = 0;
    uint64_t total_in = 0, total_out = 0;
    uint64_t total_alignments = 0;
    double enc_time = 0;
} global_stats_t;

typedef struct {
    int8_t aligned = 0;
    // A hashmap to associate reads names with ids
    std::unordered_map<std::string, int> reads_name_hm;
    uint32_t n_read = 0;
    long nxt_idx = 0;
    // A basecall array with a copy of all basecall sequences.
    // BC sequence of read 'idx' is in position 'idx'
    char **bcs_array;
    uint32_t total_sz = 0;
    uint32_t bcs_size = 0;
    char* new_ref;
    uint32_t* ref_stats;
    uint32_t ref_len = 0;
    uint32_t ** snps_stats;
} global_index_t;

static inline void add_g_idx(global_index_t &g_idx, std::string name, char *bc_seq) {
    g_idx.reads_name_hm[name] = g_idx.nxt_idx;
    g_idx.bcs_array[g_idx.nxt_idx++] = bc_seq;
    assert(g_idx.nxt_idx <= g_idx.bcs_size);
}

typedef struct {
    int8_t aligned = 0;
    // A hashmap to associate reads names with ids
    std::unordered_map<std::string, int> reads_name_hm;
    long nxt_idx = 0;

    char **bcs_array;
    uint32_t bcs_size = 0;

} reads_index_t;

/*
 * enano parameter block.
 */
typedef struct
{
    uint8_t klevel;                                  // -k level
    uint8_t llevel;                                  // -l length
    float min_cvg, min_ovlp;
    uint8_t num_threads;                             // Number of threads
    bool max_comp, upd_m, duplicate_names, verbose;  // Max compression mode, and update models in compressors.
    uint8_t blk_upd_freq;                            // Block batch size for update fase in Fast mode
    uint8_t blk_upd_thresh;                          // Number of blocks of update fase.
    uint8_t aligned;                                 // Wether we are doing AVA or REF alignment of the basecall sequences.
    int decompress;                                  // If we are decoding a file.
    int in_fd;                                       // Input file id
    int out_fd;                                      // Output file id
    uint32_t num_reads;
    char *ref_name, *paf_name, *fastq_name;  // Reference file name
    int in_paf;                              // .paf file name
    global_index_t g_idx;                    // Name - id hashmap index.
    uint32_t enc_reads;                      // Number of encoded reads.
    global_stats_t g_stats;
} enano_params;

static inline void initialize_enano_params(enano_params &p) {
    p.klevel = DEFAULT_K_LEVEL;
    p.llevel = DEFAULT_L_LEVEL;
    p.min_cvg = MIN_COVERAGE;
    p.min_ovlp = MIN_ALI_OVLP;
    p.num_threads = DEFAULT_THREADS_NUM;
    p.blk_upd_freq = DEFAULT_BLK_UPD_FREQ;
    p.blk_upd_thresh = DEFAULT_BLK_UPD_THRESH;
    p.max_comp = false;
    p.aligned = 0;
    p.decompress = 0;
    p.in_fd = 0;
    p.in_paf = 0;
    p.ref_name = 0;
    p.out_fd = 1;
    p.num_reads = 0;
    p.enc_reads = 0;
    p.duplicate_names = true;
    p.upd_m = true;
    p.verbose = false;
}

#include <sys/stat.h>

static inline void parse_reference_file(enano_params &p) {
    printf("Parsing reference file... \n");
    double start_time = omp_get_wtime();
    gzFile fp;
    kseq_t *seq;
    int l;

    fp = gzopen(p.ref_name, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        std::string r_name_key = std::string(seq->name.s, seq->name.l);
        char *bc_seq = new char[seq->seq.l];
        memcpy(bc_seq, seq->seq.s, seq->seq.l * sizeof(char));
        add_g_idx(p.g_idx, r_name_key, bc_seq);
        p.g_idx.total_sz += seq->seq.l;
    }
    kseq_destroy(seq);
    gzclose(fp);
    printf("Finished in %.2fs \n", omp_get_wtime() - start_time);
}

static inline void get_num_reads(enano_params &p) {
    //TODO: Optimize this function
    gzFile fp;
    kseq_t *seq;
    int l;

    fp = gzopen(p.fastq_name, "r");
    seq = kseq_init(fp);
    p.num_reads = 0;
    while ((l = kseq_read(seq)) > 0) {
        assert(p.num_reads < INT_FAST32_MAX);
        p.num_reads++;
    }
    kseq_destroy(seq);
    gzclose(fp);
}

static inline void init_g_idx(enano_params &p) {
    if (p.aligned) {
        p.g_idx.aligned = p.aligned;
        if (p.aligned == REF_ALI) {
            p.g_idx.bcs_array = new char *[MAX_CHROM];
            p.g_idx.bcs_size = MAX_CHROM;
            printf("Basecall alignment against reference.\n");
            parse_reference_file(p);
        } else {
            printf("Basecall alignment against stored reference.\n");
            if (!p.decompress) {
                p.g_idx.bcs_array = new char *[MAX_CHROM];
                p.g_idx.bcs_size = MAX_CHROM;
                parse_reference_file(p);
            } else {
                p.g_idx.bcs_array = new char *[1];
                p.g_idx.bcs_size = 1;
            }
            // Get the total number of reads
            get_num_reads(p);
        }
    } else
        printf("Basecall stream not aligned. \n");
}

static inline void free_g_idx(global_index_t &g_idx) {
    if (g_idx.aligned) {
        for (int i = 0; i < g_idx.nxt_idx; i++) {
            delete[] g_idx.bcs_array[i];
        }
        delete[] g_idx.bcs_array;
        std::unordered_map<std::string, int>().swap(g_idx.reads_name_hm);
    }
}

#endif