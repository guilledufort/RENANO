#include "paf_process.h"

static inline uint32_t get_int(char *&buf, char delim) {
    uint32_t res = 0;
    while (*buf != delim)
        res = res * 10 + (int)(*buf++ - '0');
    buf++;
    return res;
}

static inline std::string get_string(char *&buf, char delim) {
    int j = 0;
    char *buf_aux = buf;
    while (*buf++ != delim)
        j++;
    return std::string(buf_aux, j);
}

static inline void move_buf_to_nl(char *&buf) {
    while (*buf++ != '\n')
        ;
}

static inline uint32_t get_cs_i_pos(char *&buf, file_io_t *in_paf) {
    while (!(*buf == 'c' && *(buf + 1) == 's'))
        buf++;
    return buf - in_paf->buf;
}

//Pre: buf is assumed to be at the begining of the cs string
static inline void get_cs_string(cs_string_t &res, char *&buf) {
    buf += 5;
    char *buf_aux = buf;
    int i = 0;
    while (*buf_aux != '\n' && *buf_aux != '\r') {
        buf_aux++;
        i += 1;
    }
    res.top = i;
    res.l = new char[i];
    memcpy(res.l, buf, res.top * sizeof(char));
    buf = buf_aux;
}

bool ali_sorter_dec(alignment_t *ac, alignment_t *ar) {
    return ac->q_end > ar->q_end;
}
bool ali_sorter_inc(alignment_t *ac, alignment_t *ar) {
    return ac->q_start < ar->q_start;
}

void add_alis(std::string q_name, ali_list_t *al, alignment_t *ali_group[MAX_ALI_R], uint32_t j) {
    //encode valid alis
    uint8_t ali_qtty = 0;
    ali_info ali_i;
    ali_i.alis_pos = al->ali_top;
    for (int i = 0; i < j; i++) {
        if (ali_group[i]->aligned) {
            assert(al->ali_top < al->max_alis);
            al->alis[al->ali_top++] = *ali_group[i];
            ali_qtty++;
        } 
        else {
            //Free the cs string
            if (ali_group[i]->cs.l != NULL) {
                delete[] ali_group[i]->cs.l;
                ali_group[i]->cs.l = NULL;
            }
        }
    }
    if (ali_qtty > 0) {
        ali_i.qtty = ali_qtty;
        al->alis_info_hm[q_name] = ali_i;
    }
}

/* Retruns the final number of valid alignments*/
/* This function receives an alignment and decides wether to add it or not to the ali_lists. */
/* In case there is more than one alignment for the same read, they are preprocessed so there are no overlappings. *. */
static inline int32_t preprocess_alignments(alignment_t *ali_group[MAX_ALI_R], uint8_t j, ali_list_t *al, global_index_t *g_idx) {
    /* First we sort the alis in decending order by q_end position */
    std::sort(ali_group, ali_group + j, ali_sorter_dec);
    int32_t q_s = ali_group[0]->q_end;
    int32_t sum_targets = 0;
    for (int i = 0; i < j; i++) {
        if ((q_s - (int)ali_group[i]->q_start < MIN_ALI_LEN) ||            // If the alignment is too small we dont want to encode it
            (ali_group[i]->t_idx == -1) ||                                 // If the t_name is not in the read index
            (ali_group[i]->q_end > q_s && ali_group[i]->strand == REV)) {  // TODO: It is a possible bug that reverse alignments get cut when changing the q_end
            ali_group[i]->aligned = false;
        } else {
            ali_group[i]->q_end = MIN(q_s, ali_group[i]->q_end);
            q_s = ali_group[i]->q_start - 1;
            sum_targets += ali_group[i]->t_end - ali_group[i]->t_start;
        }
    }
    std::sort(ali_group, ali_group + j, ali_sorter_inc);
    return sum_targets;
}

static inline void read_alignment(alignment_t *ali, char *&buf, global_index_t *g_idx, bool read_cs, file_io_t *in_paf) {
    ali->q_name = get_string(buf, '\t');
    get_int(buf, '\t'); //q_len
    ali->q_start = get_int(buf, '\t');
    ali->q_end = get_int(buf, '\t');

    if (*buf == '*') {
        ali->aligned = false;
        move_buf_to_nl(buf);
        return;
    }
    ali->aligned = true;
    ali->strand = (*buf++ == '+') ? FOR : REV;
    ali->t_name = get_string(++buf, '\t');
    ali->t_idx = search_index(ali->t_name, g_idx);
    get_int(buf, '\t'); //t_len
    ali->t_start_new = ali->t_start = get_int(buf, '\t');
    ali->t_end_new = ali->t_end = get_int(buf, '\t');
    get_int(buf, '\t'); //num_matches
    get_int(buf, '\t'); //ali_len
    get_int(buf, '\t'); //map_qual
    ali->cs.i_pos = get_cs_i_pos(buf, in_paf);
    assert(in_paf->buf[ali->cs.i_pos] == 'c');
    if (read_cs)
        get_cs_string(ali->cs, buf);
    move_buf_to_nl(buf);
}

bool has_next_ali(char *buf_ptr, file_io_t *in_paf) {
    while ((buf_ptr < in_paf->buf + in_paf->in_len) && *buf_ptr != '\n')
        buf_ptr++;
    return buf_ptr < in_paf->buf + in_paf->in_len;
}

void paf_parse_read_cs_strings(std::string q_name, file_io_t *in_paf, ali_list_t *al, int ns, global_index_t *g_idx) {
    if (in_paf->in_len == 0)
        fill_buf_io(in_paf);
    std::unordered_map<std::string, ali_info>::iterator got = al->alis_info_hm.find(q_name);
    if (got != al->alis_info_hm.end()) {
        ali_info ali_i = got->second;
        uint8_t qtty = ali_i.qtty;
        alignment_t *ali_ptr = &al->alis[ali_i.alis_pos];
        //check for change in num_blk
        uint32_t cant = 0;
        for (uint32_t j = 0; j < qtty; j++) {
            if (ali_ptr[j].num_blk == al->num_blk) {
                cant++;
                if (ali_ptr[j].aligned == true) {
                    char* aux_buf = &in_paf->buf[ali_ptr[j].cs.i_pos];
                    get_cs_string(ali_ptr[j].cs, aux_buf);
                }
            }
        }
        //in case some of the alignments start in the next block
        if (cant < qtty) {
            al->num_blk++;
            uint32_t blk_lim = al->blk_lims.front();
            al->blk_lims.pop_front();
            move_remainder_io(in_paf, in_paf->buf + blk_lim);
            fill_buf_io(in_paf);
            for (uint32_t j = 0; j < qtty; j++) {
                if (ali_ptr[j].num_blk == al->num_blk) {
                    cant++;
                    if (ali_ptr[j].aligned == true) {
                        char* aux_buf = &in_paf->buf[ali_ptr[j].cs.i_pos];
                        get_cs_string(ali_ptr[j].cs,aux_buf);
                    }
                }
            }
        }
        assert(cant == qtty);
    }
}

void paf_parse_read_alignments(std::string q_name, file_io_t *in_paf, ali_list_t *al, int ns, global_index_t *g_idx) {
    alignment_t ali;
    uint32_t j;
    int sz;
    alignment_t *ali_group[MAX_ALI_R];
    bool read_cs = true;
    for (j = 0; j < MAX_ALI_R; j++)
        ali_group[j] = new alignment_t;
    j = 0;
    while ((in_paf->crt_ptr < in_paf->buf + in_paf->in_len) || (sz = fill_buf_io(in_paf, in_paf->crt_ptr)) > 0) {
        char *buf_ptr = in_paf->crt_ptr;
        if (has_next_ali(buf_ptr, in_paf)) {
            char *aux_buf = buf_ptr;
            std::string nxt_name = get_string(aux_buf, '\t');
            if (nxt_name.compare(q_name) == 0) {
                read_alignment(&ali, buf_ptr, g_idx, read_cs, in_paf);
                if (ali.aligned) {
                    if (j < MAX_ALI_R)
                        *ali_group[j++] = ali;
                }
                in_paf->crt_ptr = buf_ptr;
            } else {
                break;
            }
        } else {
            fill_buf_io(in_paf, in_paf->crt_ptr);
        }
    }
    if (j > 0) {
        preprocess_alignments(ali_group, j, al, g_idx);
        add_alis(q_name, al, ali_group, j);
    }
    for (j = 0; j < MAX_ALI_R; j++)
        delete ali_group[j];
}

uint64_t paf_parse_all_alignments(file_io_t *in_paf, ali_list_t *al, enano_params *p, global_index_t *g_idx) {
    int sz;
    alignment_t ali;
    uint32_t j;
    std::string last_q_name = "";
    alignment_t *ali_group[MAX_ALI_R];
    for (j = 0; j < MAX_ALI_R; j++)
        ali_group[j] = new alignment_t;

    uint64_t sum_targets = 0, num_blk = 0;
    bool read_cs = false; //CHANGE
    j = 0;
    while ((sz = fill_buf_io(in_paf)) > 0) {
        char *buf_ptr = in_paf->buf;
        while (has_next_ali(buf_ptr, in_paf)) {
            char *aux_buf = buf_ptr;
            std::string nxt_name = get_string(aux_buf, '\t');
            if (nxt_name.compare(last_q_name) == 0) {
                read_alignment(&ali, buf_ptr, g_idx, read_cs, in_paf);
                ali.num_blk = num_blk;
                if (j < MAX_ALI_R)
                    *ali_group[j++] = ali;
            } else if (j > 0) {
                // add the ali_group
                sum_targets += preprocess_alignments(ali_group, j, al, g_idx);
                add_alis(last_q_name, al, ali_group, j);
                j = 0;
            } else {
                read_alignment(&ali, buf_ptr, g_idx, read_cs, in_paf);
                ali.num_blk = num_blk;
                if (ali.aligned) {
                    *ali_group[j++] = ali;
                    assert(j < MAX_ALI_R);
                    last_q_name = ali.q_name;
                }
            }
        }
        num_blk++;
        al->blk_lims.push_back(buf_ptr - in_paf->buf);
        move_remainder_io(in_paf, buf_ptr);
    }
    if (j > 0) {
        sum_targets += preprocess_alignments(ali_group, j, al, g_idx);
        add_alis(last_q_name, al, ali_group, j);
    }
    for (j = 0; j < MAX_ALI_R; j++)
        delete ali_group[j];

    return sum_targets;
}
