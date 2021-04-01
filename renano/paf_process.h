#ifndef PAF_PROCESS_H
#define PAF_PROCESS_H

#include "alignments.h"

uint64_t paf_parse_all_alignments(file_io_t *in_paf, ali_list_t *al, global_index_t *g_idx);
void paf_parse_read_alignments(std::string q_name, file_io_t *in_paf, ali_list_t *al, global_index_t *g_idx);
void paf_parse_read_cs_strings(std::string q_name, file_io_t *in_paf, ali_list_t *al);

#endif