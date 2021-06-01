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

#include <fstream>
#include <iostream>

#include "enc_algs.h"

uint32_t tot_a0 = 0, equal_a0 = 0, tot_0 = 0, equal_0 = 0;

/* --------------------------------------------------------------------------
 * Main program entry.
 */
static void usage(int err) {
    printf("Renano v%d.%d Author Guillermo Dufort y Alvarez, 2020-2021\n\n",
           MAJOR_VERS, MINOR_VERS);

    printf("COMPRESSION: \n");

    printf("> Without reference:\n" 
            "\trenano [options] [input_file [output_file]]\n\n");

    printf("> With reference:\n"
            "\trenano [options] -r [ref_file [paf_file]] [input_file [output_file]]\n\n");

    printf("> With reference and making decompression independent of the reference:\n"
            "\trenano [options] -s [ref_file [paf_file]] [input_file [output_file]]\n\n");

    printf("COMPRESSION OPTIONS: \n");
    printf(
        "\t-k <length>    Base call sequence context length. Default is 7 (max "
        "13).\n\n");
    printf(
        "\t-l <lenght>    Length of the DNA neighborhood sequence. Default is "
        "6.\n\n");
    printf(
        "\t-t <num>       Maximum number of threads allowed to use by the "
        "compressor. Default is 8.\n\n");

    printf("DECOMPRESSION: \n");

    printf("> Without reference:\n"   
            "\trenano -d [options] foo.enano foo.fastq\n\n");

    printf("> With reference:\n"   
            "\trenano -d [options] -r [ref_file] foo.enano foo.fastq\n\n");

    printf("DECOMPRESSION OPTIONS: \n");
    printf(
        "\t-t <num>       Maximum number of threads allowed to use by the "
        "decompressor. Default is 8.\n\n");


    printf("CREDITS :\n");
    printf(
        "The methods used for encoding the reads identifiers, and to model "
        "frequency counters, \n");
    printf(
        "are the ones proposed by James Bonefield in FQZComp, with some "
        "modifications.\n");
    printf("The range coder is derived from Eugene Shelwien.\n");
    printf("The kseq library used to parse FASTA files is authored by Heng Li.\n\n");

    exit(err);
}

bool get_duplicate_name(char *f_in) {
    std::ifstream File(f_in);
    std::vector<std::string> Line(3);
    int counter = 0;

    while (counter < 3 && std::getline(File, Line[counter++]))
        ;
    File.close();

    if (counter < 2) {
        printf("ERROR: empty file. \n");
        exit(1);
    }

    return Line[2].size() > 1;
}

void parse_command_line_args(enano_params &p, int argc, char **argv) {
    int opt;
    while ((opt = getopt(argc, argv, "v:g:hsdark:l:t:cb:")) != -1) {
        switch (opt) {
            case 'h':
                usage(0);
                break;
            case 'd':
                p.decompress = 1;
                break;
            case 'r':
                p.aligned = REF_ALI;
                break;
            case 's':
                p.aligned = STORE_REF_ALI;
                break;
            case 'k':
                char *end;
                p.klevel = strtol(optarg, &end, H_SZ);
                if (p.klevel > 13) usage(1);
                break;
            case 'l':
                p.llevel = atoi(optarg);
                break;
            case 't':
                p.num_threads = atoi(optarg);
                break;
            case 'b':
                p.blk_upd_thresh = atoi(optarg);
                break;
            case 'c':
                p.max_comp = true;
                break;
            case 'g':
                p.min_cvg = std::atof(optarg);
                break;
            case 'v':
                p.min_ovlp = std::atof(optarg);
                break;
            default:
                usage(1);
        }
    }

    int open_flag = O_RDONLY;

    if (p.aligned) {
        if (((p.decompress == 1) && (argc - optind != 3)) || ((p.decompress == 0) && (argc - optind != 4))) {
            printf("ERROR: Not enough arguments for REF aligned de/compression.\n");
            usage(1);
        }
        p.ref_name = argv[optind++];
        
        if (!p.decompress) {
            p.paf_name = argv[optind];
            if ((p.in_paf = open(argv[optind], open_flag)) == -1) {
                perror(argv[optind]);
                exit(1);
            }
            optind++;
        }
    } else {
        if (argc - optind != 2) {
            printf("ERROR: Either input or output files are missing.\n");
            usage(1);
        }
    }

    if (optind != argc) {
        if (!p.decompress)
            p.duplicate_names = get_duplicate_name(argv[optind]);
        p.fastq_name = argv[optind];
        if ((p.in_fd = open(argv[optind], open_flag)) == -1) {
            perror(argv[optind]);
            exit(1);
        }
        optind++;
    }

    if (optind != argc) {
        p.out_fd = open(argv[optind], O_RDWR | O_CREAT | O_TRUNC, 0666);
        if (p.out_fd == -1) {
            perror(argv[optind]);
            exit(1);
        }
        optind++;
    }
    omp_set_num_threads(p.num_threads);
}

void write_header(enano_params &p) {
    unsigned char magic[H_SZ] = {'.',
                                 'e',
                                 'n',
                                 'a',
                                 MAJOR_VERS,
                                 (unsigned char)p.klevel,
                                 (unsigned char)p.llevel,
                                 (unsigned char)p.max_comp,
                                 (unsigned char)p.blk_upd_thresh,
                                 (unsigned char)p.aligned,
                                 (unsigned char)p.duplicate_names};

    if (H_SZ != write(p.out_fd, magic, H_SZ)) {
        printf("Abort: truncated write.\n");
        exit(1);
    }
}

void read_header(enano_params &p) {
    unsigned char magic[H_SZ];

    /* Check magic number */
    if (H_SZ != read(p.in_fd, magic, H_SZ)) {
        printf("Abort: truncated read.\n");
        exit(1);
    }
    if (memcmp(".ena", magic, 4) != 0) {
        printf("Unrecognised file format.\n");
        exit(1);
    }
    if (magic[4] != MAJOR_VERS) {
        printf("Unsupported file format version %d.%d\n", magic[4], magic[5]);
        exit(1);
    }

    p.klevel = magic[5] & 0x0f;
    if (p.klevel > 13 || p.klevel < 1) {
        printf("Unexpected quality compression level %d\n", p.klevel);
        exit(1);
    }

    p.llevel = magic[6] & 0x0f;
    p.max_comp = magic[7] & 0x0f;
    p.blk_upd_thresh = magic[8] & 0xff;
    p.aligned = magic[9] & 0xff;
    p.duplicate_names = magic[10] & 0xff;
}

int main(int argc, char **argv) {
    double start_time = omp_get_wtime();
    int res = 0;
    enano_params p;

    /* Initialise and parse command line arguments */
    initialize_enano_params(p);
    parse_command_line_args(p, argc, argv);

    if (p.decompress) {
        read_header(p);

        printf("Parameters - k: %d, l: %d, b: %d \n", p.klevel, p.llevel,
               p.blk_upd_thresh);

        /* If p.aligned we intialize the necessary structures */
        init_g_idx(p);

        if (p.max_comp)
            res = decode_st(p);
        else
            res = decode(p);
    } else {
        /* If p.aligned we intialize the necessary structures */
        init_g_idx(p);

        write_header(p);

        printf("Parameters - k: %d, l: %d, b: %d \n", p.klevel, p.llevel,
               p.blk_upd_thresh);

        if (p.max_comp)
            res = encode_st(p);
        else
            res = encode(p);
    }
    free_g_idx(p.g_idx);

    printf("Total time: %.2f s\n", (double)(omp_get_wtime() - start_time));
    return res;
}
