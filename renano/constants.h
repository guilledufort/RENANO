// MIT License

// Copyright (c) 2021 Guillermo Dufort y Álvarez

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

#ifndef CONSTANTS_H
#define CONSTANTS_H

// #define DEBUG_STORE_REF

//Header size
#define H_SZ 11

#define MAX_INS_LEN ((1 << 16) - 1)
#define MAX_MATCH_LEN ((1 << 16) - 1)
#define MAX_SKIP_LEN ((1 << 16) - 1)

#define BLK_SIZE 10000000  //10 MB
#define MEDIUM_BLK_SIZE 1000000  //1 MB
#define SMALL_BLK_SIZE 100000  //0.1 MB
#define MAX_READS 10000000
#define MAX_CHROM 100000
// Size of lists to store cs operations
#define MAX_ALIS_L 2000
// Average ali amount per read
#define AVG_ALI_PER_R 2
// Maximum number of aligments per read.
#define MAX_ALI_R 10
// Minimum length of alignment
#define MIN_ALI_LEN 30
//Minimum match operation length
#define MIN_MATCH_LEN 2
//Minium coverage of overlapping region for storing reference sequence
#define MIN_COVERAGE 1
#define MIN_ALI_OVLP 0

#define MATCH 0
#define INSER 1
#define FOR 0
#define REV 1

//Default parameters
#define DEFAULT_K_LEVEL 7
#define DEFAULT_L_LEVEL 6
#define DEFAULT_THREADS_NUM 8
#define DEFAULT_BLK_UPD_THRESH 32
#define DEFAULT_BLK_UPD_FREQ 4

//#define __TIMING__
//#define __ORDER_SYMBOLS__
//#define __STATS_SNPS__

#define REF_ALI 1
#define STORE_REF_ALI 2

#define MAJOR_VERS 1
#define MINOR_VERS 0

#define SEQ_ALPHA_SIZE 5
static const unsigned int pow5[] = {
    1,
    5,
    25,
    125,
    625,
    3125,
    15625,
    78125,
    390625,
    1953125,
    9765625,
    48828125,
    244140625,
    1220703125};

static const unsigned char IS_BC[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                        0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                        0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static const unsigned char REV_BC[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84, 0, 71, 0, 0, 0, 67,
                                          0, 0, 0, 0, 0, 78, 0, 0, 0, 0, 0, 0, 65, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0};

// Sequence table lookups ACGTN->0..4
static const unsigned char L[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                                    0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                                    0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//#define __GLOBAL_STATS__
//#define __DEBUG_LOG__
//#define __CONTEXT_STATS__

#define ROUND_UP(x) (x == 0 ? 0 : 1 << ((x)-1))
/* Quantize quality scores q' = ((q + ROUND_UP(QUANT_Q)) >> QUANT_Q) */

#define DIV_ROUND(q, SHIFT) ((q + ROUND_UP(SHIFT)) >> SHIFT)

/* Keep as a power of 2 */
#define QMAX 128
#define A_LOG 2

#define QUANT_Q 3

#define LOG_AVGS 4
#define LOG_AVGS_ERRS 2

#define DIF_CTX_LOG 3
#define DIF_CANT (1 << DIF_CTX_LOG)

#define QUANT_D_LOG 5
#define QUANT_D_CANT (1 << QUANT_D_LOG)
#define QUANT_D_MASK (QUANT_D_CANT - 1)
#define QUANT_D_MAX QUANT_D_MASK

#define Q_LOG (7 - QUANT_Q)
#define TOTAL_Q_LOG (DIF_CTX_LOG + Q_LOG)

#define Q_CTX (1 << TOTAL_Q_LOG)

#define BC_CTX (1 << (LOG_AVGS + LOG_AVGS_ERRS))

#define CTX_CNT (BC_CTX * Q_CTX)

#define Q_LOG_CANT (1 << Q_LOG)

#define AVG_SHIFT 4
#define TOTAL_ERR_SHIFT 4

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define ABS(N) ((N < 0) ? (-N) : (N))

#endif
