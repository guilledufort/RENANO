#ifndef CLR_H
#define CLR_H
/*
 * Note it is up to the calling code to ensure that no overruns on input and
 * output buffers occur.
 *
 * Call the input() and output() functions to set and query the current
 * buffer locations.
 */
#include <stdlib.h>

#define DO(n) for (int _ = 0; _ < n; _++)
#define TOP (1 << 24)

#define CodecName "CLR"

typedef unsigned char uc;

class RangeCoder {
    uint64_t low;
    uint32_t range, code;

    uint8_t bits_left;

   public:
    uc *in_buf;
    uc *out_buf;

    void input(char *in) { out_buf = in_buf = (uc *)in; }

    void output(char *out) { in_buf = out_buf = (uc *)out; }

    char *input(void) { return (char *)in_buf; }

    char *output(void) { return (char *)out_buf; }

    int size_out(void) { return out_buf - in_buf; }

    int size_in(void) { return in_buf - out_buf; }

    void StartEncode(void) {
        low = 0;
        range = (uint32_t)-1;
    }

    void StartEncodeBinary(void) {
        bits_left = 8;
        *out_buf = 0;
    }

    void StartDecode(void) {
        low = 0;
        range = (uint32_t)-1;
        DO(8)
        code = (code << 8) | *in_buf++;
    }

    void StartDecodeBinary(void) {
        bits_left = 8;
    }

    void FinishEncode(void) {
        DO(8)
        (*out_buf++ = low >> 56), low <<= 8;
    }

    void FinishEncodeBinary(void) {
        if (bits_left < 8)
            out_buf++;
    }

    void FinishDecode(void) {}

    void FinishDecodeBinary(void) {}

    void Encode(uint32_t cumFreq, uint32_t freq, uint32_t totFreq) {
        low += cumFreq * (range /= totFreq);
        range *= freq;

        if (cumFreq + freq > totFreq)
            abort();

        while (range < TOP) {
            // range = 0x00ffffff..
            // low/high may be matching
            //       eg 88332211/88342211 (range 00010000)
            // or differing
            //       eg 88ff2211/89002211 (range 00010000)
            //
            // If the latter, we need to reduce range down
            // such that high=88ffffff.
            // Eg. top-1      == 00ffffff
            //     low|top-1  == 88ffffff
            //     ...-low    == 0000ddee
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint32_t(low) | (TOP - 1)) - uint32_t(low));
            *out_buf++ = low >> 56, range <<= 8, low <<= 8;
        }
    }

    void EncodeBinary(uint8_t sym) {
        *out_buf = *out_buf + (sym << --bits_left);
        if (bits_left == 0) {
            bits_left = 8;
            ++out_buf;
            *out_buf = 0;
        }
    }

    uint32_t GetFreq(uint32_t totFreq) {
        return code / (range /= totFreq);
    }

    void Decode(uint32_t cumFreq, uint32_t freq) {
        uint32_t temp = cumFreq * range;
        low += temp;
        code -= temp;
        range *= freq;

        while (range < TOP) {
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint32_t(low) | (TOP - 1)) - uint32_t(low));
            code = (code << 8) | *in_buf++, range <<= 8, low <<= 8;
        }
    }

    int8_t DecodeBinary() {
        int8_t ret = (*in_buf & (1 << --bits_left)) > 0;
        if (bits_left == 0) {
            bits_left = 8;
            in_buf++;
        }
        return ret;
    }
};

#endif
