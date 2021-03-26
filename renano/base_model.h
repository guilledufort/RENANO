#ifndef BASE_MODEL_H
#define BASE_MODEL_H
/*
 * A fixed alphabet encoder for symbols 0 to 4 inclusive.
 * We have no escape symbol or resorting of data. It simply
 * accumulates stats and encodes in proportion.
 *
 * The intention is to use this per column for encoding
 * bases that differ to the consensus.
 */

// Small enough to not overflow uint16_t even after +STEP
#ifdef WSIZ
#undef WSIZ
#endif

#define M4(a) ((a)[0] > (a)[1] ? ((a)[0] > (a)[2] ? ((a)[0] > (a)[3] ? (a)[0] : (a)[3]) : ((a)[2] > (a)[3] ? (a)[2] : (a)[3])) : ((a)[1] > (a)[2] ? ((a)[1] > (a)[3] ? (a)[1] : (a)[3]) : ((a)[2] > (a)[3] ? (a)[2] : (a)[3])))

template <typename st_t>
struct BASE_MODEL {
    enum { STEP = sizeof(st_t) == 1 ? 1 : 8 };
    enum { WSIZ = (1 << 8 * sizeof(st_t)) - 2 * STEP };

    BASE_MODEL();
    BASE_MODEL(int *start);
    void reset();
    void reset(int *start);
    inline void encodeSymbol(RangeCoder *rc, uint sym, bool update);
    inline void encodeSymbolNotB(RangeCoder *rc, uint sym, uint not_b, bool update);
    inline void encodeSymbolRegular(RangeCoder *rc, uint sym);
    inline void encodeSymbolNoUpdate(RangeCoder *rc, uint sym);
    inline void updateModel(uint sym);
    inline void mix(BASE_MODEL *cm);
    inline void mix_array(void **models, uc len);
    inline uint decodeSymbol(RangeCoder *rc, bool update);
    inline uint decodeSymbolNotB(RangeCoder *rc, uint not_b, bool update);
    inline uint decodeSymbolRegular(RangeCoder *rc);
    inline uint decodeSymbolNoUpdate(RangeCoder *rc);
    inline uint getTopSym(void);
    inline uint getSummFreq(void);

    void rescaleRare();
    st_t Stats[5];
};

template <typename st_t>
BASE_MODEL<st_t>::BASE_MODEL() {
    reset();
}

template <typename st_t>
BASE_MODEL<st_t>::BASE_MODEL(int *start) {
    for (int i = 0; i < 5; i++) {
        Stats[i] = start[i];
    }
}

template <typename st_t>
void BASE_MODEL<st_t>::reset() {
    for (int i = 0; i < 4; i++)
        Stats[i] = 10 * STEP;
    Stats[4] = 1;
}

template <typename st_t>
void BASE_MODEL<st_t>::reset(int *start) {
    for (int i = 0; i < 5; i++) {
        Stats[i] = start[i];
    }
}

template <typename st_t>
void BASE_MODEL<st_t>::rescaleRare() {
    Stats[0] -= (Stats[0] >> 1);
    Stats[1] -= (Stats[1] >> 1);
    Stats[2] -= (Stats[2] >> 1);
    Stats[3] -= (Stats[3] >> 1);
    Stats[4] -= (Stats[4] >> 1);
}

#define Encode256 Encode
#define GetFreq256 GetFreq

template <typename st_t>
inline void BASE_MODEL<st_t>::encodeSymbol(RangeCoder *rc, uint sym, bool update) {
    if (update) {
        encodeSymbolRegular(rc, sym);
    } else {
        encodeSymbolNoUpdate(rc, sym);
    }
    return;
}

template <typename st_t>
inline void BASE_MODEL<st_t>::encodeSymbolRegular(RangeCoder *rc, uint sym) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    if (SummFreq >= WSIZ) {
        rescaleRare();
        SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    }

    switch (sym) {
        case 0:
            rc->Encode256(0, Stats[0], SummFreq);
            Stats[0] += STEP;
            break;
        case 1:
            rc->Encode256(Stats[0], Stats[1], SummFreq);
            Stats[1] += STEP;
            break;
        case 2:
            rc->Encode256(Stats[0] + Stats[1], Stats[2], SummFreq);
            Stats[2] += STEP;
            break;
        case 3:
            rc->Encode256((Stats[0] + Stats[1]) + Stats[2], Stats[3], SummFreq);
            Stats[3] += STEP;
            break;
        case 4:
            rc->Encode256((Stats[0] + Stats[1]) + (Stats[2] + Stats[3]), Stats[4], SummFreq);
            Stats[4] += STEP;
            break;
    }

    return;
}

template <typename st_t>
inline void BASE_MODEL<st_t>::encodeSymbolNoUpdate(RangeCoder *rc, uint sym) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];

    switch (sym) {
        case 0:
            rc->Encode256(0, Stats[0], SummFreq);
            break;
        case 1:
            rc->Encode256(Stats[0], Stats[1], SummFreq);
            break;
        case 2:
            rc->Encode256(Stats[0] + Stats[1], Stats[2], SummFreq);
            break;
        case 3:
            rc->Encode256((Stats[0] + Stats[1]) + Stats[2], Stats[3], SummFreq);
            break;
        case 4:
            rc->Encode256((Stats[0] + Stats[1]) + (Stats[2] + Stats[3]), Stats[4], SummFreq);
            break;
    }

    return;
}

template <typename st_t>
inline void BASE_MODEL<st_t>::encodeSymbolNotB(RangeCoder *rc, uint sym, uint not_b, bool update) {
    //We cant modify Stats directly because it is a shared resource
    st_t cpyStats[5];
    memcpy(cpyStats, Stats, 5);

    cpyStats[not_b] = 3 * STEP;
    uint SummFreq_not_b = (cpyStats[0] + cpyStats[1]) + (cpyStats[2] + cpyStats[3]) + cpyStats[4];

    switch (sym) {
        case 0:
            rc->Encode256(0, cpyStats[0], SummFreq_not_b);
            break;
        case 1:
            rc->Encode256(cpyStats[0], cpyStats[1], SummFreq_not_b);
            break;
        case 2:
            rc->Encode256(cpyStats[0] + cpyStats[1], cpyStats[2], SummFreq_not_b);
            break;
        case 3:
            rc->Encode256((cpyStats[0] + cpyStats[1]) + cpyStats[2], cpyStats[3], SummFreq_not_b);
            break;
        case 4:
            rc->Encode256((cpyStats[0] + cpyStats[1]) + (cpyStats[2] + cpyStats[3]), cpyStats[4], SummFreq_not_b);
            break;
    }

    if (update) 
        Stats[sym] += STEP;

    return;
}


template <typename st_t>
inline void BASE_MODEL<st_t>::updateModel(uint sym) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    if (SummFreq >= WSIZ) {
        rescaleRare();
    }

    /* known symbol */
    Stats[sym] += STEP;
}

template <typename st_t>
inline void BASE_MODEL<st_t>::mix(BASE_MODEL<st_t> *m) {
    for (int i = 0; i < 5; i++) {
        Stats[i] = MAX((Stats[i] + m->Stats[i]) >> 1, 1);
    }
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    if (SummFreq >= WSIZ) {
        rescaleRare();
    }
}

template <typename st_t>
inline void BASE_MODEL<st_t>::mix_array(void **models, uc len) {
    for (int i = 0; i < 5; i++) {
        uint fmix = 0;
        for (uc c = 0; c < len; c++) {
            fmix += ((BASE_MODEL<st_t> *)models[c])->Stats[i];
        }
        fmix = round((double)fmix / len);
        Stats[i] = fmix;
    }
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    if (SummFreq >= WSIZ) {
        rescaleRare();
    }
}

/*
 * Returns the bias of the best symbol compared to all other symbols.
 * This is a measure of how well adapted this model thinks it is to the
 * incoming probabilities.
 */
template <typename st_t>
inline uint BASE_MODEL<st_t>::getTopSym(void) {
    uint m4 = M4(Stats);
    return m4 >= Stats[4] ? m4 : Stats[4];
}

template <typename st_t>
inline uint BASE_MODEL<st_t>::getSummFreq(void) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    return SummFreq;
}

template <typename st_t>
inline uint BASE_MODEL<st_t>::decodeSymbol(RangeCoder *rc, bool update) {
    if (update) {
        return decodeSymbolRegular(rc);
    } else {
        return decodeSymbolNoUpdate(rc);
    }
}

template <typename st_t>
inline uint BASE_MODEL<st_t>::decodeSymbolRegular(RangeCoder *rc) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];

    if (SummFreq >= WSIZ) {
        rescaleRare();
        SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];
    }

    uint count = rc->GetFreq256(SummFreq);
    uint HiCount = 0;

    st_t *p = Stats;
    if ((HiCount += *p) > count) {
        rc->Decode(0, *p);
        Stats[0] += STEP;
        return 0;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        Stats[1] += STEP;
        return 1;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        Stats[2] += STEP;
        return 2;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        Stats[3] += STEP;
        return 3;
    }

    rc->Decode(HiCount, Stats[4]);
    Stats[4] += STEP;
    return 4;
}

template <typename st_t>
inline uint BASE_MODEL<st_t>::decodeSymbolNotB(RangeCoder *rc, uint not_b, bool update) {
    //We cant modify Stats directly because it is a shared resource
    st_t cpyStats[5];
    memcpy(cpyStats, Stats, 5);

    cpyStats[not_b] = 3 * STEP;
    uint SummFreq_not_b = (cpyStats[0] + cpyStats[1]) + (cpyStats[2] + cpyStats[3]) + cpyStats[4];

    uint count = rc->GetFreq256(SummFreq_not_b);
    uint HiCount = 0;
    st_t sym = 0;

    st_t *p = cpyStats;
    if ((HiCount += *p) > count) {
        rc->Decode(0, *p);
        sym = 0;
    } else if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        sym = 1;
    } else if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        sym = 2;
    } else if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        sym = 3;
    } else {
        rc->Decode(HiCount, cpyStats[4]);
        sym = 4;
    }

    if (update)
        Stats[sym] += STEP;
    
    return sym;
}

template <typename st_t>
inline uint BASE_MODEL<st_t>::decodeSymbolNoUpdate(RangeCoder *rc) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]) + Stats[4];

    uint count = rc->GetFreq256(SummFreq);
    uint HiCount = 0;

    st_t *p = Stats;
    if ((HiCount += *p) > count) {
        rc->Decode(0, *p);
        return 0;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        return 1;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        return 2;
    }

    if ((HiCount += *++p) > count) {
        rc->Decode(HiCount - *p, *p);
        return 3;
    }

    rc->Decode(HiCount, Stats[4]);
    return 4;
}

#endif
