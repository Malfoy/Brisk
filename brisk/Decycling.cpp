#include "Decycling.h"
#include <cstdint>
#include <iostream>



DecyclingSet::DecyclingSet(uint k) : k(k), unit(2 * M_PI / k), coef(4 * k, 0) {
    for (size_t i = 4; i < 4 * k; i += 4) {
        coef[i + 1] = sin(unit * (i / 4));
        coef[i + 2] = 2 * coef[i + 1];
        coef[i + 3] = 3 * coef[i + 1];
    }
}



double DecyclingSet::computeR(uint64_t seq) {
    double R = 0;
    for (size_t i = 4 * (k - 1); i > 0; i -= 4) {
        R += coef[i + (seq & 0b11)];
        seq >>= 2;
    }
    return R;
}



bool DecyclingSet::mem(uint64_t seq) {
    if (computeR(seq) > eps) {
        uint64_t rot = ((seq & 0b11) << (2 * (k - 1))) + (seq >> 2);
        return computeR(rot) < eps;
    }
    return false;
}



uint DecyclingSet::memDouble(uint64_t seq) {
    double Rseq = computeR(seq);
    if (Rseq > eps) {
        uint64_t rot = ((seq & 0b11) << (2 * (k - 1))) + (seq >> 2);
        if (computeR(rot) < eps) {
            return 0;
        }
    } else if (Rseq < -eps) {
        uint64_t rot = ((seq & 0b11) << (2 * (k - 1))) + (seq >> 2);
        if (computeR(rot) > -eps) {
            return 1;
        }
    }
    return 2;
}
