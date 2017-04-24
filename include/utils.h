#ifndef UTILS_H
#define UTILS_H

#include <complex>

typedef double ffast_real;
typedef std::complex<double> ffast_complex;

// Perform a mod operation and make sure that the result is positive
inline int positiveMod(int a, int b)
{
  const int result = a % b;

  return result >= 0 ? result : result + b;
}

// useful in WHT: maps shorter binary number to pattern indicated by longer
// required: numOnes(longer) == lenS
// [] denotes base 2 representation: 5 = [101]
// mapToInt([101],3,[1010010],7) = [1000010]
inline int mapToInt(int shorter, int lenS, int longer, int lenL) {
    int shortCounter = 0;
    int ret = 0;
    for (int longCounter =0; longCounter < lenL; longCounter++) {
        if (longer & (1 << longCounter))
            if (shorter & (1 << shortCounter++))
                ret += (1 << longCounter);
    }
    return ret;
}

// Black magic (http://stackoverflow.com/questions/14555607/explanation-required-number-of-bits-set-in-a-number)
inline int numOnes(int v) {
    v = v - ((v>>1) & 0x55555555);
    v = (v & 0x33333333) + ((v>>2) & 0x33333333);
    return ((v + (v>>4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

#endif // UTILS_H
