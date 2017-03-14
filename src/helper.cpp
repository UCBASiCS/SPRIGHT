#include "helper.h"
#include <stdlib.h>

#include <vector>
#include <iostream>


int mapToInt(int shorter, int lenS, int longer, int lenL) {
    int shortCounter = 0;
    int ret = 0;
    for (int longCounter =0; longCounter < lenL; longCounter++) {
        if (longer & (1 << longCounter))
            if (shorter & (1 << shortCounter++))
                ret += (1 << longCounter);
    }
    return ret;
}

int numOnes(int v) {
    v = v - ((v>>1) & 0x55555555);
    v = (v & 0x33333333) + ((v>>2) & 0x33333333);
    return ((v + (v>>4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

