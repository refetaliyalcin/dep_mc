#include "dSFMT.h"
#include "dSFMT-jump.h"
#include <stdlib.h>

/* copyright James Spencer 2012.
 * New BSD License, see LICENSE.txt for details.
 */

/* Utility (memory-access) functions to enable use of dSFMT from Fortran. */

void* malloc_dsfmt_t(void) {
    /* Allocate sufficient memory for a dSFMT state (ie a variable of type dsfmt_t). */
    return malloc(sizeof(dsfmt_t));
}

void free_dsfmt_t(dsfmt_t* ptr) {
    /* Free memory associated with a dSFMT state (ie a variable of type dsfmt_t). */
    free(ptr);
}



// Utils needed to make a jump at the distance of 10^20
// based on the sample1.c by Mutsuo Saito, Makoto Matsumoto.
// Modifications made by Timo Väisänen, 2017

static const char * jumppoly = "2c6b058dca1fbfb57ebf41e67fec066c8828f2bf"
"9414331d2767fa740aba89987685a79114b3543edbc83476e35fd1e52b8b2436932"
"b9d18a728d8c7d7009a31aabed9bf4646909b8138f3e2a05e611c48dd1ce58a4618"
"3adabf3314da38599af92720efca1535872e7f85ef916b2c1e41dfe8ea764730f6b"
"a2654ab287a55214bcf08cfae416906e4979108d606819b5d9e5b2f11ce028577f6"
"c3788a9c688f1d64f5ae341eb95169824954";

//IN: steps how many 10**20 are taken
void jump_dsfmt_t(dsfmt_t* ptr,const int steps){    
    for (int i = 0; i < steps; i++) {
	    dSFMT_jump(ptr, jumppoly);
    }
}

