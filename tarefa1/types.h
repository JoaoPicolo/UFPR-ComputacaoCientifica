#ifndef __TYPES__
#define __TYPES__

#include <stdint.h>

typedef struct {
   float  minimum, maximum;
   int isUnitary;
} interval;

typedef struct {
   long int  firstIndex, secondIndex;
   char op;
} operation;

typedef union
{
    int32_t i;
    float f;
    struct
    {   // Bitfields for exploration.
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    } parts;
} Float_t;

#endif