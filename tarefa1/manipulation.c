
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#include "types.h"

float getMinEqual(Float_t num) {
    if(num.f == -INFINITY) {
        return num.f;
    }

    if(num.parts.sign) {
        num.i = num.i + 1;
    }
    else {
        num.i = num.i - 1;
    }
    
    return num.f;
}

float getMaxEqual(Float_t num) {
    if(num.f == INFINITY) {
        return num.f;
    }

    if(num.parts.sign) {
        num.i = num.i - 1;
    }
    else {
        num.i = num.i + 1;
    }

    return num.f;
}

float calculaEpsilonRelativo(float num) {
    float epsilon;
    epsilon = num / 2.0f;
    while (num + epsilon / 2.0f > num)
        epsilon /= 2.0f;
    
    return epsilon;
}

interval getInterval(float num1, float num2) {
    interval result;

    // Zero is an exception
    if(fabs(num1) == 0.0f && fabs(num2) == 0.0f) {
        result.minimum = -0.0f;
        result.maximum = 0.0f;

        return result;
    }

    // Changes comparison methods
    // depending on range of number
    if(num1 >= 1.0f && num1 <= 2.0f) {
        if((num1 - num2) <= FLT_EPSILON) {
            result.minimum = num1;
            result.maximum = num2;
        }
        else {
            result.minimum = num2;
            result.maximum = num1;
        }
    }
    else if(num1 <= 1.0f) {
        if(num1 <= num2) {
            result.minimum = num1;
            result.maximum = num2;
        }
        else {
            result.minimum = num2;
            result.maximum = num1;
        }
    }
    else {
        float epsilon = calculaEpsilonRelativo(num1);
        if((num1 - num2) <= epsilon) {
            result.minimum = num1;
            result.maximum = num2;
        }
        else {
            result.minimum = num2;
            result.maximum = num1;
        }
    }

    return result;
}

int isUnitary(float a, float b) {
    Float_t num1, num2;
    
    num1.f = a;
    num2.f = b;
    // Different signs means they do not match.
    if (num1.parts.sign != num2.parts.sign)
    {
        // Check for equality to make sure +0==-0
        if (num1.f == num2.f)
            return 1;
        return 0;
    }
 
    // Find the difference in ULPs.
    int ulpsDiff = abs(num1.i - num2.i);
    
    if (ulpsDiff < 1)
        return 1;

    return 0;
}

int isValid(interval entry) {
    if(
        (entry.minimum > entry.maximum) || (isnan(entry.minimum) || isnan(entry.maximum) ||
        ((entry.minimum == -INFINITY || entry.minimum == INFINITY) && (entry.minimum == entry.maximum)))
    ) {
        return -1;
    }
    else {
        return 1;
    }
}

int isInRange(float start, float end, float num) {
    interval result = getInterval(start, end);

    if(result.minimum <= num && num <= result.maximum) {
        return 1;
    }

    return 0;
}

float defineMinimum(float mult1, float mult2, float mult3, float mult4) {
    return getInterval(getInterval(mult1, mult2).minimum, getInterval(mult3, mult4).minimum).minimum;
}

float defineMaximum(float mult1, float mult2, float mult3, float mult4) {
    return getInterval(getInterval(mult1, mult2).maximum, getInterval(mult3, mult4).maximum).maximum;
}