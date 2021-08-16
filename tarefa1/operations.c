#include <float.h>
#include <stdio.h>
#include <math.h>

#include "types.h"
#include "manipulation.h"

int sumInterval(interval firstInterval, interval secondInterval, interval *result) {
    Float_t minimum, maximum;

    minimum.f = firstInterval.minimum + secondInterval.minimum;
    maximum.f = firstInterval.maximum + secondInterval.maximum;

    result->minimum = getMinEqual(minimum);
    result->maximum = getMaxEqual(maximum);
    
    result->isUnitary = isUnitary(result->minimum, result->maximum) ? 1 : 0;

    return isValid(*result);
}

int subtractInterval(interval firstInterval, interval secondInterval, interval *result) {
    Float_t minimum, maximum;

    minimum.f = firstInterval.minimum - secondInterval.maximum;
    maximum.f = firstInterval.maximum - secondInterval.minimum;

    result->minimum = getMinEqual(minimum);
    result->maximum = getMaxEqual(maximum);

    result->isUnitary = isUnitary(result->minimum, result->maximum) ? 1 : 0;

    return isValid(*result);
}

int multiplyInterval(interval firstInterval, interval secondInterval, interval *result) {
    Float_t minimum, maximum;
    float mult1, mult2, mult3, mult4;

    mult1 = firstInterval.minimum * secondInterval.minimum;
    mult2 = firstInterval.minimum * secondInterval.maximum;
    mult3 = firstInterval.maximum * secondInterval.minimum;
    mult4 = firstInterval.maximum * secondInterval.maximum;

    minimum.f = defineMinimum(mult1, mult2, mult3, mult4);
    maximum.f = defineMaximum(mult1, mult2, mult3, mult4);

    result->minimum = getMinEqual(minimum);
    result->maximum = getMaxEqual(maximum);

    result->isUnitary = isUnitary(result->minimum, result->maximum) ? 1 : 0;

    return isValid(*result);
}

int divideInterval(interval firstInterval, interval secondInterval, interval *result) {
    Float_t minimum, maximum;

    if(!isInRange(secondInterval.minimum, secondInterval.maximum, 0)) {
        float mult1, mult2, mult3, mult4;

        mult1 = firstInterval.minimum * (1.0 / secondInterval.minimum);
        mult2 = firstInterval.minimum * (1.0 / secondInterval.maximum);
        mult3 = firstInterval.maximum * (1.0 / secondInterval.minimum);
        mult4 = firstInterval.maximum * (1.0 / secondInterval.maximum);

        minimum.f = defineMinimum(mult1, mult2, mult3, mult4);
        maximum.f = defineMaximum(mult1, mult2, mult3, mult4);

        result->minimum = getMinEqual(minimum);
        result->maximum = getMaxEqual(maximum);
    }
    else {
        result->minimum = -INFINITY;
        result->maximum = INFINITY;
    }
    
    result->isUnitary = isUnitary(result->minimum, result->maximum) ? 1 : 0;
    return isValid(*result);
}