#ifndef __OPERATIONS__
#define __OPERATIONS__

int sumInterval(interval firstInterval, interval secondInterval, interval *result);
int subtractInterval(interval firstInterval, interval secondInterval, interval *result);
int multiplyInterval(interval firstInterval, interval secondInterval, interval *result);
int divideInterval(interval firstInterval, interval secondInterval, interval *result);

#endif