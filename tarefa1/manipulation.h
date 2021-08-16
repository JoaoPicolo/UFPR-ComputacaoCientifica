#ifndef __MANIPULATION__
#define __MANIPULATION__

float getMinEqual(Float_t num);
float getMaxEqual(Float_t num);
float calculaEpsilonRelativo(float num);
interval getInterval(float num1, float num2);
int isUnitary(float a, float b);
int isValid(interval entry);
int isInRange(float start, float end, float num);
float defineMinimum(float mult1, float mult2, float mult3, float mult4);
float defineMaximum(float mult1, float mult2, float mult3, float mult4);

#endif