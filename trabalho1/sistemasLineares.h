// ###################################
// Gabriel Marckzuk Thá - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#ifndef __SISTEMAS_LINEARES_H__
#define __SISTEMAS_LINEARES_H__

#include <string.h>

int eliminacaoGaussSemPivoteamento(Matriz_t *matriz, double *tTotal);
int eliminacaoGaussComPivoteamento(Matriz_t *matriz, double *tTotal);

void calculaY(real_t **L, real_t *y, unsigned int n, int idx);
int calculaX(real_t **U, real_t *x, real_t *y, unsigned int n);
int calculaInversa(Matriz_t *matriz, double *tTotalY, double *tTotalX);

void calculaResiduo(Matriz_t *matriz, real_t *b, real_t *res, int idx);
real_t normaL2Residuo(Matriz_t *matriz, int idx);
void calculaNormas(Matriz_t *matriz);

#endif // __SISTEMAS_LINEARES_H__
