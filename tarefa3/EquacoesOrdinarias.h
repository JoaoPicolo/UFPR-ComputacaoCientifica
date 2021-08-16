#ifndef __EQ_ORDINARIAS__
#define __EQ_ORDINARIAS__

#include <stdlib.h>
#include <math.h>

#include "utils.h"

typedef double real_t;

typedef struct {
    real_t *D, *Di, *Ds, *B;
    int n;
} SLTridiag_t;

// Equação Diferencial Ordinária
typedef struct {
    int n; // número de pontos internos na malha
    real_t a, b; // intervalo
    real_t ya, yb; // condições contorno
    real_t (* p)(real_t), (* q)(real_t), (* r)(real_t);
} Edo_t;

real_t funcaoAP(real_t x);
real_t funcaoAQ(real_t x);
real_t funcaoAR(real_t x);

real_t funcaoCP(real_t x);
real_t funcaoCQ(real_t x);
real_t funcaoCR(real_t x);

SLTridiag_t* alocaSLTridiag(int n);
void liberaSLTridiag(SLTridiag_t *sl);
void imprimeSistemaTriDiagonal(SLTridiag_t *sistemaLinear);
void geraTriDiagonal(Edo_t *edoeq, SLTridiag_t *sl);

Edo_t* alocaEqOrd();
void liberaEqOrd(Edo_t * eqOrdinaria);
void preencheEqOrd(Edo_t *eqOrdinaria, real_t a, real_t b, real_t ya, real_t yb, real_t (*p)(real_t), 
                    real_t (*q)(real_t), real_t (*r)(real_t));
void executaFuncaoTridiagonal(Edo_t* eqOrdinaria, char item);

void calculaResiduoTri(SLTridiag_t *sl, real_t *solucao, real_t *residuo);
real_t normaL2ResiduoTri(SLTridiag_t *sl, real_t *solucao);
void gaussSeidel(Edo_t *edoeq, real_t *Y);

void imprimeResultadoTri(double *solucao, int n);;

#endif // __EQ_ORDINARIAS__