#ifndef __EQ_PARCIAIS__
#define __EQ_PARCIAIS__

#include "EquacoesOrdinarias.h"

typedef struct {
    real_t *D, *Di, *Di2, *Ds, *Ds2, *B;
    int n, m;
} SLPentadiag_t;

// Equação Diferencial Parcial
typedef struct {
    int n, m; // número de pontos internos na malha
    real_t aX, bX, aY, bY; // intervalo
    real_t (* uX1)(real_t), (* uX2)(real_t), (* uY1)(real_t), (* uY2)(real_t); // Condicoes de contorno
    real_t (* u)(real_t, real_t), (* f)(real_t, real_t); // Resultado
} Edp_t;

real_t funcaoBUX1(real_t x);
real_t funcaoBUX2(real_t x);
real_t funcaoBUY1(real_t y);
real_t funcaoBUY2(real_t y);
real_t funcaoB(real_t hx, real_t hy);
real_t funcaoBF(real_t x, real_t y);

real_t funcaoDUX1(real_t y);
real_t funcaoDUX2(real_t y);
real_t funcaoDUY1(real_t x);
real_t funcaoDUY2(real_t x);
real_t funcaoD(real_t hx, real_t hy);
real_t funcaoDF(real_t x, real_t y);

SLPentadiag_t* alocaSLPentadiag(int n, int m);
void liberaSLPentadiag(SLPentadiag_t *sl);
void imprimeSistemaPentaDiagonal(SLPentadiag_t *sistemaLinear);
void geraPentaDiagonal(Edp_t *edpeq, SLPentadiag_t *sl);

Edp_t* alocaEqParc();
void liberaEqParc(Edp_t * eqParcial);
void preencheEqParc(Edp_t *eqParcial, real_t aX, real_t bX, real_t aY, real_t bY, real_t (* uX1)(real_t),
                    real_t (* uX2)(real_t), real_t (* uY1)(real_t), real_t (* uY2)(real_t),
                    real_t (* u)(real_t, real_t), real_t (* f)(real_t, real_t));
void executaFuncaoPentadiagonal(Edp_t* eqOrdinaria, char item);

real_t normaL2ResiduoPenta(SLPentadiag_t *sl, real_t **solucao);
void calculaResiduoPenta(SLPentadiag_t *sl, real_t **solucao, real_t *residuo);
real_t normaL2ResiduoPenta(SLPentadiag_t *sl, real_t **solucao);
void gaussSeidelPentadiagonal(Edp_t *edpeq, real_t **Y);

void imprimeResultadoPenta(real_t **solucao, int n, int m);

#endif