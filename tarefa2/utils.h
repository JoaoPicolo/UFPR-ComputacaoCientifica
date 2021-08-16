#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "SistemasLineares.h"

typedef struct {
  int tipo; // 0 se troca de linha, 1 se operacao matematica
  int linha1, linha2; // linhas
  real_t m; // valor de m
} RegrasAplicadas_t;

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/
double timestamp(void);

int linhaVazia(char *linha);

real_t** copiaMatrizCoeficientes(real_t **matCoeficientes, unsigned int n);

real_t* copiaVetorTermos(real_t *vetTermos, unsigned int n);

void liberaCopiasLocais(real_t **matCoeficientes, real_t *vetTermos, unsigned int n);

unsigned int encontraMax(real_t **matCoeficientes, unsigned int n, unsigned int col);

void trocaLinha(real_t **matCoeficientes, real_t *vetTermos, unsigned int n, unsigned int i, unsigned int iPivo);

void retrosubstituicao(real_t **matCoeficientes, real_t *x, real_t *vetTermos, unsigned int n);

void calculaResiduo(real_t **matCoeficientes, real_t *xRefinado, real_t *vetTermos, real_t *residuo, unsigned int n);

void inicializaVetAproximacao(real_t *x, unsigned int n);

int menorQueErro(real_t *vetInicial, real_t *x, real_t erro, unsigned int n, int it);

int sistemaConverge(real_t **matCoeficientes, unsigned int n);

void aplicaRegrasResiduo(real_t *residuo, RegrasAplicadas_t *regras, int n);

int possuiSolucao(real_t **matCoeficientes, real_t *vetTermos, int n);

#endif // __UTILS_H__

