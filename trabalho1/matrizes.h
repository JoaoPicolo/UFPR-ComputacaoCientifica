// ###################################
// Gabriel Marckzuk Thá - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#ifndef __MATRIZES_H__
#define __MATRIZES_H__

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"

/*!
    \brief Estrutura de matriz a guardar os dados

    \param n Ordem da matriz.
    \param ordemLinhas Vetor que armazena em seu i-ésimo indíce a
                posição do valor 1 na i-ésima coluna da matriz identidade
    \param A Matriz lida da entrada padrão (stdin).
    \param AI Matri inversa obtida após a fatoração LU.
    \param L Vetor que em sua i-ésima posição possui um ponteiro para os elementos
                da i-ésima linha da matriz L da fatoração LU. Não inclui valores zeros.
    \param U Vetor que em sua i-ésima posição possui um ponteiro para os elementos
                da i-ésima linha da matriz U da fatoração LU. Não inclui valores zeros.
    \param normas Vetor contendo em sua i-ésima posição a norma L2 do resíduo da matriz A
                com a i-ésima coluna de AI em relação a i-ésima coluna da matriz identidade.
*/
typedef struct {
    unsigned int n;
    int *ordemLinhas;
    real_t **A, **AI;
    real_t **L, **U;
    real_t *normas;
} Matriz_t;

Matriz_t* alocaMatriz(unsigned int n);
void liberaMatriz(Matriz_t *matriz);
real_t** copiaMatriz(real_t **matriz, unsigned int n);
void liberaCopia(real_t **matriz, unsigned int n);

Matriz_t* lerMatriz();

unsigned int encontraMax(real_t **matriz, unsigned int n, unsigned int col);
void trocaLinha(real_t **matriz, real_t **L, int *ordemLinhas, 
    unsigned int n, unsigned int i, unsigned int iPivo);

#endif // __MATRIZES__
