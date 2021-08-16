// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#ifndef __SISTEMA_LINEAR_H__
#define __SISTEMA_LINEAR_H__

#include <math.h>
#include <string.h>
#include <stdlib.h>

/*!
    \brief Estrutura do sistema linear utilizado.

    \param n Ordem do sistema linear.
    \param U Matriz de coeficientes do sistema que será montada de forma a satisfazer os requerimentos
            dos métodos de interpolação polinomial e ajuste de curvas. Será alterada durante a triangularização
            como forma de representar a matriz U do Método de Fatoração LU.
    \param L Vetor que representa a matriz formada pelos coeficientes obtidos durante a execução do Método
            de Gauss. Não inclui os valores zeros como forma de poupar o espaço utilizado.
    \param b Vetor de termos independentes do sistema que será montado de forma a satisfazer os requerimentos
            dos métodos de interpolação polinomial e ajuste de curvas.
*/
typedef struct {
    unsigned int n;
    real_t *U, *L, *b;
} Sistema_t;

Sistema_t* alocaSistema(unsigned int n);
void liberaSistema(Sistema_t *sistema);

void preencheSistemaInterpolacao(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas);
void preencheSistemaAjuste(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas);
void preencheSistemaAjusteOtimizada(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas);
unsigned int encontraMax(real_t *matriz, unsigned int n, unsigned int col);
void trocaLinha(Sistema_t *sistema, unsigned int n, unsigned int i, unsigned int iPivo);
int triangularizacaoSistema(Sistema_t *sistema);
int triangularizacaoSistemaOtimizada(Sistema_t *sistema);
void calculaY(real_t *L, real_t *y, real_t *b, unsigned int n);
int calculaX(real_t *U, real_t *x, real_t *y, unsigned int n);
int resolveSistema(Sistema_t *sistema, real_t *resultado);

#endif // __SISTEMA_LINEAR_H__