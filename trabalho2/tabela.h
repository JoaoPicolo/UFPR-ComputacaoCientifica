// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#ifndef __TABELA_H__
#define __TABELA_H__

#include <stdlib.h>
#include <string.h>

#include "utils.h"

typedef struct {
    unsigned int n, m;
    real_t *x, *f;
} Tabela_t;

Tabela_t* alocaTabela(unsigned int n, unsigned int m);
void liberaTabela(Tabela_t *tabela);
Tabela_t* lerTabela();

int executaFuncoesPolinomiais(Tabela_t *tabela);

#endif // __TABELA_H__