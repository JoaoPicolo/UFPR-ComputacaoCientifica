#ifndef __COMUNICACAO_IO_H__
#define __COMUNICACAO_IO_H__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "matrizes.h"

void processaLinhaComando(int *pivoteamento, char **saida, int argc, char **argv);
void imprimeResultado(Matriz_t *matriz, double *tTriangularizacao, 
                double *tCalculoY, double* tCalculoX, char *arquivoSaida);

#endif // __COMUNICACAO_IO_H__
