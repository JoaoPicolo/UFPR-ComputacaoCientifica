#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    SistLinear_t *sistemaLinear;
    real_t *vetSolucao, *residuo, normaL2;
    double *tempo;
    int iteracoes;
    int sistemaCount = 0;
    char proximaLinha;

    while((proximaLinha = getchar()) != EOF) {
        if(proximaLinha != '\n') {
            ungetc(proximaLinha, stdin);

            if((sistemaLinear = lerSistLinear()) != NULL) {
                vetSolucao = malloc(sistemaLinear->n * sizeof(real_t));
                residuo = malloc(sistemaLinear->n * sizeof(real_t));
                tempo = malloc(sizeof(double));

                if(vetSolucao != NULL && residuo != NULL && tempo != NULL) {
                    sistemaCount++;
                    printf("\n\n***** Sistema %d --> n = %d, erro: %f\n", sistemaCount, sistemaLinear->n, sistemaLinear->erro);

                    // Gauss
                    if((iteracoes = eliminacaoGauss(sistemaLinear, vetSolucao, tempo)) != -2) { // -2 indica que o sistema não tem solução
                        if(iteracoes != -1) { // -1 indica que houve erro durante a execução e Gauss
                            printf("===> Eliminação Gauss: %f ms\n", *tempo);
                            printf("  --> X: ");
                            prnVetor(vetSolucao, sistemaLinear->n);
                            if((normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo)) >= 0) {
                                printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                if(normaL2 > 5.0) {
                                    iteracoes = refinamento(sistemaLinear, vetSolucao, tempo);
                                    printf("\n===> Refinamento: %f ms --> %d iteracoes\n", *tempo, iteracoes);
                                    printf("  --> X: ");
                                    prnVetor(vetSolucao, sistemaLinear->n);
                                    normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo);
                                    printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                }
                            }
                        }

                        // Jacobi
                        if((iteracoes = gaussJacobi(sistemaLinear, vetSolucao, tempo)) >= 0) {
                            printf("\n===> Jacobi: %f ms --> %d iteracoes\n", *tempo, iteracoes);
                            printf("  --> X: ");
                            prnVetor(vetSolucao, sistemaLinear->n);
                            if((normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo)) >= 0) {
                                printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                if(normaL2 > 5.0) {
                                    iteracoes = refinamento(sistemaLinear, vetSolucao, tempo);
                                    printf("\n===> Refinamento: %f ms --> %d iteracoes\n", *tempo, iteracoes);
                                    printf("  --> X: ");
                                    prnVetor(vetSolucao, sistemaLinear->n);
                                    normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo);
                                    printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                }
                            }
                        }

                        // Gauss - Seidel
                        if((iteracoes = gaussSeidel(sistemaLinear, vetSolucao, tempo)) >= 0) {
                            printf("\n===> Gauss-Seidel: %f ms --> %d iteracoes\n", *tempo, iteracoes);
                            printf("  --> X: ");
                            prnVetor(vetSolucao, sistemaLinear->n);
                            if((normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo)) >= 0) {
                                printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                if(normaL2 > 5.0) {
                                    iteracoes = refinamento(sistemaLinear, vetSolucao, tempo);
                                    printf("\n===> Refinamento: %f ms --> %d iteracoes\n", *tempo, iteracoes);
                                    printf("  --> X: ");
                                    prnVetor(vetSolucao, sistemaLinear->n);
                                    normaL2 = normaL2Residuo(sistemaLinear, vetSolucao, residuo);
                                    printf("  --> Norma L2 do residuo: %1.8e\n", normaL2);
                                }
                            }
                        }
                    }

                    // Libera alocações antes da proxima iteração
                    free(tempo);
                    free(residuo);
                    free(vetSolucao);
                    liberaSistLinear(sistemaLinear);
                }
                else {
                    fprintf(stderr, "%s", "Erro ao alocar variáveis auxiliares. Desconsiderando entrada...\n"); 
                }
            }
            else {
                fprintf(stderr, "%s", "Erro ao alocar sistema. Desconsiderando entrada...\n"); 
            }
        }
    }
}

