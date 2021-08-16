// ###################################
// Gabriel Marckzuk Thá - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################


#include "utils.h"
#include "matrizes.h"
#include "comunicacaoIO.h"
#include "sistemasLineares.h"

int main(int argc, char** argv) {
    int pivoteamentoParcial = 0;
    int execucaoGauss, execucaoInversa;
    char *arquivoSaida = NULL;
    char ordemSistema = linhaVazia;
    double *tTriangularizacao = (double*)malloc(sizeof(double));
    double *tCalculoY = (double*)malloc(sizeof(double));
    double *tCalculoX = (double*)malloc(sizeof(double));

    Matriz_t *matriz;

    processaLinhaComando(&pivoteamentoParcial, &arquivoSaida, argc, argv);

    while((ordemSistema = getchar()) != EOF) {
        // Verifica se é inicio de novo sistema e devolve pro stdin
        if(ordemSistema != linhaVazia) {
            ungetc(ordemSistema, stdin);
            matriz = lerMatriz();

            // Só procede se não houve erro durante a leitura
            if(matriz != NULL) {
                *tTriangularizacao = 0.0;
                *tCalculoX = 0.0;
                *tCalculoY = 0.0;

                if(!pivoteamentoParcial) {
                    execucaoGauss = eliminacaoGaussSemPivoteamento(matriz, tTriangularizacao);
                }
                else {
                    execucaoGauss = eliminacaoGaussComPivoteamento(matriz, tTriangularizacao);
                }

                // Só procede se houve sucesso na eliminação de Gauss
                if(!execucaoGauss) {
                    execucaoInversa = calculaInversa(matriz, tCalculoY, tCalculoX);

                    // Só procede se houve sucesso no cálculo da matriz inversa
                    if(!execucaoInversa) {
                        calculaNormas(matriz);
                        imprimeResultado(matriz, tTriangularizacao, tCalculoY, tCalculoX, arquivoSaida);
                    }
                }

                liberaMatriz(matriz);
            }
        }
    }

    free(tTriangularizacao);
    free(tCalculoY);
    free(tCalculoX);
}