// ###################################
// Gabriel Marckzuk Thá - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "matrizes.h"
#include "sistemasLineares.h"

/*!
    \brief Método de Gauss para triangularização da
            matriz sem pivoteamento parcial.

    \param matriz Estrutura contendo a matriz a ser triangularizada.
    \param tTotal Ponteiro para a variável que guardará o tempo
                    da triangularização, em milissegundos.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
            Dois caso haja erro de alocação da matriz auxiliar.
*/
int eliminacaoGaussSemPivoteamento(Matriz_t *matriz, double *tTotal) {
    unsigned int n = matriz->n;
    real_t **matrizCopia;
    matrizCopia = copiaMatriz(matriz->A, n);

    if(matrizCopia == NULL) {
        fprintf(stderr, "%s", "\nErro ao copiar matriz durante execução do método de Gauss. Cancelando método..."); 
        return 2;
    }

    double tempo = timestamp();

    // Trata a primeira linha da matriz L
    matriz->L[0][0] = 1.0;

    // Execução Método de Gauss
    for(int i = 0; i < n; ++i) {
        for(int k = i + 1; k < n; ++k) {
            real_t denominador = matrizCopia[i][i];
            if(denominador == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do método de Gauss. Cancelando método...", i + 1);
                liberaCopia(matrizCopia, n);
                return 1;
            }

            real_t m = matrizCopia[k][i] / denominador;
            matrizCopia[k][i] = 0.0;

            matriz->L[k][i] = m;
            matriz->L[k][i + 1] = 1.0;

            for(int j = i + 1; j < n; ++j) {
                matrizCopia[k][j] -= matrizCopia[i][j] * m;
            }
        }
    }

    tempo = timestamp() - tempo;
    *tTotal = tempo;

    // Copia matriz escalonada para U
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < (n - i); j++) {
            matriz->U[i][j] = matrizCopia[i][j + i];
        }
    }

    liberaCopia(matrizCopia, n);
    return 0;
}

/*!
    \brief Método de Gauss para triangularização da
            matriz com pivoteamento parcial.

    \param matriz Estrutura contendo a matriz a ser triangularizada.
    \param tTotal Ponteiro para a variável que guardará o tempo
                    da traingularização em milissegundos.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
            Dois caso haja erro de alocação da matriz auxiliar.
*/
int eliminacaoGaussComPivoteamento(Matriz_t *matriz, double *tTotal) {
    unsigned int n = matriz->n;
    real_t **matrizCopia;
    matrizCopia = copiaMatriz(matriz->A, n);

    if(matrizCopia == NULL) {
        fprintf(stderr, "%s", "\nErro ao copiar matriz durante execução do método de Gauss. Cancelando método..."); 
        return 2;
    }

    double tempo = timestamp();

    // Trata a primeira linha da matriz L
    matriz->L[0][0] = 1.0;

    // Execução Método de Gauss
    for(int i = 0; i < n; ++i) {
        unsigned int iPivo = encontraMax(matrizCopia, n, i);
        if(i != iPivo) {
            trocaLinha(matrizCopia, matriz->L, matriz->ordemLinhas, n, i, iPivo);
        }

        for(int k = i + 1; k < n; ++k) {
            real_t denominador = matrizCopia[i][i];
            real_t m = matrizCopia[k][i] / denominador;
            if(denominador == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do método de Gauss. Cancelando método...", i + 1);
                liberaCopia(matrizCopia, n);
                return 1;
            }
            matrizCopia[k][i] = 0.0;

            matriz->L[k][i] = m;
            matriz->L[k][i + 1] = 1.0;

            for(int j = i + 1; j < n; ++j) {
                matrizCopia[k][j] -= matrizCopia[i][j] * m;
            }
        }
    }

    tempo = timestamp() - tempo;
    *tTotal = tempo;

    // Copia matriz escalonada para U
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < (n - i); j++) {
            matriz->U[i][j] = matrizCopia[i][j + i];
        }
    }

    liberaCopia(matrizCopia, n);
    return 0;
}

/*!
    \brief Cálculo do vetor y da fatoração LU.

    \param L Matriz L da fatoração LU, definida em "matrizes.h".
    \param y Ponteiro para o vetor y que será calculado.
    \param n Ordem da matriz L.
    \param idx Índice que contém o valor 1 na coluna da matriz
            identidade.
*/
void calculaY(real_t **L, real_t *y, unsigned int n, int idx) {
    real_t *b = (real_t*)malloc(n * sizeof(real_t));
   
    for(int i = 0; i < n; i++) {
        b[i] = 0.0;
    }
    b[idx] = 1.0;

    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;
        for(int j = 0; j < i; j++) {
            soma += L[i][j] * y[j];
        }
        y[i] = b[i] - soma;
    }

    free(b);
}

/*!
    \brief Cálculo do vetor x da fatoração LU.

    \param U Matriz U da fatoração LU, definida em "matrizes.h".
    \param x Ponteiro para o vetor x que será calculado.
    \param n Ordem da matriz U.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
*/
int calculaX(real_t **U, real_t *x, real_t *y, unsigned int n) {
    for(int i = n - 1; i >= 0; i--) {
        real_t soma = 0.0;
        for(int j = i + 1; j < n; j++) {
            soma += U[i][j - i] * x[j];
        }

        real_t denominador = U[i][0];
        if(denominador == 0.0) {
            fprintf(stderr, "\nDivisão por zero durante cálculo do vetor X. Cancelando método...");
            return 1;
        }
        x[i] = (y[i] - soma) / denominador;
    }

    return 0;
}

/*!
    \brief Método de para o cálcula da inversa da matriz de entrada
            a partir das matrizes L e U obtidas durante a fatoração.

    \param matriz Estrutura contendo os dados para o cálculo da
            inversa da matriz de entrada.
    \param tTotalY Ponteiro para a variável que guardará a média de 
            tempo do cálculo de todos os vetores y da fatoração LU,
            em milissegundos.
    \param tTotalX Ponteiro para a variável que guardará a média de 
            tempo do cálculo de todos os vetores x da fatoração LU,
            em milissegundos.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
*/
int calculaInversa(Matriz_t *matriz, double *tTotalY, double *tTotalX) {
    unsigned int n = matriz->n;
    double tempoY = 0.0, tempoX = 0.0, tempo;
    int execucaoX;

    real_t *y = (real_t*)malloc(n * sizeof(real_t));
    real_t *x = (real_t*)malloc(n * sizeof(real_t));
    
    for(int i = 0; i < n; i++) {
        tempo = timestamp();
        calculaY(matriz->L, y, n, matriz->ordemLinhas[i]);
        tempoY += timestamp() - tempo;

        tempo = timestamp();
        execucaoX = calculaX(matriz->U, x, y, n);
        tempoX += timestamp() - tempo;
        if(execucaoX != 0) {
            return 1;
        }

        for(int j = 0; j < n; j++) {
            matriz->AI[j][i] = x[j];
        }
    }

    free(y);
    free(x);

    *tTotalY = tempoY / n;
    *tTotalX = tempoX / n;

    return 0;
}

/*!
    \brief Cálculo do resíduo em relação a uma determinada coluna
            das matrizes.

    \param matriz Estrutura contendo os dados para o cálculo do resíduo.
    \param b Coluna da matriz identidade a ser utilizada no cálculo.
    \param res Ponteiro para a variável que guardará o valor dos resíduos.
    \param idx Índice da coluna da matriz que será utilizada para calcular
            o resíduo.
*/
void calculaResiduo(Matriz_t *matriz, real_t *b, real_t *res, int idx) {
    unsigned int n = matriz->n;
    real_t *x = (real_t*)malloc(n * sizeof(real_t));
    for(int i = 0; i < n; i++) {
        x[i] = matriz->AI[i][idx];
    }
    
    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;

        for(int j = 0; j < n; j++) {
            soma += matriz->A[i][j] * x[j];
        }

        res[i] = b[i] - soma;
    }

    free(x);
}

/*!
    \brief Cálculo da norma L2 em relação a uma determinada
            coluna das matrizes.

    \param matriz Estrutura definida em "matrizes.h" contendo os dados para o cálculo das normas.
    \param idx Índice da coluna das matrizes inversa e identidade que serão
                utilizada para o cálculo da norma.

    \return Norma L2 da idx-ésima coluna.
*/
real_t normaL2Residuo(Matriz_t *matriz, int idx) {
    unsigned int n = matriz->n;
    real_t *res = (real_t*)malloc(n * sizeof(real_t));
    real_t *b = (real_t*)malloc(n * sizeof(real_t));

    for(int i = 0; i < n; i++) {
        b[i] = 0.0;
    }
    b[idx] = 1.0;

    real_t normaL2 = 0.0;

    calculaResiduo(matriz, b, res, idx);

    for(int i = 0; i < n; i++) {
        normaL2 += pow(res[i], 2);
    }
    normaL2 = sqrt(normaL2);
    normaL2 = fabs(normaL2);

    free(b);
    free(res);
    return normaL2;
}

/*!
    \brief Cálculo da norma L2 de cada coluna da matriz inversa
        em relação a coluna correspondente da matriz identidade.

    \param matriz Estrutura contendo os dados para o cálculo das normas.
*/
void calculaNormas(Matriz_t *matriz) {
    unsigned int n = matriz->n;

    for(int i = 0; i < n; i++) {
        matriz->normas[i] = normaL2Residuo(matriz, matriz->ordemLinhas[i]);
    }
}