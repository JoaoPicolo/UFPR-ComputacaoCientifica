// ###################################
// Gabriel Marckzuk Thá - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "utils.h"
#include "matrizes.h"

/*!
    \brief Alocação dinâmica da estrutura Matriz_t definida em "matrizes.h".

    \param n Ordem da matriz a ser alocada.

    \return Ponteiro para a matriz. NULL se houve erro de alocação.
*/
Matriz_t* alocaMatriz(unsigned int n) {
    Matriz_t* matriz;

    // Aloca estrutura
    if((matriz = (Matriz_t*)malloc(sizeof(Matriz_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar estrutura de armazenamento.\n"); 
        return NULL;
    }

    // Aloca espaço para matriz
    if((matriz->A = (real_t**)malloc(n * sizeof(real_t*))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de armazenamento da matriz.\n"); 
        return NULL;
    }
    for(int i = 0; i < n; i++) {
        if((matriz->A[i] = (real_t*)malloc(n * sizeof(real_t))) == NULL) {
            fprintf(stderr, "%s", "Erro ao alocar linha da matriz.\n"); 
            return NULL;
        }
    }

    // Aloca espaço para matriz inversa
    if((matriz->AI = (real_t**)malloc(n * sizeof(real_t*))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de armazenamento da matriz inversa.\n"); 
        return NULL;
    }
    for(int i = 0; i < n; i++) {
        if((matriz->AI[i] = (real_t*)malloc(n * sizeof(real_t))) == NULL) {
            fprintf(stderr, "%s", "Erro ao alocar linha da matriz inversa.\n"); 
            return NULL;
        }
    }

    // Aloca espaço para matriz L
    if((matriz->L = (real_t**)malloc(n * sizeof(real_t*))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de armazenamento da matriz L.\n"); 
        return NULL;
    }
    for(int i = 0; i < n; i++) {
        if((matriz->L[i] = (real_t*)malloc((i + 1) * sizeof(real_t))) == NULL) {
            fprintf(stderr, "%s", "Erro ao alocar linha da matriz L.\n"); 
            return NULL;
        }
    }

    // Aloca espaço para matriz U
    if((matriz->U = (real_t**)malloc(n * sizeof(real_t*))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de armazenamento da matriz U.\n"); 
        return NULL;
    }
    for(int i = 0; i < n; i++) {
        if((matriz->U[i] = (real_t*)malloc((n - i) * sizeof(real_t))) == NULL) {
            fprintf(stderr, "%s", "Erro ao alocar linha da matriz U.\n"); 
            return NULL;
        }
    }

    // Aloca espaço para vetor de ordem das linhas
    if((matriz->ordemLinhas = (int*)malloc(n * sizeof(int))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de estrutura auxiliar.\n"); 
        return NULL;
    }
    
    // Aloca espaço para vetor de normas
    if((matriz->normas = (real_t*)malloc(n * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço de estrutura auxiliar.\n"); 
        return NULL;
    }

    return matriz;
}

/*!
    \brief Libera matriz na memória.

    \param matriz Ponteiro para a matriz a ser liberada.
*/
void liberaMatriz(Matriz_t *matriz) {
    unsigned int n = matriz->n;

    for(int i = 0; i < n; i++) {
        free(matriz->A[i]);
        free(matriz->AI[i]);
        free(matriz->L[i]);
        free(matriz->U[i]);
    }

    free(matriz->ordemLinhas);
    free(matriz->normas);
    free(matriz);
}

/*!
    \brief Cria e retorna a cópia de uma matriz.

    \param matriz Matriz a ser copiada.
    \param n Ordem da matriz a ser copiada.

    \return Ponteiro para a copia da matriz. NULL se houve erro de alocação.
*/
real_t** copiaMatriz(real_t **matriz, unsigned int n) {
    real_t **copia;

    copia = (real_t**)malloc(n * sizeof(real_t*));
    memcpy(copia, matriz, n * sizeof(real_t*));

    for(int i = 0; i < n; i++) {
        copia[i] = (real_t*)malloc(n * sizeof(real_t));
        memcpy(copia[i], matriz[i], n * sizeof(real_t));
    }

    return copia;
}

/*!
    \brief Libera cópia da matriz na memória.

    \param matriz Ponteiro para a matriz a ser liberada.
*/
void liberaCopia(real_t **matriz, unsigned int n) {
    for(int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

/*!
    \brief Leitura da matriz a partir de entrada padrão (stdin).

    \return Ponteiro para a estrutura Matriz_t definida em "matrizes.h".
            NULL se houve erro de alocação.
*/
Matriz_t* lerMatriz() {
    unsigned int n;
    
    Matriz_t *matriz;

    scanf("%d", &n);

    matriz = alocaMatriz(n);

    if(matriz != NULL) {
        matriz->n = n;

        for (int i = 0; i < n; i++) {
            matriz->ordemLinhas[i] = i; // Guarda ordem inicial das linhas
            for (int j = 0; j < n; j++)
                scanf("%f", &matriz->A[i][j]);
        }   
    }

    return matriz;
}

/*!
    \brief Encontra o valor máximo da coluna de uma matriz.

    \param matriz Matriz a ser analisada.
    \param n Ordem da matriz a ser analisada.
    \param col Índice da coluna a ser analisada.

    \return Índice do maior valor da coluna analisada.
*/
unsigned int encontraMax(real_t **matriz, unsigned int n, unsigned int col) {
    unsigned int indexMaior = col;

    for(unsigned int i = col + 1; i < n; i++) {
        if(fabs(matriz[i][col]) > fabs(matriz[indexMaior][col])) {
            indexMaior = i;
        }
    }

    return indexMaior;
}

/*!
    \brief Troca as linhas de uma matriz, da matriz L da fatoração LU
            e altera a ordem das linhas da matriz identidade através
            da estrutura auxiliar definida em "matrizes.h".

    \param matriz Matriz a ter as linhas trocadas.
    \param L Matriz L da fatoração LU, a ter seus elementos trocados.
    \param ordemLinhas Estrutura auxiliar que define as linhas da matriz identidade.
    \param n Ordem da matriz a ter as linhas trocadas.
    \param i Índice da primeira linha a ser trocada.
    \param iPivo Índice da segunda linha a ser trocada.
*/
void trocaLinha(real_t **matriz, real_t **L, int *ordemLinhas, 
    unsigned int n, unsigned int i, unsigned int iPivo) {
    
    // Troca linhas na matriz
    real_t aux;
    for(int j = 0; j < n; j++) {
        aux = matriz[i][j];
        matriz[i][j] = matriz[iPivo][j];
        matriz[iPivo][j] = aux;
    }

    // Troca linhas das matriz L se já possui valores salvos
    if(i > 0) {
        real_t temp;
        // Como i corresponde a coluna atual da iteração sabemos que
        // L foi preenchida até i - 1, por isso trocamos todas as linhas
        // até este índice
        for(int j = 0; j <= i - 1; j++) {
            temp = L[i][j];
            L[i][j] = L[iPivo][j];
            L[iPivo][j] = temp;
        }
    }

    // Troca ordem das linhas na estrutura que define a matriz identidade
    int index1 = -1;
    int index2 = -1;
    int count = 0;
    while(index1 == -1 || index2 == -1) {
        if(ordemLinhas[count] == i) {
            index1 = count;
        }
        else if(ordemLinhas[count] == iPivo) {
            index2 = count;
        }
        count++;
    }

    aux = ordemLinhas[index1];
    ordemLinhas[index1] = ordemLinhas[index2];
    ordemLinhas[index2] = aux;
}
