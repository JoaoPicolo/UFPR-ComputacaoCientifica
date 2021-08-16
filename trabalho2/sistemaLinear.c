// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "utils.h"
#include "sistemaLinear.h"

/*!
    \brief Alocação dinâmica da estrutura Sistema_t definida em "sistemaLinear.h".

    \param n Ordem do sistema.

    \return Ponteiro para o sistema. NULL se houve erro de alocação.
*/
Sistema_t* alocaSistema(unsigned int n) {
    Sistema_t *sistema;

    // Aloca sistema
    if((sistema = (Sistema_t*)malloc(sizeof(Sistema_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar sistema linear.\n"); 
        return NULL;
    }

    // Aloca matriz de coeficientes
    if((sistema->U = (real_t*)malloc(n * n * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar matriz de coeficientes.\n"); 
        return NULL;
    }

    // Aloca matriz auxiliar L
    // Aloca somente espaços necessários, L tem tamanho da soma
    // dos n primeiros naturais
    int tamanho = ceil((n*n + n) / 2);
    if((sistema->L = (real_t*)malloc(tamanho * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar matriz U.\n"); 
        return NULL;
    }

    // Aloca vetor de termos independentes
    if((sistema->b = (real_t*)malloc(n * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar vetor de termos independentes.\n"); 
        return NULL;
    }

    return sistema;
}

/*!
    \brief Libera sistema na memória.

    \param sistema Ponteiro para o sistema a ser liberado.
*/
void liberaSistema(Sistema_t *sistema) {
    free(sistema->U);
    free(sistema->L);
    free(sistema->b);
    free(sistema);
}

/*!
    \brief Preenche o sistema linear utilizando a lógica de preenchimento para
            o método de interpolação polinomial.

    \param sistema Ponteiro para o sistema que irá armazenar os valores.
    \param valoresTabelados Ponteiro para um vetor contendo os valores X tabelados.
    \param funcoesTabeladas Ponteiro para um vetor contendo os valores de f(X) tabelados.
*/
void preencheSistemaInterpolacao(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas) {
    unsigned int n = sistema->n;
    
    memcpy(sistema->b, funcoesTabeladas, n * sizeof(real_t));
    
    for(int i = 0; i < n; i++) {
        real_t valor = valoresTabelados[i];

        for(int j = 0; j < n; j++) {
            sistema->U[i*n + j] = pow(valor, j);
        }
    }
}

/*!
    \brief Preenche o sistema linear utilizando a lógica de preenchimento para
            o método de ajuste de curvas.

    \param sistema Ponteiro para o sistema que irá armazenar os valores.
    \param valoresTabelados Ponteiro para um vetor contendo os valores X tabelados.
    \param funcoesTabeladas Ponteiro para um vetor contendo os valores de f(X) tabelados.
*/
void preencheSistemaAjuste(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas) {
    unsigned int n = sistema->n;
    
    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;
        for(int j = 0; j < n; j ++) {
            soma += funcoesTabeladas[j] * pow(valoresTabelados[j], i);   
        }
        sistema->b[i] = soma;
    }

    // Monta a primeira linha
    for(int i = 0; i < n; i++) {
        // Cuida do somatório
        real_t soma = 0.0;
        for(int j = 0; j < n; j++) {
            soma += pow(valoresTabelados[j], i);
        }

        sistema->U[i] = soma;
    }

    // Monta linhas seguintes
    // Só precisa calcular o último valor, o resto repete da anterior
    for(int i = 1; i < n; i++) {
        // Preenche até n - 1 colunas
        // Preenche com linha anterior e coluna seguinte
        for(int j = 0; j < n - 1; j++) {
            sistema->U[i*n + j] = sistema->U[(i-1)*n + (j+1)];
        }

        // Calcula a última coluna
        real_t soma = 0.0;
        for(int j = 0; j < n; j++) {
            soma += pow(valoresTabelados[j], n-1) * pow(valoresTabelados[j], i);
        }
        sistema->U[i*n + (n-1)] = soma;
    }
}

/*!
    \brief Preenche o sistema linear utilizando a lógica de preenchimento para
            o método de ajuste de curvas utilizando otimização no código.

    \param sistema Ponteiro para o sistema que irá armazenar os valores.
    \param valoresTabelados Ponteiro para um vetor contendo os valores X tabelados.
    \param funcoesTabeladas Ponteiro para um vetor contendo os valores de f(X) tabelados.
*/
void preencheSistemaAjusteOtimizada(Sistema_t *sistema, real_t *valoresTabelados, real_t *funcoesTabeladas) {
    unsigned int n = sistema->n;

    real_t *lookupPow = (real_t*) malloc(n * n * sizeof(real_t));

    // Preenche primeira coluna da primeira linha
    // assim como o primeiro elemento do vetor de termos independentes
    real_t soma_b = 0.0;
    real_t soma_U = 0.0;
    real_t potencia = 1.0;
    for(int j = 0; j < n; j ++) {
        lookupPow[j] = potencia;
        soma_b += funcoesTabeladas[j];
        soma_U += potencia;
    }
    sistema->b[0] = soma_b;
    sistema->U[0] = soma_U;
    
    // Preenche colunas restantes da primeira linha
    // assim como os elementos seguintes do vetor de termos independentes
    for(int i = 1; i < n; i++) {
        soma_b = 0.0;
        soma_U = 0.0;
        for(int j = 0; j < n; j ++) {
            potencia = lookupPow[(i-1)*n + j] * valoresTabelados[j];
            lookupPow[(i*n) + j] = potencia;
            soma_b += funcoesTabeladas[j] * potencia;
            soma_U += potencia;
        }
        sistema->b[i] = soma_b;
        sistema->U[i] = soma_U;
    }

    // Monta linhas seguintes
    // Só precisa calcular o último valor, o resto repete da anterior
    int ultimaCol = n - 1;
    for(int i = 1; i < n; i++) {
        // Preenche até n - 1 colunas
        // Preenche com linha anterior e coluna seguinte
        int linhaAnt = (i - 1) * n;
        int linhaAtual = i * n;
        for(int j = 0; j < n - 1; j++) {
            sistema->U[linhaAtual + j] = sistema->U[linhaAnt + (j+1)];
        }

        // Calcula a última coluna
        real_t soma = 0.0;
        for(int j = 0; j < n; j++) {
            soma += lookupPow[(ultimaCol*n) + j] * lookupPow[(i*n) + j];
        }
        sistema->U[linhaAtual + ultimaCol] = soma;
    }

    free(lookupPow);
}


/*!
    \brief Encontra o valor máximo da coluna de uma matriz.

    \param matriz Matriz U da fatoração LU a ser analisada.
    \param n Ordem da matriz a ser analisada.
    \param col Índice da coluna a ser analisada.

    \return Índice do maior valor da coluna analisada.
*/
unsigned int encontraMax(real_t *matriz, unsigned int n, unsigned int col) {
    unsigned int indexMaior = col;

    for(unsigned int i = col + 1; i < n; i++) {
        if(fabs(matriz[i*n + col]) > fabs(matriz[indexMaior*n + col])) {
            indexMaior = i;
        }
    }

    return indexMaior;
}

/*!
    \brief Troca as linhas de uma matriz U, da matriz L da fatoração LU
            e do vetor b de termos independentes.

    \param sistema Sistema contendo as informações das matrizes L, U
            e do vetor de termos independentes b que deverão ter suas linhas trocadas.
    \param n Ordem da matriz a ter as linhas trocadas.
    \param i Índice da primeira linha a ser trocada.
    \param iPivo Índice da segunda linha a ser trocada.
*/
void trocaLinha(Sistema_t *sistema, unsigned int n, unsigned int i, unsigned int iPivo) {
    
    // Troca linhas na matriz U
    real_t aux;
    for(int j = 0; j < n; j++) {
        aux = sistema->U[i*n + j];
        sistema->U[i*n + j] = sistema->U[iPivo*n + j];
        sistema->U[iPivo*n + j] = aux;
    }

    // Troca termos independentes
    aux = sistema->b[i];
    sistema->b[i] = sistema->b[iPivo];
    sistema->b[iPivo] = aux;

    // Troca linhas das matriz L se já possui valores salvos
    if(i > 0) {
        int baseI = ceil((i*i + i) / 2);
        int baseIPivo = ceil((iPivo*iPivo + iPivo) / 2);
        // Como i corresponde a coluna atual da iteração sabemos que
        // L foi preenchida até i - 1, por isso trocamos todas as linhas
        // até este índice
        for(int j = 0; j <= i - 1; j++) {
            aux = sistema->L[baseI + j];
            sistema->L[baseI + j] = sistema->L[baseIPivo + j];
            sistema->L[baseIPivo + j] = aux;
        }
    }
}

/*!
    \brief Triangularização do sistema linear através do método de Gauss,
            para posterior resolução com fatoração LU.

    \param sistema Ponteiro para o sistema que deverá ser triangularizado.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
*/
int triangularizacaoSistema(Sistema_t *sistema) {
    unsigned int n = sistema->n;

    // Preenche a primeira diagonal de L
    sistema->L[0] = 1.0;

    // Execução Método de Gauss
    for(int i = 0; i < n; ++i) {
        unsigned int iPivo = encontraMax(sistema->U, n, i);
        if(i != iPivo) {
            trocaLinha(sistema, n, i, iPivo);
        }

        for(int k = i + 1; k < n; ++k) {
            real_t denominador = sistema->U[i*n + i];

            real_t m = sistema->U[k*n + i] / denominador;
            sistema->U[k*n + i] = 0.0;

            int tamanho = ceil((k*k + k) / 2);
            sistema->L[tamanho + i] = m;
            sistema->L[tamanho + i + 1] = 1.0;

            for(int j = i + 1; j < n; ++j) {
                sistema->U[k*n + j] -= sistema->U[i*n + j] * m;
            }
        }
    }

    return 0;
}

/*!
    \brief Triangularização do sistema linear através do método de Gauss,
            para posterior resolução com fatoração LU. Difere da anterior por
            estar utilizando técnicas de otimização.

    \param sistema Ponteiro para o sistema que deverá ser triangularizado.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
*/
int triangularizacaoSistemaOtimizada(Sistema_t *sistema) {
    unsigned int n = sistema->n;

    // Preenche a primeira diagonal de L
    sistema->L[0] = 1.0;
    int stride = 8;
    int limite = n - (n % stride);

    // Execução Método de Gauss
    for(int i = 0; i < n; ++i) {
        int linhaAnt = i * n;
        unsigned int iPivo = encontraMax(sistema->U, n, i);
        if(i != iPivo) {
            trocaLinha(sistema, n, i, iPivo);
        }

        for(int k = i + 1; k < n; ++k) {
            int linhaAtual = k * n;
            real_t denominador = sistema->U[linhaAnt + i];

            real_t m = sistema->U[linhaAtual + i] / denominador;
            sistema->U[linhaAtual + i] = 0.0;

            int tamanho = ceil((k*k + k) / 2);
            sistema->L[tamanho + i] = m;
            sistema->L[tamanho + i + 1] = 1.0;

            int j = 0;
            for(j = i + 1; j < limite; j += stride) {
                sistema->U[linhaAtual + j] -= sistema->U[linhaAnt + j] * m;
                sistema->U[linhaAtual + j + 1] -= sistema->U[linhaAnt + j + 1] * m;
                sistema->U[linhaAtual + j + 2] -= sistema->U[linhaAnt + j + 2] * m;
                sistema->U[linhaAtual + j + 3] -= sistema->U[linhaAnt + j + 3] * m;
                sistema->U[linhaAtual + j + 4] -= sistema->U[linhaAnt + j + 4] * m;
                sistema->U[linhaAtual + j + 5] -= sistema->U[linhaAnt + j + 5] * m;
                sistema->U[linhaAtual + j + 6] -= sistema->U[linhaAnt + j + 6] * m;
                sistema->U[linhaAtual + j + 7] -= sistema->U[linhaAnt + j + 7] * m;
            }

            for(; j < n; ++j) {
                sistema->U[linhaAtual + j] -= sistema->U[linhaAnt + j] * m;
            }
        }
    }

    return 0;
}

/*!
    \brief Cálculo do vetor y da fatoração LU, respeitando o seguinte: Ly = b.

    \param L Matriz L da fatoração LU, definida em "sistemaLinear.h".
    \param y Ponteiro para o vetor y que será calculado.
    \param b Ponteiro para o vetor b que será utilizado no cálculo.
    \param n Ordem da matriz L.
*/
void calculaY(real_t *L, real_t *y, real_t *b, unsigned int n) {
    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;
        int tamanho = ceil((i*i + i) / 2);
        for(int j = 0; j < i; j++) {
            soma += L[tamanho + j] * y[j];
        }
        y[i] = b[i] - soma;
    }
}

/*!
    \brief Cálculo do vetor x da fatoração LU, respeitando o seguinte: Ux = y.

    \param U Matriz U da fatoração LU, definida em "sistemaLinear.h".
    \param x Ponteiro para o vetor x que será calculado.
    \param y Ponteiro para o vetor y que será utilizado no cálculo.
    \param n Ordem da matriz U.

    \return Zero em caso de sucesso. Um em caso de divisão por zero.
*/
int calculaX(real_t *U, real_t *x, real_t *y, unsigned int n) {
    for(int i = n - 1; i >= 0; i--) {
        real_t soma = 0.0;
        for(int j = i + 1; j < n; j++) {
            soma += U[i*n + j] * x[j];
        }

        real_t denominador = U[i*n + i];
        if(denominador == 0.0) {
            fprintf(stderr, "\nDivisão por zero durante execução %d do cálculo de X. Cancelando método...", i + 1);
            return 1;
        }
        x[i] = (y[i] - soma) / denominador;
    }

    return 0;
}

/*!
    \brief Resolução do sistema linear através do método de fatoração LU.

    \param sistema Ponteiro para o sistema que deverá ser resolvido.
    \param resultado Ponteiro para vetor contendo o resultado x provido
                    pela fatoração LU.

    \return Zero em caso de sucesso. Um em caso de divisão por zero no cálculo de X.
*/
int resolveSistema(Sistema_t *sistema, real_t *resultado) {
    unsigned int n = sistema->n;
    real_t *y = (real_t*)malloc(n * sizeof(real_t));

    calculaY(sistema->L, y, sistema->b, n);
    int saida = calculaX(sistema->U, resultado, y, n);

    free(y);

    return saida;
}