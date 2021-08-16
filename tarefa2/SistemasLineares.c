#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
    real_t normaL2 = 0.0;

    calculaResiduo(SL->A, x, SL->b, res, SL->n);

    for(int i = 0; i < SL->n; i++) {
        normaL2 += pow(res[i], 2);
    }
    normaL2 = sqrt(normaL2);
    normaL2 = fabs(normaL2);

    return normaL2;
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
    double tempo = timestamp();

    // Cria copias locais
    unsigned int n = SL->n;
    real_t **matCoeficientes;
    matCoeficientes = copiaMatrizCoeficientes(SL->A, n);
    real_t *vetTermos;
    vetTermos = copiaVetorTermos(SL->b, n);
    if(matCoeficientes == NULL || vetTermos == NULL) {
        fprintf(stderr, "%s", "\nErro de alocação durante execução do método de Gauss. Cancelando método..."); 
        return -1;
    }

    // Execução Método de Gauss
    for(int i = 0; i < n; ++i) {
        unsigned int iPivo = encontraMax(matCoeficientes, n, i);
        if(i != iPivo) {
            trocaLinha(matCoeficientes, vetTermos, n, i, iPivo);
        }

        for(int k = i + 1; k < n; ++k) {
            real_t denominator = matCoeficientes[i][i];
            if(denominator == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do método de Gauss. Cancelando método...", i + 1);
                liberaCopiasLocais(matCoeficientes, vetTermos, n);
                return -1;
            }

            real_t m = matCoeficientes[k][i] / denominator;
            matCoeficientes[k][i] = 0.0;
            for(int j = i + 1; j < n; ++j) {
                matCoeficientes[k][j] -= matCoeficientes[i][j] * m;
            }

            vetTermos[k] -= vetTermos[i] * m;
        }
    }

    if(!possuiSolucao(matCoeficientes, vetTermos, n)) {
        fprintf(stderr, "%s", "\nSistema não possui solução. Cancelando análise do sistema..."); 
        return -2;
    }

    retrosubstituicao(matCoeficientes, x, vetTermos, n);

    // Libera espaço de cópias locais
    liberaCopiasLocais(matCoeficientes, vetTermos, n);

    tempo = timestamp() - tempo;
    *tTotal = tempo;

    return 0;
}

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{
    double tempo = timestamp();
    int iteracao = 0;

    // Cria copias locais
    unsigned int n = SL->n;
    real_t erro = SL->erro;

    // Inicializa vetor de aproximacoes
    real_t *vetInicial;
    vetInicial = malloc(n * sizeof(real_t));
    inicializaVetAproximacao(x, n);
    if(vetInicial == NULL) {
        fprintf(stderr, "%s", "\nErro de alocação durante execução do método de Jacobi. Cancelando método..."); 
        return -3;
    }

    if(!sistemaConverge(SL->A, n)) {
        fprintf(stderr, "%s", "\nSistema não converge durante execução do método de Jacobi.");
    }

    // Execução Método de Jacobi
    while(!menorQueErro(vetInicial, x, erro, n, iteracao) && iteracao < MAXIT) {
        // Passa resultado da iteracao antiga para vetInicial
        memcpy(vetInicial, x, n * sizeof(real_t));
        for(int i = 0; i < n; i++) {
            real_t soma = 0.0;
            real_t denominator = SL->A[i][i];
            if(denominator == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do método de Jacobi. Cancelando método...", i + 1);
                free(vetInicial);
                return -3;
            }

            for(int j = 0; j < n; j++) {
                if(i != j) {
                    soma += SL->A[i][j] * vetInicial[j];
                }
            }

            if(isinf(soma)) {
                fprintf(stderr, "\nOverflow na soma durante execução %d do método de Jacobi. Cancelando método...", i + 1);
                free(vetInicial);
                return -3;
            }

            x[i] = (SL->b[i] - soma) / denominator;
        }
        iteracao++;
    }

    // Libera espaço
    free(vetInicial);

    tempo = timestamp() - tempo;
    *tTotal = tempo;

    return iteracao;
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{
    double tempo = timestamp();
    int iteracao = 0;

    // Cria copias locais
    unsigned int n = SL->n;
    real_t erro = SL->erro;

    // Inicializa vetor de aproximacoes
    real_t *vetInicial;
    vetInicial = malloc(n * sizeof(real_t));
    inicializaVetAproximacao(x, n);

    if(vetInicial == NULL) {
        fprintf(stderr, "%s", "\nErro de alocação durante execução do método de Gauss-Seidel. Cancelando método..."); 
        return -3;
    }

    if(!sistemaConverge(SL->A, n)) {
        fprintf(stderr, "%s", "\nSistema não converge durante execução do método de Gauss-Seidel.");
    }

    // Execução Método de Gauss - Seidel
    while(!menorQueErro(vetInicial, x, erro, n, iteracao) && iteracao < MAXIT) {
        // Passa resultado da iteracao antiga para vetInicial
        memcpy(vetInicial, x, n * sizeof(real_t));
        for(int i = 0; i < n; i++) {
            real_t soma = 0.0;
            real_t denominator = SL->A[i][i];
            if(denominator == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do método de Gauss-Seidel. Cancelando método...", i + 1);
                free(vetInicial);
                return -3;
            }
            
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    soma += SL->A[i][j] * x[j];
                }
            }

            if(isinf(soma)) {
                fprintf(stderr, "\nOverflow na soma durante execução %d do método de Gauss-Seidel. Cancelando método...", i + 1);
                free(vetInicial);
                return -3;
            }

            x[i] = (SL->b[i] - soma) / denominator;
        }
        iteracao++;
    }

    // Libera espaço de cópias locais
    free(vetInicial);

    tempo = timestamp() - tempo;
    *tTotal = tempo;

    return iteracao;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{
    RegrasAplicadas_t *regras;
    regras = malloc(SL->n * SL->n * sizeof(RegrasAplicadas_t));
    int countRegras = 0;
    int iteracao = 0;

    // Cria copias locais
    unsigned int n = SL->n;
    real_t **matCoeficientes;
    matCoeficientes = copiaMatrizCoeficientes(SL->A, n);
    real_t *vetTermos;
    vetTermos = copiaVetorTermos(SL->b, n);

    if(matCoeficientes == NULL || vetTermos == NULL) {
        fprintf(stderr, "%s", "\nErro de alocação durante execução do refinamento. Cancelando método..."); 
        return -1;
    }

    // Execução Método de Gauss e guarda regras aplicadas
    for(int i = 0; i < n; ++i) {
        unsigned int iPivo = encontraMax(matCoeficientes, n, i);
        if(i != iPivo) {
            trocaLinha(matCoeficientes, vetTermos, n, i, iPivo);
            regras[countRegras].tipo = 0;
            regras[countRegras].linha1 = i;
            regras[countRegras].linha2 = iPivo;
            countRegras++;
        }

        for(int k = i + 1; k < n; ++k) {
            real_t denominator = matCoeficientes[i][i];
            if(denominator == 0.0) {
                fprintf(stderr, "\nDivisão por zero durante execução %d do refinamento. Cancelando método...", i + 1);
                liberaCopiasLocais(matCoeficientes, vetTermos, n);
                return -1;
            }

            real_t m = matCoeficientes[k][i] / denominator;
            matCoeficientes[k][i] = 0.0;
            
            for(int j = i + 1; j < n; ++j) {
                matCoeficientes[k][j] -= matCoeficientes[i][j] * m;
            }
            
            vetTermos[k] -= vetTermos[i] * m;
            regras[countRegras].tipo = 1;
            regras[countRegras].linha1 = k;
            regras[countRegras].linha2 = i;
            regras[countRegras].m = m;
            countRegras++;
        }
    }

    real_t *residuo = malloc(SL->n * sizeof(real_t));
    real_t *w = malloc(SL->n * sizeof(real_t));
    int condicao = 1;
    while(condicao && iteracao < MAXIT) {
        iteracao++;
        calculaResiduo(SL->A, x, SL->b, residuo, SL->n);
        aplicaRegrasResiduo(residuo, regras, countRegras);
        retrosubstituicao(matCoeficientes, w, residuo, n);
        for(int i = 0; i < n; i++) {
            x[i] = x[i] + w[i];
        }
        condicao = (normaL2Residuo(SL, x, residuo) >= SL->erro);
    }

    // Libera espaços
    free(w);
    free(residuo);
    liberaCopiasLocais(matCoeficientes, vetTermos, n);

    return iteracao;
}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n)
{
    // Aloca sistema linear
    SistLinear_t *sistemaLinear;
    if((sistemaLinear = malloc(sizeof(SistLinear_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar o sistema linear\n"); 
        return NULL;
    }

    // Aloca Matriz de Coeficientes
    if((sistemaLinear->A = malloc(n * sizeof(real_t*))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar a matriz de coeficientes\n"); 
        return NULL;
    }
    for(int i = 0; i < n; i++) {
        if((sistemaLinear->A[i] = malloc(n * sizeof(real_t))) == NULL) {
            fprintf(stderr, "%s", "Erro ao alocar linha da matriz de coeficientes\n"); 
            return NULL;
        }
    }

    // // Aloca Vetor de Termos Independentes
    if((sistemaLinear->b = malloc(n * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar o vetor de termos independentes\n"); 
        return NULL;
    }

    return sistemaLinear;
}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{
    unsigned int n = SL->n;

    // Desaloca Matriz de Coeficientes
    for(int i = 0; i < n; i++) {
        free(SL->A[i]);
    }
    free(SL->A);

    // Desaloca Vetor de Termos Independentes
    free(SL->b);
}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{
    long int ordemSistema;
    char linha[1024];
    SistLinear_t *sistemaLinear;

    if((fgets(linha, sizeof(linha), stdin) != NULL) && (!linhaVazia(linha))) {
        // Aloca sistema e guarda ordem
        ordemSistema = strtol(linha, (char **)NULL, 10);
        if((sistemaLinear = alocaSistLinear(ordemSistema)) == NULL) {
            return NULL;
        }
        sistemaLinear->n = ordemSistema;

        // Le erro e guarda
        if((fgets(linha, sizeof(linha), stdin) != NULL) && (!linhaVazia(linha))) {
            sistemaLinear->erro = strtof(linha, NULL);
        }
        else {
            fprintf(stderr, "%s", "Erro ao ler valor do critério de parada\n"); 
            return NULL;
        }

        // Le matriz de coeficientes
        for(int i = 0; i < ordemSistema; i++) {
            // Le linhas da matriz e processa
            if((fgets(linha, sizeof(linha), stdin) != NULL) && (!linhaVazia(linha))) {
                char* token;
                char* rest = linha;
            
                for(int j = 0; j < ordemSistema; j++) {
                    token = strtok_r(rest, " ", &rest);
                    // Retorna se o valor de alguma coluna esta faltando
                    if(token == NULL) {
                        fprintf(stderr, "%s", "Erro ao ler valor da coluna da matriz de coeficientes\n"); 
                        return NULL;
                    }

                    sistemaLinear->A[i][j] = (real_t)strtof(token, NULL);
                }
            }
            else {
                fprintf(stderr, "%s", "Erro ao ler linha da matriz de coeficientes\n"); 
                return NULL;
            }
        }

        // Le vetor de termos independentes
        if((fgets(linha, sizeof(linha), stdin) != NULL) && (!linhaVazia(linha))) {
            char* token;
            char* rest = linha;
        
            for(int i = 0; i < ordemSistema; i++) {
                token = strtok_r(rest, " ", &rest);
                // Retorna se o valor de alguma coluna esta faltando
                if(token == NULL) {
                    fprintf(stderr, "%s", "Erro ao ler valor do vetor de termos independentes\n"); 
                    return NULL;
                }

                sistemaLinear->b[i] = strtof(token, NULL);
            }
        }
        else {
            fprintf(stderr, "%s", "Erro ao ler vetor de termos independentes\n"); 
            return NULL;
        }
    }
    else {
        fprintf(stderr, "%s", "Erro ao ler ordem do sistema linear\n"); 
        return NULL;
    }

    return sistemaLinear;
}


// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL)
{
    unsigned int n = SL->n;

    // Imprime ordem do sistema linear
    printf("%d\n", n);

    // Imprime erro do sistema linear
    printf("%f\n", SL->erro);

    // Imprime matriz de coeficientes
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%1.8e ", SL->A[i][j]);
        }
        printf("\n");
    }

    // Imprime vetor de termos independentes
    for(int i = 0; i < n; i++) {
        printf("%1.8e ", SL->b[i]);
    }
    printf("\n");
}

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
    for(int i = 0; i < n; i++) {
        printf("%1.8e ", v[i]);
    }
    printf("\n");
}

