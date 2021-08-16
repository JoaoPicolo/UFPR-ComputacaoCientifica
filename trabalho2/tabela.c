// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "utils.h"
#include "tabela.h"
#include "sistemaLinear.h"
#include "comunicacaoIO.h"

/*!
    \brief Alocação dinâmica da estrutura Tabela_t definida em "tabela.h".

    \param n Número de valores tabelados.
    \param m Número de funções tabeladas.

    \return Ponteiro para a tabela. NULL se houve erro de alocação.
*/
Tabela_t* alocaTabela(unsigned int n, unsigned int m) {
    Tabela_t *tabela;

    // Aloca tabela
    if((tabela = (Tabela_t*)malloc(sizeof(Tabela_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar estrutura da tabela.\n"); 
        return NULL;
    }

    // Aloca vetor de valores tabelados
    if((tabela->x = (real_t*)malloc(n * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço para vetor de valores tabelados.\n"); 
        return NULL;
    }

    // Aloca vetor de funções tabeladas
    if((tabela->f = (real_t*)malloc(n * m * sizeof(real_t))) == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço para vetor de funções tabeladas.\n"); 
        return NULL;
    }

    return tabela;
}

/*!
    \brief Libera tabela na memória.

    \param matriz Ponteiro para a tabela a ser liberada.
*/
void liberaTabela(Tabela_t *tabela) {
    free(tabela->x);
    free(tabela->f);
    free(tabela);
}


/*!
    \brief Leitura dos valores tabelados a partir de entrada padrão (stdin).

    \return Ponteiro para a estrutura Tabela_t definida em "tabela.h".
            NULL se houve erro de alocação.
*/
Tabela_t* lerTabela() {
    unsigned int n, m;
    int result = 0;
    if(result) result = 0;
    
    Tabela_t *tabela;

    result = scanf("%d", &n);
    result = scanf("%d", &m);

    tabela = alocaTabela(n, m);

    if(tabela != NULL) {
        tabela->n = n;
        tabela->m = m;

        // Lê valores tabelados
        for (int i = 0; i < n; i++)
            result = scanf("%lf", &tabela->x[i]);

        // Lê valores das funções tabeladas
        for (int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                result = scanf("%lf", &tabela->f[i*n + j]);
    }

    return tabela;
}

/*!
    \brief Responsável pelo processamento da Tabela de entrada de forma a realizar os
        cálculos necessários para a Interpolação Polinomial e o Ajuste de Curvas.
        Após o preenchimento do sistema de cada etapa a partir da Tabela lida, executa a
        triangularização, calcula o resultado e imprime na saída padrão (stdout).

    \param tabela Ponteiro para a estrutura Tabela_t definida em "tabela.h".

    \return Zero em caso de sucesso. Um em caso de erro de alocação.
*/
int executaFuncoesPolinomiais(Tabela_t *tabela) {
    unsigned int n = tabela->n;
    unsigned int m = tabela->m;
    real_t *funcoes = (real_t*)malloc(n * sizeof(real_t));
    real_t *resultado = (real_t*)malloc(n * sizeof(real_t));
    if(funcoes == NULL || resultado == NULL) {
        fprintf(stderr, "%s", "Erro ao alocar espaço para vetores auxiliares.\n"); 
        return 1;
    }

    Sistema_t *sistema = alocaSistema(n);
    if(sistema == NULL) {
        return 1;
    }
    sistema->n = n;

    double tempo  = 0.0;
    double tempo_p = 0.0;
    double tempo_t  = 0.0;
    int i;
    int triangularizacao, solucao;
    for(i = 0; i < m; i++) {
        memcpy(funcoes, tabela->f + (n * i), n*sizeof(real_t));

        // Resolve utilizando interpolação polinomial
        preencheSistemaInterpolacao(sistema, tabela->x, funcoes);
        triangularizacao = triangularizacaoSistema(sistema);

        if(!triangularizacao) { // Se triangularização ocorreu certo, resolve
            solucao = resolveSistema(sistema, resultado);

            if(!solucao) { // Se resolveu corretamente, imprime
                imprimeResultado(resultado, n);
            }
        }
  
        LIKWID_MARKER_INIT;

        // Resolve utilizando ajuste de curvas
        LIKWID_MARKER_START("preencheSistemaAjuste");
            tempo = timestamp();
            preencheSistemaAjusteOtimizada(sistema, tabela->x, funcoes);
            tempo_p += timestamp() - tempo;
        LIKWID_MARKER_STOP("preencheSistemaAjuste");

        LIKWID_MARKER_START("triangularizacaoOtimizada");
            tempo = timestamp();
            triangularizacao = triangularizacaoSistemaOtimizada(sistema);
            tempo_t += timestamp() - tempo;
        LIKWID_MARKER_STOP("triangularizacaoOtimizada");

        LIKWID_MARKER_CLOSE;

        if(!triangularizacao) { // Se triangularizaçãoo ocorreu certo, resolve
            solucao = resolveSistema(sistema, resultado);

            if(!solucao) { // Se resolveu corretamente, imprime
                imprimeResultado(resultado, n);
            }
        }
    }

    printf("\n%d --- T Preenche Sistema: %lf", n, tempo_p/i);
    printf("\n%d --- T Triangularizacao: %lf", n, tempo_t/i);

    free(funcoes);
    free(resultado);
    liberaSistema(sistema);

    return 0;
}
