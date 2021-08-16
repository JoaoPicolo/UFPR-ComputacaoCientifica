#include "comunicacaoIO.h"

/*!
    \brief Processa as informações da linha de comando.

    \param pivoteamento Ponteiro para variável que indica se deverá
                    ser realizado ou não o pivoteamento parcial.
    \param saida Ponteiro para char contendo o nome do arquivo no qual a saída
                    deverá ser escrita. NULL caso deva escrever na saída padrão (stdout).
    \param argc Quantidade de argumentos do programa.
    \param argv Vetor de ponteiros para cada string dos argumentos do programa.
*/
void processaLinhaComando(int *pivoteamento, char **saida, int argc, char **argv) {
    for(int i = 0; i < argc; i++) {
        if(strcmp(argv[i], "-p") == 0) {
            *pivoteamento = 1;
        }
        else if(strcmp(argv[i], "-o") == 0) {
            size_t tamanho = strlen(argv[i + 1]) + 1;
            *saida = (char*)malloc(tamanho * sizeof(char));
            strcpy(*saida, argv[i + 1]);
        }
    }

    // Garante que o arquivo de saída esteja vazio para posterior
    // impressão nele
    if(*saida != NULL) {
        FILE *fp;
        fp = fopen(*saida, "w");
        fclose(fp);
    }
}


/*!
    \brief Imprime resultado obtido a partir de uma matriz de entrada.

    \param matriz Estrutura definida em "matrizes.h", contendo as informações que deverão
                    ser impressas na tela.
    \param tTriangularizacao Tempo necessário para a triangularização, em milissegundos.
    \param tCalculoY Tempo médio para o calculo de todos os vetores y da fatoração LU,
                    em milissegundos.
    \param tCalculoX Tempo médio para o calculo de todos os vetores x da fatoração LU,
                    em milissegundos.
    \param arquivoSaida Ponteiro para char contendo o nome do arquivo no qual a saída
                    deverá ser escrita. NULL caso deva escrever na saída padrão (stdout).
*/
void imprimeResultado(Matriz_t *matriz, double *tTriangularizacao, 
                double *tCalculoY, double* tCalculoX, char *arquivoSaida) {
    
    FILE *fp = stdout;
    if(arquivoSaida != NULL)
        fp = fopen(arquivoSaida, "a");

    // Imprime ordem da matiz
    unsigned int n = matriz->n;
    fprintf(fp, "%d\n", n);

    // Imprime matriz de entrada
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(fp, "%f ", matriz->A[i][j]);
        }
        fprintf(fp, "\n");
    }

    // Imprime matriz inversa
    fprintf(fp, "#\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(fp, "%f ", matriz->AI[i][j]);
        }
        fprintf(fp, "\n");
    }

    // Imprime informações de tempo
    fprintf(fp, "############\n");
    fprintf(fp, "# Tempo Triangularização: %f\n", *tTriangularizacao);
    fprintf(fp, "# Tempo cálculo de Y: %f\n", *tCalculoX);
    fprintf(fp, "# Tempo cálculo de X: %f\n", *tCalculoY);
    fprintf(fp, "# Norma L2 do Residuo: ");
    for(int i = 0; i < n; i++) {
        fprintf(fp, "%f ", matriz->normas[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "############\n\n");

    if(arquivoSaida != NULL)
        fclose(fp);
}