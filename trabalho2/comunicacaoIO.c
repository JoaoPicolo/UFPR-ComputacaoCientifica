// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "comunicacaoIO.h"

/*!
    \brief Imprime o vetor de resultados passado por parâmetro na saída padrão (stdout).

    \param resultado Ponteiro para vetor do tipo real_t (definido em "utils.h") contendo
                    os valores resultantes das operações aplicadas.
    \param n Tamanho do vetor de resultado.
*/
void imprimeResultado(real_t *resultado, unsigned int n) {
    for(int i = 0; i < n; i++)
        printf("%f ", resultado[i]);

    printf("\n");
}