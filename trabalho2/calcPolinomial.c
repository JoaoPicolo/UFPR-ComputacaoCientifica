// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "utils.h"
#include "tabela.h"

int main(int argc, char** argv) {
    Tabela_t *tabela;
    tabela = lerTabela();

    // Só executa se não houve erro de alocação
    if(tabela != NULL) {
        executaFuncoesPolinomiais(tabela);
        liberaTabela(tabela);
    }
}