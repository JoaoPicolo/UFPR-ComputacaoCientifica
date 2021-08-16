#include "EquacoesOrdinarias.h"
#include "EquacoesParciais.h"
#include "utils.h"

int main() {
    Edo_t *eqOrdinaria;
    eqOrdinaria = alocaEqOrd();

    Edp_t *eqParcial;
    eqParcial = alocaEqParc();

    // Função A
    eqOrdinaria->n = 5;
    preencheEqOrd(eqOrdinaria, 0.0, 12.0, 0.0, 0.0, funcaoAP, funcaoAQ, funcaoAR);
    executaFuncaoTridiagonal(eqOrdinaria, 'a');
    eqOrdinaria->n = 10;
    executaFuncaoTridiagonal(eqOrdinaria, 'a');

    // Função B
    eqParcial->m = 3;
    eqParcial->n = 5;
    preencheEqParc(eqParcial, 0.0, 6.0, 0.0, 8.0, funcaoBUX1, funcaoBUX2, funcaoBUY1, funcaoBUY2, funcaoB, funcaoBF);
    executaFuncaoPentadiagonal(eqParcial, 'b');
    eqParcial->n = 10;
    executaFuncaoPentadiagonal(eqParcial, 'b');

    // Função C
    eqOrdinaria->n = 5;
    preencheEqOrd(eqOrdinaria, 0.0, 1.0, 0.0, 1.0, funcaoCP, funcaoCQ, funcaoCR);
    executaFuncaoTridiagonal(eqOrdinaria, 'c');
    eqOrdinaria->n = 10;
    executaFuncaoTridiagonal(eqOrdinaria, 'c');

    // Função D
    eqParcial->m = 3;
    eqParcial->n = 5;
    preencheEqParc(eqParcial, 0.0, M_PI, 0.0, M_PI / 2.0, funcaoDUX1, funcaoDUX2, funcaoDUY1, funcaoDUY2, funcaoD, funcaoDF);
    executaFuncaoPentadiagonal(eqParcial, 'd');
    eqParcial->n = 10;
    executaFuncaoPentadiagonal(eqParcial, 'd');

    liberaEqParc(eqParcial);
    liberaEqOrd(eqOrdinaria);
}