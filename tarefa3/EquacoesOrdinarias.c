#include "utils.h"
#include "EquacoesOrdinarias.h"

// Função A
real_t funcaoAP(real_t x) {
    return 0.0;
}

real_t funcaoAQ(real_t x) {
    return 0.0;
}

real_t funcaoAR(real_t x) {
    return (6 * x - 0.5 * x * x);
}

// Função C
real_t funcaoCP(real_t x) {
    return 0.0;
}

real_t funcaoCQ(real_t x) {
    return 1.0;
}

real_t funcaoCR(real_t x) {
    return 0.0;
}

SLTridiag_t* alocaSLTridiag(int n) {
    SLTridiag_t *sl;
    sl = (SLTridiag_t*)malloc(sizeof(SLTridiag_t));

    sl->D = (real_t*)malloc(n * sizeof(real_t));
    sl->Di = (real_t*)malloc(n * sizeof(real_t));
    sl->Ds = (real_t*)malloc(n * sizeof(real_t));
    sl->B = (real_t*)malloc(n * sizeof(real_t));
    for(int i = 0; i < n; i++) {
        sl->D[i] = 0.0;
        sl->Di[i] = 0.0;
        sl->Ds[i] = 0.0;
        sl->B[i] = 0.0;
    }
    sl->n = n;
    
    return sl;
}

void liberaSLTridiag(SLTridiag_t *sl) {
    free(sl->D);
    free(sl->Di);
    free(sl->Ds);
    free(sl->B);
    free(sl);
}

void imprimeSistemaTriDiagonal(SLTridiag_t *sistemaLinear) {
    int n = sistemaLinear->n;

    for(int i = 0; i < n - 1; i++) {
        printf("%f ", sistemaLinear->Ds[i]);
    }
    printf("\n");

    for(int i = 0; i < n; i++) {
        printf("%f ", sistemaLinear->D[i]);
    }
    printf("\n");

    for(int i = 1; i < n; i++) {
        printf("%f ", sistemaLinear->Di[i]);
    }
    printf("\n");

    for(int i = 0; i < n; i++) {
        printf("%f ", sistemaLinear->B[i]);
    }
    printf("\n");
}

void geraTriDiagonal(Edo_t *edoeq, SLTridiag_t *sl)
{
    real_t xi, h;
    h = (edoeq->b - edoeq->a) / (edoeq->n + 1.0);
    
    for (int i=0; i < edoeq->n; ++i) {
        xi = edoeq->a + (i + 1) * h; // ponto da malha
        sl->Di[i] = 1 - h * edoeq->p(xi)/2.0; // diagonal inferior
        sl->D[i] = -2 + h * h * edoeq->q(xi); // diagonal principal
        sl->Ds[i] = 1 + h * edoeq->p(xi)/2.0; // diagonal superior
        sl->B[i] = h * h * edoeq->r(xi); // termo independente
    }
    
    // Condições de contorno subtraídas do 1º e último termos independentes
    sl->B[0] -= edoeq->ya * (1 - h * edoeq->p(edoeq->a + h) / 2.0);
    sl->B[edoeq->n-1] -= edoeq->yb * (1 + h * edoeq->p(edoeq->b - h) / 2.0);
}

Edo_t* alocaEqOrd() {
    Edo_t *eqOrdinaria;
    eqOrdinaria = (Edo_t*)malloc(sizeof(Edo_t));
    
    return eqOrdinaria;
}

void liberaEqOrd(Edo_t * eqOrdinaria) {
    free(eqOrdinaria);
}

void preencheEqOrd(Edo_t *eqOrdinaria, real_t a, real_t b, real_t ya, real_t yb, real_t (*p)(real_t), 
                    real_t (*q)(real_t), real_t (*r)(real_t)) {
    eqOrdinaria->a = a;
    eqOrdinaria->b = b;
    eqOrdinaria->ya = ya;
    eqOrdinaria->yb = yb;
    eqOrdinaria->p = p;
    eqOrdinaria->q = q;
    eqOrdinaria->r = r;
}

void executaFuncaoTridiagonal(Edo_t* eqOrdinaria, char item) {
    double tempo;
    int n = eqOrdinaria->n;
    real_t normaL2 = 0.0;

    real_t *solucao;
    solucao = (real_t*)malloc(n * sizeof(real_t));
    for(int i = 0; i < n; i++) {
        solucao[i] = 0.0;
    }

    SLTridiag_t *sistemaLinear;
    sistemaLinear = alocaSLTridiag(n);

    printf("\n\n***** item (%c): n = %d, H = %f", item, n, (eqOrdinaria->b - eqOrdinaria->a) / (n + 1));

    geraTriDiagonal(eqOrdinaria, sistemaLinear);
    printf("\nSL:\n");
    imprimeSistemaTriDiagonal(sistemaLinear);

    tempo = timestamp();
    gaussSeidel(eqOrdinaria, solucao);
    tempo = timestamp() - tempo;
    printf("\nY: ");
    imprimeResultadoTri(solucao, n);

    normaL2 = normaL2ResiduoTri(sistemaLinear, solucao);
    printf("Norma L2: %f, ", normaL2);
    printf("Tempo: %f ms\n", tempo);

    free(solucao);
    liberaSLTridiag(sistemaLinear);
}

void calculaResiduoTri(SLTridiag_t *sl, real_t *solucao, real_t *residuo) {
    int n = sl->n;

    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;

        if(i == 0) {
            soma += sl->D[i] * solucao[i];
            soma += sl->Ds[i] * solucao[i + 1];
        }
        else if(i == n - 1) {
            soma += sl->D[i] * solucao[i];
            soma += sl->Di[i - 1] * solucao[i - 1];
        }
        else {
            soma += sl->Ds[i] * solucao[i + 1];
            soma += sl->D[i] * solucao[i];
            soma += sl->Di[i - 1] * solucao[i - 1];
        }

        residuo[i] = sl->B[i] - soma;
    }
}

real_t normaL2ResiduoTri(SLTridiag_t *sl, real_t *solucao) {
    real_t normaL2 = 0.0;
    real_t *residuo = malloc(sl->n * sizeof(real_t));

    calculaResiduoTri(sl, solucao, residuo);

    for(int i = 0; i < sl->n; i++) {
        normaL2 += pow(residuo[i], 2);
    }
    normaL2 = sqrt(normaL2);
    normaL2 = fabs(normaL2);

    free(residuo);

    return normaL2;
}

void gaussSeidel(Edo_t *edoeq, real_t *Y)
{
    int n = edoeq->n, k, i;
    real_t h, xi, bi, d, di, ds;

    h = (edoeq->b - edoeq->a) / (n + 1); // Largura do passo da malha

    for (k = 0; k < MAXIT; ++k) {
        for (i = 0; i < n; ++i) { // Para cada equação do SL
            xi = edoeq->a + (i + 1) * h; // valor xi da malha
            bi = h * h * edoeq->r(xi); // termo independente
            di = 1 - h * edoeq->p(xi) / 2.0; // diagonal inferior
            d = -2 + h * h * edoeq->q(xi); // diagonal principal
            ds = 1 + h * edoeq->p(xi) / 2.0; // diagonal superior
            
            if (i == 0) bi -= ds * Y[i + 1] + edoeq->ya * (1 - h * edoeq->p(edoeq->a + h) / 2.0);
            else if (i == n-1) bi -= di * Y[i - 1] + edoeq->yb * (1 + h * edoeq->p(edoeq->b - h) / 2.0);
            else bi -= ds * Y[i + 1] + di * Y[i - 1];

            Y[i] = bi / d;
        }
    }
}

void imprimeResultadoTri(double *solucao, int n) {
    for(int i = 0; i < n; i++) {
        printf("%f ", solucao[i]);
    }

    printf("\n");
}