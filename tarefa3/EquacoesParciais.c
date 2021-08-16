#include "EquacoesParciais.h"

// Função B
real_t funcaoBUX1(real_t x) {
    return 20.0;
}

real_t funcaoBUX2(real_t x) {
    return 45.0;
}

real_t funcaoBUY1(real_t y) {
    return 0.0;
}

real_t funcaoBUY2(real_t y) {
    return 100.0;
}

real_t funcaoB(real_t hx, real_t hy) {
    return hx * hx * hy * hy * -1.0;
}

real_t funcaoBF(real_t x, real_t y) {
    return sin(x) * sin(x);
}

// Função D
real_t funcaoDUX1(real_t y) {
    return cos(y);
}

real_t funcaoDUX2(real_t y) {
    return -cos(y);
}

real_t funcaoDUY1(real_t x) {
    return cos(x);
}

real_t funcaoDUY2(real_t x) {
    return 0.0;
}

real_t funcaoD(real_t hx, real_t hy) {
    return 0.0;
}

real_t funcaoDF(real_t x, real_t y) {
    return -cos(x + y) - cos(x-y);
}

SLPentadiag_t* alocaSLPentadiag(int n, int m) {
    SLPentadiag_t *sl;
    sl = (SLPentadiag_t*)malloc(sizeof(SLPentadiag_t));

    int tam = n * m;

    sl->D = (real_t*)malloc(tam * sizeof(real_t));
    sl->Di = (real_t*)malloc(tam * sizeof(real_t));
    sl->Di2 = (real_t*)malloc(tam * sizeof(real_t));
    sl->Ds = (real_t*)malloc(tam * sizeof(real_t));
    sl->Ds2 = (real_t*)malloc(tam * sizeof(real_t));
    sl->B = (real_t*)malloc(tam * sizeof(real_t));
    for(int i = 0; i < tam; i++) {
        sl->D[i] = 0.0;
        sl->Di[i] = 0.0;
        sl->Di2[i] = 0.0;
        sl->Ds[i] = 0.0;
        sl->Ds2[i] = 0.0;
        sl->B[i] = 0.0;
    }
    sl->n = n;
    sl->m = m;
    
    return sl;
}

void liberaSLPentadiag(SLPentadiag_t *sl) {
    free(sl->D);
    free(sl->Di);
    free(sl->Di2);
    free(sl->Ds);
    free(sl->Ds2);
    free(sl->B);
    free(sl);
}

void imprimeSistemaPentaDiagonal(SLPentadiag_t *sistemaLinear) {
    int t = sistemaLinear->n * sistemaLinear->m;

    for(int i = 0; i < t - sistemaLinear->n; i++) {
        printf("%f ", sistemaLinear->Ds2[i]);
    }
    printf("\n");

    for(int i = 0; i < t - 1; i++) {
        printf("%f ", sistemaLinear->Ds[i]);
    }
    printf("\n");

    for(int i = 0; i < t; i++) {
        printf("%f ", sistemaLinear->D[i]);
    }
    printf("\n");

    for(int i = 1; i < t; i++) {
        printf("%f ", sistemaLinear->Di[i]);
    }
    printf("\n");

    for(int i = sistemaLinear->n; i < t; i++) {
        printf("%f ", sistemaLinear->Di2[i]);
    }
    printf("\n");

    for(int i = 0; i < t; i++) {
        printf("%f ", sistemaLinear->B[i]);
    }
    printf("\n");
}

void geraPentaDiagonal(Edp_t *edpeq, SLPentadiag_t *sl)
{
    int n = edpeq->n, m = edpeq->m, i, j, k = 0;
    real_t hx, hy, yi, xj, bi, d, di, di2, ds, ds2;

    hx = (edpeq->bX - edpeq->aX) / (n + 1);
    hy = (edpeq->bY - edpeq->aY) / (m + 1);

    for (i = 1; i <= m; i++) {
        yi = edpeq->aY + i * hy; // valor yi da malha

        for(j = 1; j <= n; j++) {
            xj = edpeq->aX + j * hx; // valor xj da malha
            bi = hx * hx * hy * hy * edpeq->f(xj, yi); // termo independente
            d = -2.0 * ((hx * hx) + (hy * hy)) + edpeq->u(hx, hy); // diagonal principal
            di = hy * hy; // diagonal inferior
            di2 = hx * hx; // diagonal inferior seg
            ds = hy * hy; // diagonal superior
            ds2 = hx * hx; // diagonal superior seg

            // Trata Y
            if(i == 1) {
                bi -= di2 * edpeq->uY1(yi);
                di2 = 0.0;
            }
            else if(i == m) {
                bi -= ds2 * edpeq->uY2(yi);
                ds2 = 0.0;
            }

            // Trata X
            if(j == 1) {
                bi -= di * edpeq->uX1(xj);
                di = 0.0;
            }
            else if(j == n) {
                bi -= ds * edpeq->uX2(xj);
                ds = 0.0;
            }

            sl->D[k] = d;
            sl->Di[k] = di;
            sl->Di2[k] = di2;
            sl->Ds[k] = ds;
            sl->Ds2[k] = ds2;
            sl->B[k] = bi;
            k++;
        }
    }
}

Edp_t* alocaEqParc() {
    Edp_t *eqParcial;
    eqParcial = (Edp_t*)malloc(sizeof(Edp_t));
    
    return eqParcial;
}

void liberaEqParc(Edp_t * eqParcial) {
    free(eqParcial);
}

void preencheEqParc(Edp_t *eqParcial, real_t aX, real_t bX, real_t aY, real_t bY, real_t (* uX1)(real_t),
                    real_t (* uX2)(real_t), real_t (* uY1)(real_t), real_t (* uY2)(real_t), 
                    real_t (* u)(real_t, real_t), real_t (* f)(real_t, real_t)) {
    eqParcial->aX = aX;
    eqParcial->bX = bX;
    eqParcial->aY = aY;
    eqParcial->bY = bY;
    eqParcial->uX1 = uX1;
    eqParcial->uX2 = uX2;
    eqParcial->uY1 = uY1;
    eqParcial->uY2 = uY2;
    eqParcial->u = u;
    eqParcial->f = f;
}

void executaFuncaoPentadiagonal(Edp_t* eqParcial, char item) {
    double tempo;
    int n = eqParcial->n;
    int m = eqParcial->m;
    real_t normaL2 = 0.0;

    real_t **solucao;
    solucao = (real_t**)malloc((m + 2) * sizeof(real_t*));
    for(int i = 0; i < (m + 2); i++) {
        solucao[i] = (real_t*)malloc((n + 2) * sizeof(real_t));
    }

    for(int i = 0; i < (m + 2); i++) {
        for(int j = 0; j < (n + 2); j++) {
            solucao[i][j] = 0.0;
        }
    }

    SLPentadiag_t *sistemaLinear;
    sistemaLinear = alocaSLPentadiag(n, m);

    printf("\n\n***** item (%c): L = %f, W = %f, n = %d, m = %d, Hx = %f, Hy = %f", item, eqParcial->bX, eqParcial->bY, n, m, (eqParcial->bX - eqParcial->aX) / (n + 1), (eqParcial->bY - eqParcial->aY) / (m + 1));
    
    geraPentaDiagonal(eqParcial, sistemaLinear);
    printf("\nSL:\n");
    imprimeSistemaPentaDiagonal(sistemaLinear);

    tempo = timestamp();
    gaussSeidelPentadiagonal(eqParcial, solucao);
    tempo = timestamp() - tempo;
    printf("\nT: ");
    imprimeResultadoPenta(solucao, n, m);

    normaL2 = normaL2ResiduoPenta(sistemaLinear, solucao);
    printf("Norma L2: %f, ", normaL2);
    printf("Tempo: %f ms\n", tempo);

    for(int i = 0; i < (m + 2); i++) {
        free(solucao[i]);
    }
    free(solucao);
    liberaSLPentadiag(sistemaLinear);
}

void calculaResiduoPenta(SLPentadiag_t *sl, real_t **solucao, real_t *residuo) {
    int n = sl->n;
    int m = sl->m;
    int count = 0;

    real_t *solucaoMapeada = NULL;
    solucaoMapeada = (real_t*)malloc(n * m * sizeof(real_t));

    for(int i = 1; i <= m; i++) {
        for(int j = 1; j <= n; j++) {
            solucaoMapeada[count] = solucao[i][j];
            count++;
        }
    }

    for(int i = 0; i < (n * m); i++) {
        real_t soma = 0.0;

        if(i == 0) {
            soma += sl->D[i] * solucaoMapeada[i];
            soma += sl->Ds[i] * solucaoMapeada[i + 1];
            soma += sl->Ds2[i] * solucaoMapeada[i + 2];
        }
        else if (i >= 1 && i <= 4) {
            soma += sl->Di[i - 1] * solucaoMapeada[i - 1];
            soma += sl->D[i] * solucaoMapeada[i];
            soma += sl->Ds[i] * solucaoMapeada[i + 1];
            soma += sl->Ds2[i] * solucaoMapeada[i + 2];
        }
        else if(i >= 5 && i <= 9) {
            soma += sl->Di2[i - 1] * solucaoMapeada[i - 2];
            soma += sl->Di[i - 1] * solucaoMapeada[i - 1];
            soma += sl->D[i] * solucaoMapeada[i];
            soma += sl->Ds[i] * solucaoMapeada[i + 1];
            soma += sl->Ds2[i] * solucaoMapeada[i + 2];
        }
        else if(i >= 10 && i <= 13) {
            soma += sl->Di2[i - 1] * solucaoMapeada[i - 2];
            soma += sl->Di[i - 1] * solucaoMapeada[i - 1];
            soma += sl->D[i] * solucaoMapeada[i];
            soma += sl->Ds[i] * solucaoMapeada[i + 1];
        }
        else {
            soma += sl->Di2[i - 1] * solucaoMapeada[i - 2];
            soma += sl->Di[i - 1] * solucaoMapeada[i - 1];
            soma += sl->D[i] * solucaoMapeada[i];
        }

        residuo[i] = sl->B[i] - soma;
    }

    free(solucaoMapeada);
}

real_t normaL2ResiduoPenta(SLPentadiag_t *sl, real_t **solucao) {
    int t = sl->n * sl->m;
    real_t normaL2 = 0.0;
    real_t *residuo = malloc(t * sizeof(real_t));

    calculaResiduoPenta(sl, solucao, residuo);

    for(int i = 0; i < t; i++) {
        normaL2 += pow(residuo[i], 2);
    }
    normaL2 = sqrt(normaL2);
    normaL2 = fabs(normaL2);

    free(residuo);

    return normaL2;
}

void gaussSeidelPentadiagonal(Edp_t *edpeq, real_t **Y)
{
    int n = edpeq->n, m = edpeq->m, k, i, j;
    real_t hx, hy, yi, xj, bi, d, di, di2, ds, ds2;

    hx = (edpeq->bX - edpeq->aX) / (n + 1);
    hy = (edpeq->bY - edpeq->aY) / (m + 1);

    d = -2.0 * ((hx * hx) + (hy * hy)) + edpeq->u(hx, hy); // diagonal principal
    di = hy * hy; // diagonal inferior
    di2 = hx * hx; // diagonal inferior seg
    ds = hy * hy; // diagonal superior
    ds2 = hx * hx; // diagonal superior seg


    for (k = 0; k < MAXIT; k++) {
        for (i = 1; i <= m; i++) {
            yi = edpeq->aY + i * hy; // valor yi da malha

            for(j = 1; j <= n; j++) {
                xj = edpeq->aX + j * hx; // valor xj da malha
                bi = hx * hx * hy * hy * edpeq->f(xj, yi); // termo independente

                // Trata Y
                if(i == 1) {
                    bi -= di2 * edpeq->uY1(yi) + ds2 * Y[i + 1][j];
                }
                else if(i == m) {
                    bi -= ds2 * edpeq->uY2(yi) + di2 * Y[i - 1][j];
                }
                else {
                    bi -= ds2 * Y[i + 1][j] + di2 * Y[i - 1][j];
                }

                // Trata X
                if(j == 1) {
                    bi -= di * edpeq->uX1(xj) + ds * Y[i][j + 1];
                }
                else if(j == n) {
                    bi -= ds * edpeq->uX2(xj) + di * Y[i][j - 1];
                }
                else {
                    bi -= ds * Y[i][j + 1] + di * Y[i][j - 1];
                }

                Y[i][j] = bi / d;
            }
        }
    }
}

void imprimeResultadoPenta(real_t **solucao, int n, int m) {
    for(int i = 1; i <= m; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%f ", solucao[i][j]);
        }
    }

    printf("\n");
}