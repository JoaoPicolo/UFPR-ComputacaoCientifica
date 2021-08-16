#include "utils.h"

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/
double timestamp(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

/*  Retorna se linha lida é vazia */
int linhaVazia(char *linha) {
    if(linha[0] == '\n') {
        return 1;
    }

    return  0;
}

real_t** copiaMatrizCoeficientes(real_t **matCoeficientes, unsigned int n) {
    real_t **copia;

    copia = malloc(n * sizeof(real_t*));
    memcpy(copia, matCoeficientes, n * sizeof(real_t*));

    for(int i = 0; i < n; i++) {
        copia[i] = malloc(n * sizeof(real_t));
        memcpy(copia[i], matCoeficientes[i], n * sizeof(real_t));
    }

    return copia;
}

real_t* copiaVetorTermos(real_t *vetTermos, unsigned int n) {
    real_t *copia;

    copia = malloc(n * sizeof(real_t));
    memcpy(copia, vetTermos, n * sizeof(real_t));

    return copia;
}

void liberaCopiasLocais(real_t **matCoeficientes, real_t *vetTermos, unsigned int n) {
    // Desaloca Matriz de Coeficientes
    for(int i = 0; i < n; i++) {
        free(matCoeficientes[i]);
    }
    free(matCoeficientes);

    // Desaloca Vetor de Termos Independentes
    free(vetTermos);
}

unsigned int encontraMax(real_t **matCoeficientes, unsigned int n, unsigned int col) {
    unsigned int indexMaior = col;

    for(unsigned int i = col + 1; i < n; i++) {
        if(fabs(matCoeficientes[i][col]) > fabs(matCoeficientes[indexMaior][col])) {
            indexMaior = i;
        }
    }

    return indexMaior;
}

void trocaLinha(real_t **matCoeficientes, real_t *vetTermos, 
    unsigned int n, unsigned int i, unsigned int iPivo) {
    
    // Troca linhas na matriz
    real_t aux;
    for(int j = 0; j < n; j++) {
        aux = matCoeficientes[i][j];
        matCoeficientes[i][j] = matCoeficientes[iPivo][j];
        matCoeficientes[iPivo][j] = aux;
    }

    // Troca linhas no vetor de termos indepedentes
    aux = vetTermos[i];
    vetTermos[i] = vetTermos[iPivo];
    vetTermos[iPivo] = aux;
}

void retrosubstituicao(real_t **matCoeficientes, real_t *x, real_t *vetTermos, unsigned int n) {
    for(int i = n - 1; i >= 0; i--) {
        real_t soma = 0.0;
        for(int j = i + 1; j < n; j++) {
            soma += matCoeficientes[i][j] * x[j];
        }

        x[i] = (vetTermos[i] - soma) / matCoeficientes[i][i];
    }
}

void calculaResiduo(real_t **matCoeficientes, real_t *xRefinado,
    real_t *vetTermos, real_t *residuo, unsigned int n) {

    for(int i = 0; i < n; i++) {
        real_t soma = 0.0;

        for(int j = 0; j < n; j++) {
            soma += matCoeficientes[i][j] * xRefinado[j];
        }

        residuo[i] = vetTermos[i] - soma;
    }

}

void inicializaVetAproximacao(real_t *x, unsigned int n) {
    for(int i = 0; i < n; i++) {
        x[i] = 0.0;
    }
}

int menorQueErro(real_t *vetInicial, real_t *x, real_t erro, unsigned int n, int it) {
    if(it > 0) {
        real_t maior, diff;
        maior = fabs(x[0] - vetInicial[0]);

        for(int i = 1; i < n; i++) {
            diff = fabs(x[i] - vetInicial[i]);
            if(diff > maior) {
                maior = diff;
            }
        }

        if(maior <= erro) {
            return 1;
        }
    }

    return 0;
}

int sistemaConverge(real_t **matCoeficientes, unsigned int n) {
    for(int i = 0; i < n; i++) {
        real_t elemento = fabs(matCoeficientes[i][i]);
        real_t soma = 0.0;
        for(int j = 0; j < n; j++) {
            if(i != j) {
                soma += fabs(matCoeficientes[i][j]);
            }
        }

        if(elemento < soma) {
            return 0;
        }
    }

    return 1;
}

void aplicaRegrasResiduo(real_t *residuo, RegrasAplicadas_t *regras, int n) {
    real_t aux;

    for(int i = 0; i < n; i++) {
        RegrasAplicadas_t regra = regras[i];
        if(regra.tipo == 0) { // Troca de linha
            aux = residuo[regra.linha1];
            residuo[regra.linha1] = residuo[regra.linha2];
            residuo[regra.linha2] = aux;
        }
        else { // Operacao matematica
            residuo[regra.linha1] -= residuo[regra.linha2] * regra.m;
        }
    }
}

int possuiSolucao(real_t **matCoeficientes, real_t *vetTermos, int n) {
    for(int i = 0; i < n; i++) {
        int isAllZero = 1;
        for(int j = 0; j < n; j++) {
            if(matCoeficientes[i][j] != 0.0) {
                isAllZero = 0;
            }
        }

        if(isAllZero && (vetTermos[i] != 0.0)) {
            return 0;
        }
    }

    return 1;
}
