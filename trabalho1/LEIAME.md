# Trabalho 01 - Introdução a Computação Científica

## Alunos

- Gabriel Marckzuk Thá - GRR20186070 
- João Pedro Picolo    - GRR20182659

## Descrição

O presente trabalho busca solucionar o problema do **cálculo da matriz inversa** a partir de uma matriz de entrada através do **Método de Fatoração LU**.

Neste documento procuraremos dar uma breve explicação sobre a execução do programa e as estruturas de dados utilizadas. A descrição de cada função implementada pode ser encontrada diretamente no código fonte.


## Compilação e Execução

Para compilar o programa basta utilizar o comando make.

```bash=
make
```

Após isto, será gerado o executável *matrizInv*, que poderá ser executado com os argumentos presentes na [especificação do trabalho](https://moodle.c3sl.ufpr.br/mod/assign/view.php?id=28417).

Para a realização de testes, os arquivo *exemplo.txt* foi colocado sob o diretório *testes*. Desta forma temos as seguintes opções para execução:

- Execução sem pivoteamento com escrita na saída padrão:

```bash=
./matrizInv < testes/exemplos.txt
```

- Execução com pivoteamento com escrita na saída padrão:

```bash=
./matrizInv -p < testes/exemplos.txt
```

- Execução sem pivoteamento com escrita no arquivo *out.txt* :

```bash=
./matrizInv -o out.txt < testes/exemplos.txt
```

- Execução com pivoteamento com escrita no arquivo *out.txt* :

```bash=
./matrizInv -p -o out.txt < testes/exemplos.txt
```

Vale ressaltar que os exemplos fornecidos seguem a entrada descrita na [especificação do trabalho](https://moodle.c3sl.ufpr.br/mod/assign/view.php?id=28417), assim como a saída fornecida pelo programa.


## Estrutura de dados

Para a execução do programa foi desenvovilda uma estrutura de dados chamada *Matriz_t* que guarda todas as informações necessárias para a execução do programa. Definida desta forma:

```c=
typedef float real_t;

typedef struct {
    unsigned int n;
    int *ordemLinhas;
    real_t **A, **AI;
    real_t **L, **U;
    real_t *normas;
} Matriz_t;
```

Para esta estrutura temos a seguinte configuração:

- n: ordem das matrizes (entrada, inversa e identidade);

- ordemLinhas: vetor auxiliar de tamaho *n* que armazena em seu i-ésimo índice a posição do valor 1 na i-ésima coluna da matriz identidade. Escolheu-se esta estrutura para a representação da matriz identidade por acreditar-se ser mais eficiente sua alocação na memória. Por exemplo, sejam n = 3 e ordemLinhas = [1 , 2, 0] a matriz identididade será:

$$
 \begin{bmatrix}
   0 & 0 & 1 \\
   1 & 0 & 0 \\
   0 & 1 & 0
  \end{bmatrix} \
$$

- A: matriz de entrada;

- AI: matriz inversa calculada;

- L e U: matrizes L e U da fatoração LU. Assim como no caso do vetor ordemLinhas, as estruturas escolhidas aqui diferem de uma matriz padrão por questões de eficiência relacionadas a memória. Por exemplo:

Seja A a matriz de entrada:

$$
 \begin{bmatrix}
   25 & 5 & 1 \\
   64 & 8 & 1 \\
   144 & 12 & 1
  \end{bmatrix} \
$$

Se realizarmos pivoteamento sobre essa matriz: obteremos primeramente a matriz identidade descrita acima. Também obteremos a matriz L:

$$
 \begin{bmatrix}
   1 &  &  \\
   0.1736 & 1 &  \\
   0.4444 & 0.9144 & 1
  \end{bmatrix} \
$$

Note que nesta estrutura não armazenamos valores que não serão utilizados no cálculo da matriz inversa. Comportamento semelhante ocorre na matriz U:

$$
 \begin{bmatrix}
   144 & 12 & 1 \\
    & 2.9168 & 0.8264 \\
    &  & -0.2001
  \end{bmatrix} \
$$

- normas: vetor auxiliar de tamanho *n* que armazena em sua i-ésima posição a normal L2 do resíduo da matriz *A* com a i-ésima coluna da matriz *AI* em relação a i-ésima coluna da matriz identidade.


## Algoritmo

Como citado anteriormente, para a resolução do problema utilizou-se o **Método de Fatoração LU**.

Dada a matriz de entrada aplica-se primeiramente o Método de Gauss (com ou sem pivoteamento parcial) como forma de construir as matrizes L e U necessárias.

Após obtidas, usamos o método da retrosubstituição para calcular os vetore *y* a partir de L e os vetores *x* a partir de U. Podemos utilizar este método pois ambas as matrizes são triangulares, ver exemplos da seção "Estrutura de Dados".
