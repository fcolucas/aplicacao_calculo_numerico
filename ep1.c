#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <string.h>
#define ERRO 0.00000001

//Variáveis Globais
double **M; //armazena a matriz
double *p;  //armazena a matriz de ponteiros
double *x; //armazena a solução de um SL
int n;  //armazena a ordem da matriz
int *val; //armazena os coeficientes de um SL - usado no método de Jordan

double *alocaFuncao(int n){
    /*A função lê valores para um vetor de double, alocada
    dinamicamente, que são coeficientes de uma equação algé-
    brica de grau n.*/
    int i;
    double *v = malloc(sizeof(double)*n);
    if (v != NULL){
        for(i=n-1; i>=0; i--){
            printf("Digite o coeficiente a%d: ", i);
            scanf("%lf", &v[i]);
        }
        return v;
    }
    else return NULL;
}

double **alocaMatriz(int l, int c){
	/*Se houver memoria disponivel, aloca
	uma matriz de double com l linhas e c
	colunas e devolve um ponteiro para a matriz.
	Caso contrario, devolve o ponteiro NULL.*/
	double **M;
	int i, j;
	M = malloc(sizeof(double *)*l);
	if (M == NULL) return NULL; /*falta de memória*/
	for (i=0; i<l; i++){
		M[i] = malloc(sizeof(double)*c);
		if(M[i] == NULL){ /*falta de memória*/
			for(j=0; j<i; j++){
				free(M[j]); }
			free(M);
			return NULL;
		}
	}
	return M;
} /*Fim alocaMatriz()*/

void leMatriz(){
    /*A função lê um arquivo de texto contendo um sistema linear de n equações e n variáveis */
    int i, j;
	FILE *arquivo = NULL;
	char caminhoDoArquivo[50];

	while(arquivo == NULL){
		printf("Insira o caminho do arquivo e pressione Enter: \n");
		fgets(caminhoDoArquivo, sizeof(caminhoDoArquivo), stdin);
		//Remove o último caractere do caminho, pois o fgets armazena a quebra de linha '\n'
		char *p_chr = strchr(caminhoDoArquivo, '\n');
		if(p_chr != NULL)
			*p_chr = '\0';
		arquivo = fopen(caminhoDoArquivo, "r");
	}

	fscanf(arquivo, "%d", &n); //Lê a primeira linha do arquivo

	M = alocaMatriz(n, n + 1);
	p = malloc(sizeof(double)*n); //Matriz de ponteiros
	if(M == NULL || p == NULL)
		printf("Memória insuficiente!");

	for(i = 0; i < n; i++){
		for(j = 0; j < n + 1; j++){
			fscanf(arquivo, "%lf", &M[i][j]);
		}
	}
	fclose(arquivo);
}//lerMatriz

void imprimeMatriz(double **M, int l, int c){
	/*Imprime o conteudo de uma matriz de double,
	alocada dinamicamente, com l linhas e c colunas. */
	int i, j;
	for(i=0; i<l; i++){
		for(j=0; j<c; j++){
			printf("%10.3lf ", M[i][j]);
		}
		printf("\n");
	}
}/*fim imprimeMatriz()*/

int *alocaVariavel(int n){
    /*A função aloca um vetor de variáveis de um sistema linear.
    */
    int i, *val = malloc(sizeof(int)*n);
    for(i=0; i<n; i++){
        val[i] = i+1;
    }
    return val;
}

void conversao(double num){
	/*A função converte um número no sistema decimal para os sistemas
	binário, octal e hexadecimal. */
    int i = 0;  //variavel auxiliar
	int base[3] = {2,8,16};	 //Array de possibilidades
	int precisao = 20;

	for(i=0 ; i<3 ; i++){

		int quociente = (int)num, aux = 0, j = 0;   //variavel usada nas iteracoes da parte inteira
		double parteInteira = 0;
		double parteFracionaria = num-(double) quociente;

		if((base[i] == 2) || (base[i] == 8)){
			while(quociente >= base[i]){
				parteInteira += quociente%base[i]*pow(10, j);
				quociente = quociente/base[i];
				j++;
			}
			parteInteira += quociente%base[i]*pow(10, j);
			if(base[i] == 2){
				printf("Binário: %d.", (int)parteInteira);
			}else{
				printf("Octal: %d.", (int)parteInteira);
			}
		}

		else if(base[i] == 16){
			parteFracionaria = num-(int)num;
			printf("Hexadecimal: %X.", (int)num);
		}

		for(j = 1; j <= precisao; j++){    //Imprime a parte fracionária
			parteFracionaria *= base[i];
			aux = (int)parteFracionaria;

			if(parteFracionaria > 0){
				if((base[i] == 2) || (base[i] == 8))
					printf("%d", (int)parteFracionaria);
				else if(base[i] == 16)
					printf("%X", (int)parteFracionaria);

				parteFracionaria -= aux;
			}
			else
				printf("0");

			if(parteFracionaria == 0)
				break;
		}
		printf("\n");
	}
}

int sretro(double **M, int n){
	/* Algoritmo de substituição retroativa.
	Recebe m, a matriz aumentada de um SL TS com n variaveis.
	Se o sistema linear for determinado, armazena em x a solução
	no SL e devolve 0.
	Se for indeterminado, armazena em x uma solução do SL e
	devolve 1. */
	int i,j, tipo=0;
	double soma;
	for(i=n-1; i>=0; i--){
		soma=0;
		for(j=i+1; j<n; j++){
			soma += M[i][j] * x[val[j]-1];
		}
		if(M[i][i] == 0){
			if(fabs(M[i][n] - soma) < ERRO /*m[i][n] == soma*/){
				x[i]=0;	/*variavel livre*/
				tipo = 1;
			}else{
				return 2; /*SL incompativel*/
			}
		}else{
			x[val[i]-1] = (M[i][n] - soma)/M[i][i];
			}
	}
	return tipo;
} /*Fim sretro()*/

void solucao(double **M, int n, double *x){
    /*A função fornece a solução de um sistema linear.*/
    int i, tipo = sretro(M, n);
	if(tipo == 2){
		printf("\nSL imcompativel!\n");
	}else{
		printf("\nSL %sdeterminado\n", tipo?"In":"");
		for(i=0; i<n; i++){
			printf("x[%d] = %10.3lf\n", i+1, x[i]);
		}
	}
}

void metJordan(double **M, int n){
    /*A função aplica o metodo de Jordan na matriz M, ou seja, transforma a matriz em
    uma matriz diagonal. No final, calcula a solução desse istema.
    */
    int l, c, k, i, j;
    double mult, aux, num;
    x = malloc(sizeof(double)*n);
    val = alocaVariavel(n);

    //Verificando se os elementos da diagonal é 0!
    for(j=0; j<n; j++){
        if(fabs(M[j][j]) < ERRO){
            c = j+1;
            while(c < n && fabs(M[j][c]) < ERRO){
                c++;
            }
            if(c < n){
                num = val[j];
                val[j] = val[c];
                val[c] = num;
                for(i = 0; i < n; i++){
                    aux = M[i][j];
                    M[i][j] = M[i][c];
                    M[i][c] = aux;
                }
            }
            else{
                for(i = 0; i < n; i++){
                    M[i][j] = 0;
                }
            }
        }
        if(M[j][j] != 0){
            for(l = 0; l<n; l++){
                if(l != j){
                    mult =  (- M[l][j])/M[j][j];
                    M[l][j] = 0;
                    for(k = j+1; k<=n; k++){
                        M[l][k] += (mult * M[j][k]);
                        if(fabs(M[l][k]) < ERRO) M[l][k] = 0;
                    }
                }
            }
        }
    }
    printf("Matriz triangularizada - Metodo de Jordan \n");
    imprimeMatriz(M, n, n+1);
    solucao(M, n, x);
}

double calculaFuncao(int n, double v[], double x){
    /*A função retorna o resultado de uma funcao em um ponto x, ou
    seja, seja f a função e x o ponto, a função calcula f(x).*/
    int i;
    /* n = grau do polinomio
        v = vetor de coeficientes
        x = variavel da equação*/
    double resp = 0;
    for(i=n; i>=0; i--){
        resp += v[i]*pow(x, i);
    }
    return resp;
}

void metBissecao(int n, double *v){
    /*A função aplica o método da bisseção em uma equação algébrica
    dado um intervalo de observação.*/
    int bolzano(double a, double b, int n, double *v){
        /*A função verifica se entre os intervalos dados, existe
        uma quantidade ímpar de raízes numa equação algébrica. Se
        caso positivo, retorna 1. Se não, retorna 0.*/
        double Fa, Fb;
        Fa = calculaFuncao(n, v, a);
        Fb = calculaFuncao(n, v, b);
        if(Fa * Fb > 0){
            return 0;
        }
        else return 1;
    }

    double a, b, m, erro, Fa, Fb, Fm = 1, i = 0;
    printf("Defina o intervalo: ");
    scanf("%lf%lf", &a, &b);
    fflush(stdin);
    erro = (b - a)/2;
    if(bolzano(a, b, n, v)){
        printf("  a       b      m      f(a)     f(b)      f(m)      Erro    \n");
        printf("-----------------------------------------------------\n");
        while(i < 1000 && erro > ERRO && Fm != 0){
            m = (b + a)/2;
            erro = (b - a)/2;
            Fa = calculaFuncao(n, v, a);
            Fb = calculaFuncao(n, v, b);
            Fm = calculaFuncao(n, v, m);
            printf("%7.8lf %7.8lf %7.8lf %7.8lf %7.8lf %8.8lf %8.8lf \n", a, b, m, Fa, Fb, Fm, erro);
            if(Fa * Fm > 0) a = m;
            else if (Fm == 0) printf("Raiz da equacao = %lf\n", Fm);
            else b = m;
            i++;
        }
        printf("Raiz aproximada do polinomio = %7.8lf \n", m);
    }
    else printf("O numero de raízes no intervalo [%lf, %lf] e par\n", a, b);
}

int main()
{
    setlocale(LC_ALL, "Portuguese");
    int menu = 1;
    char op;
    double *v, x, resp, num;
    while(menu){
        printf("\n-----EP1-----\n");
        printf("C - Conversao\n");
        printf("S - Sistema Linear\n");
        printf("E - Equacao Algebrica\n");
        printf("F - Finalizar\n");

        printf("-> Digite sua opcao: ");
        scanf("%c", &op);
        fflush(stdin);

        switch(op){
            case 'C':
                printf("\n***Conversao***\n");
                printf("-> Digite o numero a ser convertido: ");
                scanf("%lf", &num);
                fflush(stdin);
                conversao(num);
                break;

            case 'S':
                printf("\n***Sistema Linear***\n");
                leMatriz();
                printf("Matriz: \n");
                imprimeMatriz(M, n, n+1);
                printf("\n");
                metJordan(M, n);
                break;

            case 'E':
                printf("\n***Equacao Algebrica***\n");
                printf("Digite o grau do polinomio: ");
                scanf("%d", &n);
                v = alocaFuncao(n+1);
                metBissecao(n, v);
                break;

            case 'F':
                printf("Finalizando...\n");
                menu = 0;
                break;

            default:
                printf("Opcao Invalida!\n");
                break;
        }
    }
    return 0;
}
