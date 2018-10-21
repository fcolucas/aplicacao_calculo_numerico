#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <string.h>
#include <ctype.h>
#define ERRO 0.00000001

//Variáveis globais
double **M; //armazena a matriz
double *x; //armazena a solução de um SL
int n;  //armazena a ordem da matriz
int *val; //armazena os coeficientes de um SL - usado no método de Jordan
double *funcao; //armazena os coeficientes de uma euação algébrica
int grau; //armazena o grau de uma equação algébrica
double num; //armazena o número na base decimal a ser convertido

double *alocaFuncao(int grau){
    /*A função lê valores para um vetor de double, alocada
    dinamicamente, que são coeficientes de uma equação algé-
    brica de grau n.*/
    int i;
    double *v = malloc(sizeof(double)*grau);
    if (v != NULL){
        for(i=grau-1; i>=0; i--){
            printf("Digite o coeficiente a%d: ", i);
            scanf("%lf", &v[i]);
        }
        return v;
    }
    else return NULL;
}//alocaFuncao()

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
}//alocaMatriz()

void leMatriz(double **M, int n){
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

	M = alocaMatriz(n, n+1);
	if(M == NULL) printf("Memória insuficiente!");
	else {
		for(i = 0; i < n; i++){
			for(j = 0; j < n + 1; j++){
				fscanf(arquivo, "%lf", &M[i][j]);
			}
		}
	}
	fclose(arquivo);
}//leMatriz()

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
}//imprimeMatriz()

int *alocaVariavel(int n){
    /*A função aloca um vetor de variáveis de um sistema linear.
    */
    int i;
    val = malloc(sizeof(int)*n);
    for(i=0; i<n; i++){
        val[i] = i+1;
    }
    return val;
}//alocaVariavel()

void conversao(double num){
	/*A função converte um número no sistema decimal para os sistemas
	binário, octal e hexadecimal. */
    int i = 0, j = 1;  //variavel auxiliar
	int base[3] = {2,8,16};	 //Array de possibilidades
	int precisao = 20;

	for(i=0 ; i<3 ; i++){
		int quociente = (int)num, aux = 0, j = 0;   //variavel usada nas iteracoes da parte inteira
		double parteInteira = 0, parteFracionaria = num - (double)quociente;

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

		while(j <= precisao && parteFracionaria != 0){    //Imprime a parte fracionária
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
			j++;
		}
	}
}//conversao()

int sretro(double **M, int n){
	/* Algoritmo de substituição retroativa.
	Recebe m, a matriz aumentada de um SL TS com n variaveis.
	Se o sistema linear for determinado, armazena em x a soluÃ§Ã£o
	no SL e devolve 0.
	Se for indeterminado, armazena em x uma solução do SL e
	devolve 1. */
	int i, j, tipo = 0;
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
}//sretro()

void solucaoSistema(double **M, int n, double *x){
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
}//solucaoSistema()

void metJordan(double **M, int n){
    /*A função aplica o metodo de Jordan na matriz M, ou seja,
	transforma a matriz em uma matriz diagonal. No final,
	calcula a solução desse SL.
    */
    int l, c, k, i, j;
    double mult, aux, num;
    x = malloc(sizeof(double)*n);
    val = alocaVariavel(n);

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
    printf("Matriz diagonalizada - Metodo de Jordan \n");
    imprimeMatriz(M, n, n+1);
    solucaoSistema(M, n, x);
}//metJordan

double calculaFuncao(int grau, double *equacao, double x){
    /*O método retorna o resultado de uma função em um ponto x, ou
    seja, seja f a função e x o ponto, o método calcula f(x).*/
    int i;

    double resp = 0;
    for(i=grau; i>=0; i--){
        resp += equacao[i]*pow(x, i);
    }
    return resp;
}//calculaFuncao()

void lagrange(int grau, double *coeficientes){
	double encontraK(int grau, double *coeficientes){
	    int i;
	    double aux = 0;
	    for(i = 0; i < grau ; i++){// Achar raiz superior
	        if(!(coeficientes[i] >= 0)){
	            if(aux < grau - i){
	                aux = grau - i;
	           }
	        }
	    }
	    return aux;
	}

	double encontraB(int grau, double *coeficientes){
		int i;
		double aux1 = 0;
		double modulo = -1;
		for(i = 0; i <= grau ; i++){
			if(!(coeficientes[i] >= 0)){
		    	if(aux1 > coeficientes[i]){
		        	aux1 = coeficientes[i];
		   		}
			}
		}
		return aux1*modulo;
	}

	void imprimeVetor(int n, double *v){
		int i;
		printf("[");
		for(i=n; i>=0; i--){
			printf("%4.3lf ", v[i]);
		}
		printf("]\n");
	}

	void inverteSinal(int n, double *v){
		int j;
		for(j=n; j>=0; j--){
			v[j] = - v[j];
		}
	}

	int i;
    double grau_double = grau, K[4], B[4], an[4], limites[4];
    double *coeficientes_inverso = (double*)malloc(sizeof(double)*(grau+1)); //equação com os indices invertidos
    double *coeficientes_ex_invertido = (double*)malloc(sizeof(double)*(grau+1));//equação com os expoentes impares invertidos
    double *coeficientes_ex_invertido_ex = (double*)malloc(sizeof(double)*(grau+1));//equação invertida com os expoentes impares invertidos

    //obtém a equação com os indices invertidos
    for(i = 0; i <= grau ; i++){
        coeficientes_inverso[i] = coeficientes[grau-i];
    }

    //obtém a equação com os expoentes impares invertidos
    for(i = 0; i <= grau ; i++){
        coeficientes_ex_invertido[i] = coeficientes[i];
        if(i%2 != 0 && i != 0){
            coeficientes_ex_invertido[i] = -coeficientes_ex_invertido[i];
        }
    }

    //obtém a equação invertida com os expoentes impares invertidos
    for(i = 0; i <= grau ; i++){
        coeficientes_ex_invertido_ex[i] = coeficientes_ex_invertido[grau-i];
    }

    if(coeficientes_inverso[0] < ERRO) inverteSinal(grau, coeficientes_inverso);

    if(coeficientes_ex_invertido[0] < ERRO) inverteSinal(grau, coeficientes_ex_invertido);

	if(coeficientes_ex_invertido_ex[0] < ERRO) inverteSinal(grau, coeficientes_ex_invertido_ex);

	imprimeVetor(grau, coeficientes);
	imprimeVetor(grau, coeficientes_inverso);
	imprimeVetor(grau, coeficientes_ex_invertido);
	imprimeVetor(grau, coeficientes_ex_invertido_ex);

    //calculando limite positivo
    // econtrando limite superior
    K[0] = encontraK(grau,coeficientes);
    B[0] = encontraB(grau,coeficientes);
    an[0] = coeficientes[0];
    limites[0] = 1.0 + pow(B[0]/an[0], 1.0/(grau_double-K[0]));

    // econtrando limite inferior
    K[1] = encontraK(grau,coeficientes_inverso);
    B[1] = encontraB(grau,coeficientes_inverso);
    an[1] = coeficientes_inverso[0];
    limites[1] = 1.0 / (1 + pow(B[1]/an[1], 1.0/(grau_double-K[1])));

    //calculando limite negativo
    // econtrando limite inferior
    K[2] = encontraK(grau,coeficientes_ex_invertido);
    B[2] = encontraB(grau,coeficientes_ex_invertido);
    an[2] = coeficientes_ex_invertido[0];
    limites[2] = -1.0*(1 + pow(B[2]/an[2], 1.0/(grau_double-K[2])));

    // econtrando limite superior
    K[3] = encontraK(grau,coeficientes_ex_invertido);
    B[3] = encontraB(grau,coeficientes_ex_invertido);
    an[3] = coeficientes_ex_invertido_ex[0];
    limites[3] = -1.0/(1 + pow(B[3]/an[3], 1.0/(grau_double-K[3])));

    imprimeVetor(4, limites);

    //Mostra a respota
    printf("\n\nLimite Positivo: %10.4lf %c%c X+ %c%c %10.4lf \n", limites[1], 60, 61, 60, 61, limites[0]);
    printf("Limite Negativo: %10.4lf %c%c X- %c%c %10.4lf \n", limites[2], 60, 61, 60, 61, limites[3]);

}

int bolzano(double a, double b, int grau, double *funcao){
    /*A função verifica se entre os intervalos dados, existe
    uma quantidade ímpar de raízes numa equação algébrica.
	Caso positivo, retorna 1. Se não, retorna 0.*/
    double Fa, Fb;
    Fa = calculaFuncao(grau, funcao, a);
    Fb = calculaFuncao(grau, funcao, b);
    if(Fa * Fb > 0){
        return 0;
    }
    else return 1;
}//bolzano()

void metBissecao(int grau, double *funcao){
    /*A função aplica o método da bisseção em uma equação algébrica
    dado um intervalo de observação.*/

    double a, b, m, erro, Fa, Fb, Fm = 1, i = 0;
    printf("Defina o intervalo: ");
    scanf("%lf%lf", &a, &b);
    fflush(stdin);
    erro = (b - a)/2;
    if(bolzano(a, b, grau, funcao)){
        printf("  a       b      m      f(a)     f(b)      f(m)      Erro    \n");
        printf("-----------------------------------------------------\n");
        while(i < 1000 && erro > ERRO && Fm != 0){
            m = (b + a)/2;
            erro = (b - a)/2;
            Fa = calculaFuncao(grau, funcao, a);
            Fb = calculaFuncao(grau, funcao, b);
            Fm = calculaFuncao(grau, funcao, m);
            printf("%10.8lf %10.8lf %10.8lf %10.8lf %10.8lf %10.8lf %10.8lf \n", a, b, m, Fa, Fb, Fm, erro);
            if(Fa * Fm > 0) a = m;
            else b = m;
            i++;
        }
        if (Fm == 0) printf("Raiz da equação = %lf\n", m);
        else printf("Raiz aproximada da equação = %7.8lf \n", m);
    }
    else printf("O numero de raízes no intervalo [%lf, %lf] e par\n", a, b);
}//metBissecao()

int main() {
    setlocale(LC_ALL, "Portuguese");
    int menu = 1;
    char op;
    double num;

    while(menu){
        printf("--------EP1--------\n");
        printf("C - Conversão\n");
        printf("S - Sistema Linear\n");
        printf("E - Equação Algébrica\n");
        printf("F - Finalizar\n");

        printf("-> Digite sua opção: ");
        scanf("%c", &op);
        fflush(stdin);

        switch(toupper(op)){
            case 'C':
                printf("\n***Conversão***\n");
                printf("-> Digite o número decimal a ser convertido: ");
                scanf("%lf", &num);
                fflush(stdin);
                conversao(num);
               	printf("\n");
                break;

            case 'S':
                printf("\n***Sistema Linear***\n");
                leMatriz(M, n);
                if(M == NULL){
                	break;
				}
                printf("Matriz: \n");
                imprimeMatriz(M, n, n+1);
                printf("\n");
                metJordan(M, n);
                printf("\n");
                break;

            case 'E':
                printf("\n***Equação Algébrica***\n");
                printf("Digite o grau do polinomio: ");
                scanf("%d", &grau);
                funcao = alocaFuncao(grau+1);
                if((funcao[0] > 0 & funcao[n] != 0) & funcao != NULL){
	                lagrange(grau, funcao);
	                metBissecao(grau, funcao);
	            }
	            else if(funcao == NULL){
	            	printf("Memória insuficiente!\n");
	            	break;
				}
	            else printf("A função deve ter o a%d maior que 0 e a0 diferente de 0", grau);
                printf("\n");
                break;

            case 'F':
                printf("Finalizando...\n");
                menu = 0;
                printf("\n");
                break;

            default:
                printf("Opção inválida!\n");
                printf("\n");
                break;
        }
    }
    return 0;
}//main()
