#include <stdio.h>
#include <stdlib.h>

double **alocaMatriz(int linhas, int colunas);
void leMatriz(double **matriz, int linhas, int colunas);
void imprimeMatriz(double **matriz, int linhas, int colunas);
void Gauss (double **M, int n);

int main(){
	
	int i, j, n;
	double **matriz;
	printf("qntdd de variaveis: ");
	scanf("%d", &n);
	
	matriz = alocaMatriz(n, n+1);
	
	if(matriz == NULL){
		printf("Deu pal");
		return 0;
	}

	leMatriz(matriz,n, n+1);
	imprimeMatriz(matriz ,n, n+1);
	
	Gauss(matriz, n);
	imprimeMatriz(matriz, n, n+1);
	
	return 0;
}

double **alocaMatriz(int linhas, int colunas){

	int i, j;
	double **matriz;
	matriz = malloc(linhas * sizeof(matriz));
	if(matriz == NULL){
		return 0;
	}
	for(i = 0; i < linhas; i++){
		matriz[i] = malloc(sizeof(matriz) * linhas);
		if(matriz[i] == NULL){
			for(j = 0; j < i; j++){
				free(matriz[j]);
			}
			free(matriz);
			return 0;
		}		
	}
	return matriz;
}

void leMatriz(double **matriz, int linhas, int colunas){
		
		int i, j;
		for(i = 0; i < linhas; i++){
			for(j = 0; j < colunas; j++){
				printf("matriz[%d][%d]: ", i+1, j+1);
				scanf("%lf", &matriz[i][j]);
			}
		}
		printf("\n");
	}
	
void imprimeMatriz(double **matriz, int linhas, int colunas){
		int i, j;
		for(i = 0; i < linhas; i++){
			for(j = 0; j < colunas; j++){
				printf("%.1lf ", matriz[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
void Gauss (double **M, int n){
	
	int i, j, k;
	double mult, *aux;
	
	for(i = 0; i < n; i++){
		if(M[i][i] == 0){
			j = i+1;
			while(j < n && M[j][i] == 0){
				j++;
			}
			if(j < n){
				aux = M[i];
				M[i] = M[j];
				M[j] = aux;
			}
		}
		if(M[i][i] != 0){
			for(j = i+1; j < n; j++){
				mult = -M[j][i]/M[i][i];
				M[j][i] = 0;
				for(k = i+1; k <= n; k++){
					M[j][k] = M[j][k] + mult*M[i][k];
				}
			}
		}
	}
}