#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#define ERRO 0.00000001

double *funcao(int n){
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

double calcular_funcao(int n, double v[], double x){
    int i;
    /* n = grau do polinomio
        v = vetor de coeficientes
        x = variavel da equa��o*/
    double resp = 0;
    for(i=n; i>=0; i--){
        resp += v[i]*pow(x, i);
    }
    return resp;
}

void met_bissecao(int n, double *v){
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
            Fa = calcular_funcao(n, v, a);
            Fb = calcular_funcao(n, v, b);
            Fm = calcular_funcao(n, v, m);
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

int bolzano(double a, double b, int n, double *v){
    double Fa, Fb;
    Fa = calcular_funcao(n, v, a);
    Fb = calcular_funcao(n, v, b);
    if(Fa * Fb > 0){
        return 0;
    }
    else return 1;
}

int main()
{
    setlocale(LC_ALL, "");
    int n, i, menu = 1;
    char op;
    double *v, x, resp;
    while(menu){
        printf("\n-----EP1-----\n");
        printf("C - Conversao\n");
        printf("S - Sistema Linear\n");
        printf("E - Equacao Algebrica\n");
        printf("F - Finalizar\n");

        printf("Digite sua opcao: ");
        scanf("%c", &op);
        fflush(stdin);
        switch(op){
            case 'C':
                printf("Conversao\n");
                break;
            case 'S':
                printf("Sistema Linear\n");
                break;
            case 'E':
                printf("Equação Algebrica\n");
                printf("Digite o grau do polinomio: ");
                scanf("%d", &n);
                v = funcao(n+1);
                /*printf("\nDigite o X da funcao que deseja calcular: ");
                scanf("%lf", &x);
                resp = calcular_funcao(n, v, x);
                printf("Valor de f(%lf) = %lf", x, resp); */
                met_bissecao(n, v);
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
