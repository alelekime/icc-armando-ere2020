#include "header.h"
#include "utils.h"

void alocaEDO(SL * sistemalinear) {
  sistemalinear -> matriz -> D = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
  sistemalinear -> matriz -> Di = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
  sistemalinear -> matriz -> Ds = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
  sistemalinear -> matriz -> B = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
  sistemalinear -> Y = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
  sistemalinear -> Xi = (double*)malloc(sistemalinear-> equacaoDiferencial -> n * sizeof(double));
}

void preencheEDO(SL *sistemalinear, int tipo) {
  switch (tipo) {
    case 1:
    sistemalinear -> equacaoDiferencial -> n = 5;
    sistemalinear -> equacaoDiferencial -> a = 0;
    sistemalinear -> equacaoDiferencial -> b = 12;
    sistemalinear -> equacaoDiferencial -> ya = 0;
    sistemalinear -> equacaoDiferencial -> yb = 0;
    break;
    case 2:
    sistemalinear -> equacaoDiferencial -> n = 10;
    sistemalinear -> equacaoDiferencial -> a = 0;
    sistemalinear -> equacaoDiferencial -> b = 12;
    sistemalinear -> equacaoDiferencial -> ya = 0;
    sistemalinear -> equacaoDiferencial -> yb = 0;
    break;
    case 3:
    sistemalinear -> equacaoDiferencial -> n = 5;
    sistemalinear -> equacaoDiferencial -> a = 0;
    sistemalinear -> equacaoDiferencial -> b = 1;
    sistemalinear -> equacaoDiferencial -> ya = 0;
    sistemalinear -> equacaoDiferencial -> yb = 1;
    break;
    case 4:
    sistemalinear -> equacaoDiferencial -> n = 10;
    sistemalinear -> equacaoDiferencial -> a = 0;
    sistemalinear -> equacaoDiferencial -> b = 1;
    sistemalinear -> equacaoDiferencial -> ya = 0;
    sistemalinear -> equacaoDiferencial -> yb = 1;
    break;
  }
}
double resolveEDO_R(double xi, int tipo) {
  double resultado;
  if (tipo == 1) {
    resultado = 6*xi - 0.5 * xi*xi;
  }else {
    resultado = 0;
  }
  return resultado;
}

double resolveEDO_Q(double xi, int tipo) {
  double resultado;
  if (tipo == 1) {
    resultado = 0;
  } else {
    resultado = 1;
  }
  return resultado;
}

double resolveEDO_P(double xi) {
  return 0;
}


void calculaDiagonais(SL* sistemalinear) {

  alocaEDO(sistemalinear);
  double a, b,c;
  double h = (sistemalinear-> equacaoDiferencial -> b - sistemalinear-> equacaoDiferencial -> a)/(sistemalinear-> equacaoDiferencial -> n +1);
  for (int i = 0; i < sistemalinear-> equacaoDiferencial -> n; i++) {

    if ( i == 0) {
      sistemalinear -> Xi[i] = sistemalinear -> equacaoDiferencial -> a + ( i + 1) * h;
      if (sistemalinear -> equacaoDiferencial -> b ==12) {
        a = resolveEDO_R(sistemalinear -> Xi[i],1);
        c = resolveEDO_Q(sistemalinear -> Xi[i],1);
      } else {
        a = resolveEDO_R(sistemalinear -> Xi[i],2);
        c = resolveEDO_Q(sistemalinear -> Xi[i],2);
      }
        b = resolveEDO_P(sistemalinear -> Xi[i]);

      sistemalinear -> matriz -> B[i] = h*h*a - sistemalinear -> equacaoDiferencial -> ya * (1 - (h*b)/2);  //q
      sistemalinear -> matriz -> D[i] = -2 + h*h*c;
      sistemalinear -> matriz -> Ds[i] = 1 + (h*b)/2; //p

    } else if (i == sistemalinear-> equacaoDiferencial -> n - 1) {
      sistemalinear -> Xi[i] = sistemalinear -> equacaoDiferencial -> b - h;
      if (sistemalinear -> equacaoDiferencial -> b ==12) {
        a = resolveEDO_R(sistemalinear -> Xi[i],1);
        c = resolveEDO_Q(sistemalinear -> Xi[i],1);
      } else {
        a = resolveEDO_R(sistemalinear -> Xi[i],2);
        c = resolveEDO_Q(sistemalinear -> Xi[i],2);
      }
        b = resolveEDO_P(sistemalinear -> Xi[i]);
      sistemalinear -> matriz -> B[i] = h*h*a - sistemalinear -> equacaoDiferencial -> yb * (1 + (h*b)/2);  //q
      sistemalinear -> matriz -> D[i] = -2 + h*h*c;
      sistemalinear -> matriz -> Di[i] = 1 - (h*b)/2; //p

    } else {
      sistemalinear -> Xi[i] = sistemalinear -> equacaoDiferencial -> a + ( i + 1) * h;
      if (sistemalinear -> equacaoDiferencial -> b ==12) {
        a = resolveEDO_R(sistemalinear -> Xi[i],1);
        c = resolveEDO_Q(sistemalinear -> Xi[i],1);
      } else {
        a = resolveEDO_R(sistemalinear -> Xi[i],2);
        c = resolveEDO_Q(sistemalinear -> Xi[i],2);
      }
        b = resolveEDO_P(sistemalinear -> Xi[i]);
      sistemalinear -> matriz -> B[i] =  h * h * a;
      sistemalinear -> matriz -> Di[i] = 1 - (h*b)/2;  //p
      sistemalinear -> matriz -> D[i] = -2 + h*h*c; //q
      sistemalinear -> matriz -> Ds[i] = 1 + (h*b)/2; //p

    }


  }
  sistemalinear-> equacaoDiferencial -> h = h;

}

void GaussSeidel(SL* sistemalinear) {
  double *interacaoAtual = (double*)malloc(sistemalinear -> equacaoDiferencial -> n * sizeof(double));
  double *interacaoAnterior = (double*)malloc(sistemalinear -> equacaoDiferencial -> n * sizeof(double));
  int i = 0;
  int interacoes = 0;
  for (int i = 0; i < sistemalinear -> equacaoDiferencial -> n; i++) {
    interacaoAtual[i] = 0;
    interacaoAnterior[i] = 0;
  }

  while (1) {
    int i = 0;
    interacaoAtual[i] = (sistemalinear -> matriz -> B[i] - sistemalinear -> matriz -> Ds[i] * interacaoAnterior[i + 1])/ sistemalinear -> matriz -> D[i];

    for (i = 0; i < sistemalinear -> equacaoDiferencial -> n - 1; i++) {
      interacaoAtual[i] = (sistemalinear -> matriz -> B[i] - (sistemalinear -> matriz -> Di[i] * interacaoAtual[i - 1])- (sistemalinear -> matriz -> Ds[i] * interacaoAnterior[i + 1]))/ sistemalinear -> matriz -> D[i];
    }
    interacaoAtual[i] = (sistemalinear -> matriz -> B[i] - sistemalinear -> matriz -> Di[i] * interacaoAtual[i -1])/ sistemalinear -> matriz -> D[i];
    interacoes++;
    if (erro(interacaoAtual, interacaoAnterior,sistemalinear -> equacaoDiferencial -> n) < 0.0001 || interacoes == 50) {
      printf("%f\n", erro(interacaoAtual, interacaoAnterior,sistemalinear -> equacaoDiferencial -> n) );
      for (int k = 0; k < sistemalinear -> equacaoDiferencial -> n; k++) {
        sistemalinear -> Y[k] = interacaoAtual[k];

      }
      break;
    } else {
      for (int k = 0; k < sistemalinear -> equacaoDiferencial -> n; k++) {
        interacaoAnterior[k] = interacaoAtual[k];
      }
    }

  }
}

void normaL2(SL* sistemalinear) {
{
int i,n;
double norma;
double *res;

n = SL->n;
res = (double *) malloc(sizeof(double)*n);

res[0] = SL->B[0] - SL->D[0]*SL->Y[0] - SL->Ds[0]*SL->Y[1];

for ( i = 1 ; i < n-1 ; i++ )
res[i] = SL->B[i] - SL->Di[i]*SL->Y[i-1] - SL->D[i]*SL->Y[i] - SL->Ds[i]*SL->Y[i+1];

res[n-1] = SL->B[n-1] - SL->Di[n-1]*SL->Y[n-2] - SL->D[n-1]*SL->Y[n-1];

norma = 0.0;
for ( i = 0 ; i < n ; i++ ) norma += res[i]*res[i];

free(res);
return sqrt(norma);

}

double erro(double *interacaoAtual,double *interacaoAnterior, int n){
  double maior;
  maior = interacaoAnterior[0] - interacaoAtual[0];
  for (int i = 1; i < n; i++){
    if (fabs(interacaoAnterior[i] - interacaoAtual[i]) > maior)
    maior = fabs(interacaoAnterior[i] - interacaoAtual[i]); // máximum das diferenças
    }
  return maior;
}
void imprimeEDO(SL* sistemalinear) {
    printf("***** item (a) Y\" = 6x - 0.5 x²: ");

    printf("n = %d, H = %f\n",sistemalinear-> equacaoDiferencial -> n,sistemalinear-> equacaoDiferencial -> h);

    printf("\nSistema Linear:\n\n");

    printf("%f %f ",sistemalinear-> matriz ->D[0],sistemalinear-> matriz -> Ds[0]);
    for (int j = 2 ; j < sistemalinear -> equacaoDiferencial -> n ; j++ ){
      printf("0 ");
    }
    printf("| %f\n\n",sistemalinear-> matriz -> B[0]);

    for (int i = 1 ; i < sistemalinear -> equacaoDiferencial -> n-1 ; i++ ){
      for (int j = 0 ; j < i-1 ; j++){
        printf("0 ");
      }

      printf("%f %f %f ",sistemalinear-> matriz -> Di[i],sistemalinear-> matriz -> D[i],sistemalinear-> matriz -> Ds[i]);
      for (int j = i+2 ; j < sistemalinear -> equacaoDiferencial -> n ; j++ ) {
        printf("0 ");
      }
      printf("| %f\n\n",sistemalinear-> matriz ->B[i]);
    }

    for (int j = 0 ; j < sistemalinear -> equacaoDiferencial -> n-2 ; j++ ) {
      printf("0 ");
    }
    printf("%f %f ",sistemalinear-> matriz ->Di[sistemalinear -> equacaoDiferencial -> n-1],sistemalinear-> matriz ->D[sistemalinear -> equacaoDiferencial -> n-1]);
    printf("| %f",sistemalinear-> matriz ->B[sistemalinear -> equacaoDiferencial -> n-1]);


    printf("\n");
    printf("\nY: ");
    for (int i = 0 ; i < sistemalinear -> equacaoDiferencial -> n ; i++ ) {
      printf("%9.5f ",sistemalinear-> Y[i]);
    }
    printf("\n");

    printf("Norma L2: , Tempo: %f ms\n\n",sistemalinear-> t);

}
