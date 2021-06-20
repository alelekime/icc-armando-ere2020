#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


#ifndef EQDIFERENCIAIS
#define EQDIFERENCIAIS

// Matriz tridiagonal
typedef struct {
  double *D, *Di, *Ds, *B;
  int n;
}SL_Tridiag;

typedef struct Edo{


  int n; // número de pontos internos na malha
  double a, b; // intervalo
  double ya, yb; // condições contorno
  double h;
}Edo;

typedef struct SL{

  Edo *equacaoDiferencial;
  SL_Tridiag *matriz;
  double *Y;
  double *Xi;
  double t;
  double *norma;

}SL;


void alocaEDO(SL * sistemalinear);
void preencheEDO(SL *sistemalinear, int tipo);
void calculaDiagonais(SL* sistemalinear);
double resolveEDO_Q(double xi, int tipo);
double resolveEDO_R(double xi, int tipo);
double resolveEDO_P(double xi);
double erro(double *interacaoAtual,double *interacaoAnterior, int n);
void GaussSeidel(SL* sistemalinear) ;
void imprimeEDO(SL* sistemalinear);
#endif
