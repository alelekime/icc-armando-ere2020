#include "header.h"
#include "utils.h"

int main(int argc, char const *argv[]) {
  double t;
  SL *sistemalinear = (SL*)malloc(sizeof(SL));
  sistemalinear -> equacaoDiferencial = (Edo*)malloc(sizeof(Edo));
  sistemalinear -> matriz = (SL_Tridiag*)malloc(sizeof(SL_Tridiag));


  preencheEDO(sistemalinear, 1);
  calculaDiagonais(sistemalinear);
  t = timestamp();
  GaussSeidel(sistemalinear);
  t = timestamp() - t;
  sistemalinear ->t = t;
  imprimeEDO(sistemalinear);
  printf("----------------------------------------\n" );
  t=0;
  preencheEDO(sistemalinear, 2);
  calculaDiagonais(sistemalinear);
  t = timestamp();
  GaussSeidel(sistemalinear);
  t = timestamp() - t;
  sistemalinear ->t = t;
  imprimeEDO(sistemalinear);
  printf("----------------------------------------\n" );
  t=0;

  preencheEDO(sistemalinear, 3);
  calculaDiagonais(sistemalinear);
  t = timestamp();
  GaussSeidel(sistemalinear);
  t = timestamp() - t;
  sistemalinear ->t = t;
  imprimeEDO(sistemalinear);
printf("----------------------------------------\n" );
t=0;

  preencheEDO(sistemalinear, 4);
  calculaDiagonais(sistemalinear);
  t = timestamp();
  GaussSeidel(sistemalinear);
  t = timestamp() - t;
  sistemalinear ->t = t;
  imprimeEDO(sistemalinear);
  printf("----------------------------------------\n" );

  free(sistemalinear);
  return 0;
}
