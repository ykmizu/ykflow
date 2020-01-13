#include "struct_def.h"
//------------------------------------------------------------------------------
// Struc_def.c Contains all the tools for organization and array creation
//
// Yukiko Shimizu
// August 10, 2016
// Error Estimation for Chaotic Systems
//------------------------------------------------------------------------------

void initArray(Array *a, int size){
  int i;
  a->array = (double *)malloc(size * sizeof(double));
  for (i=0; i<size; i++){
    a->array[i]=0;
  }
  a->size = size;
}

void initArrayInt(intArray *a, int size){
  int i;
  a->array = (int *)malloc(size * sizeof(int));
  for (i=0; i<size; i++){
    a->array[i]=0;
  }
  a->count = size;
} 

void initArray2D(Array2D *a, int m, int n){
  int i, j;
  double **array2temp = (double**)malloc(m * sizeof(double*));
  for (i=0; i< m; i++){
    array2temp[i]=(double*)malloc(n*sizeof(double));
  }
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      array2temp[i][j]=0;
    }
  }
  a->array2=array2temp;
  a->xlength = m;
  a->ylength = n;
}

void initIs_it(Is_it *reduced){
  reduced->actualJbar = 0;
  reduced->delJ_psi = 0;
  reduced->delJ_dpsi = 0;
  reduced->delJ_dpsi_H1 = 0;
  reduced->tradDelJ = 0;
  reduced->dss = 0;
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  reduced->nTimeSegs = 0;
  reduced->restart = 0;
  reduced->reducedLSS = 0;
  reduced->innerStart = 0;
  reduced->nSampleNodes = 0;
  reduced->nBasisFuncs = 0;
  reduced->nBasisFuncsRJ = 0;
  reduced->nSnapshotsRJ = 0;
  reduced->outerStart = 0;
}

void initUtype(Galaxy *U){
  U->j_bar = 0;
  U->burnT_f = 0;
  U->space.dx = 0;
  U->space.elem.count = 1;
  U->space.node.count = 1;
  U->space.x_0 = 0;
  U->space.x_f = 0;
  U->basis.p = 0;
}

void delArray(Array *a){
  free(a->array);
  a->array = NULL;
}

void delArrayInt(intArray *a){
  free(a->array);
  a->array = NULL;
}

void delArray2D(Array2D *a){
  int i;
  for (i = 0; i < a->xlength; i++)
    free(a->array2[i]);
  free(a->array2);
  a->array2 = NULL;
}
