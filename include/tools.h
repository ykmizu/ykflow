/*
 ============================================================================
 Name        : LorenzError.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
*/

#ifndef TOOLS_H_
#define TOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include "struct_def.h"
#include <math.h>
#include "petsc.h"
//#include "ks_function.h"
//#include "ks_copyUtype.h"
#include "system_initialization.h"
//#include "ks_copyUtype.h"
#include <time.h>


void yk_kron(Mat A,  Vec B, double ***kron);
double dot(double* a, double*b, int size);
void projection(Universe eqn, Vec U, Vec vp);
void getSegArray(Universe eqn, Galaxy *U, double *R, int seg, Vec vp);
void data2Array(Universe eqn, Galaxy *U, double *u);
void array2Data(Universe eqn, Galaxy *U, double *u);
void array2Vec(Universe eqn, Galaxy *U, Vec v_petsc);
void vec2Array(Universe eqn, Galaxy *U, Vec vcon);
void ks_MatRow2Vec(Mat M, PetscInt m, PetscInt idxm[], Vec V, PetscInt rowi,
                   InsertMode iora);
void ks_MatCol2Vec(Mat M, PetscInt m, PetscInt idxm[], Vec V, PetscInt coli,
                   InsertMode iora);
void ks_Vec2MatCol(Mat M, PetscInt m, PetscInt idxm[], PetscInt col, Vec V,
                   InsertMode iora);
void ks_Vec2MatRow(Mat M, PetscInt m, PetscInt idcm[], PetscInt row, Vec V,
                   InsertMode iora);
int *createMeshMap(Mesh *original);
void solveEquations(Mat A, Vec b, Vec x);
int minInt(int A, int B);
void ks_Mat2MatBlock(Mat A, PetscInt m, PetscInt idxm[], PetscInt n,
                     PetscInt idxn[], Mat ABlock, InsertMode iora);
/* void findCognateState(Multiverse multiquation, Cluster *object, int node, */
/*                       Is_it *reduced); */
void MatMatMultBAIJ(Mat A, Mat B, Mat C);
void MatBlockMatMultDense(Mat A, Mat B, Mat C);
void MatBlockMatMultBAIJ(Mat A, Mat B, Mat C, Galaxy *primal, Is_it *reduced);
void inverseMatLU(Mat A, Mat *invA);
void MatBlockMult(Mat A, Vec x, Vec b);

#endif
