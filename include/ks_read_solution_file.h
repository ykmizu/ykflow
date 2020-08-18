
#ifndef KS_READ_SOLUTION_FILE_H_
#define KS_READ_SOLUTION_FILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct_def.h"
#include "petsc.h"
#include "ks_greedyAlgorithm.h"
int ks_checkSolutionExist(Universe eqn, Galaxy *U, int node);
void printSegSolution(Universe eqn, Galaxy *U, int segmentNum, int node);
void readSegSolution(Universe eqn, Galaxy *U, int segmentNum, int node);
void ks_printSolution(Universe eqn, Galaxy* U, int node);
void ks_readSolution(Universe eqn, Galaxy* U, int node);
void ks_removeSolution(Universe eqn, Galaxy *U, int node);
void ks_printReducedSolution(Universe eqn, Galaxy *U, int node);
void yk_readLSSInitialCon(char *filename, double R[]);
void yk_printLSSInitialCon(int size, double R[], int out_i);
/* void ks_readMatrix(char *filename, int row, Mat *sResidual); */
Mat *ks_readMatrix(char *filename, PetscInt systemSize, PetscInt nTime, Mat *snap,
		   PetscInt *RJnum, Is_it *reduced);
void yk_fgets(FILE *stream, PetscInt row, double **snapshotArray, PetscInt *count);

#endif
