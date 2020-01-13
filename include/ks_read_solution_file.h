
#ifndef KS_READ_SOLUTION_FILE_H_
#define KS_READ_SOLUTION_FILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct_def.h"
#include "petsc.h"

void printSegSolution(Universe eqn, Galaxy *U, int segmentNum, int node);
void readSegSolution(Universe eqn, Galaxy *U, int segmentNum, int node);
void ks_printSolution(Universe eqn, Galaxy* U, int node);
void ks_readSolution(Universe eqn, Galaxy* U, int node);
void ks_removeSolution(Universe eqn, Galaxy *U, int node);
void ks_printReducedSolution(Universe eqn, Galaxy *U, int node);
void yk_readLSSInitialCon(char *filename, double R[]);
void yk_printLSSInitialCon(int size, double R[], int out_i);
void ks_readMatrix(char *filename, int row, Mat *sResidual);
void yk_fgets(FILE *stream, int row, double **snapshotArray, int *count);

#endif
