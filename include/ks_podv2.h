
#ifndef KS_PODV2_H_
#define KS_PODV2_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "ks_copyUtype.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "slepcsvd.h"  

void createSnapshotState(Universe equation, Galaxy *primal,
                         Mat *snapshot, int dss);
/* void approxState(Universe equation, Galaxy *primal, int dss, int snapshot_i, */
/*                  Mat reducedOrderBasis); */
void moorePenrosePseudoInv(Mat A, int rowSize, int colSize, Mat *Aplus);
void properOrthogonalDecompose(Mat snapshot, int systemSize,
                               PetscInt *numSingularValues, Mat *A);
/* void reduceMatRows(Universe equation, Galaxy *primal, Mat A, int column, */
/*                    Mat *Ahat, Mesh *altered); */
#endif
