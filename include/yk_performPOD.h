
#ifndef YK_PERFORMPOD_H_
#define YK_PERFORMPOD_H_

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
#include "yk_solverCFD.h"

Mat * yk_createSnapshotState(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, Mat *snapshot, int snapshotStep,
			    Is_it *reduced);

void yk_properOrthogonalDecompose(Mat *snapshot, int systemSize,
				  PetscInt *index,
				  PetscInt *numSingularValues,
				  PetscScalar engyValues,
				  Mat *A);

void yk_properOrthogonalDecomposeGroup(Mat *snapshot, int numSnapshots,
                                       int systemSize,
                                       PetscInt *index,
				       PetscInt *numSingularValues,
                                       PetscScalar engyValues,
				       Mat *A);

void yk_createTemporalSnapshotState(yk_PrimalSolver *ykflow,
                                    Multiverse *multiquaiton, Cluster *primal,
                                    Mat *snapshot, Mat *s_mu, int snapshotStep,
                                    Is_it *reduced);

void yk_createSpaceTimeBasis(Mat *spacebasis, Mat *timebasis, Mat *spaceTime);

void yk_formatSpaceTimeBasis(Cluster *primal, Mat *subspaceTime,
                             Mat *spaceTime, Mat *spaceTime_i, Is_it *reduced);
/* void yk_createUltimateSpaceTimeBasis(yk_PrimalSolver *ykflow, */
/*                                      Multiverse *multiquation, Cluster *primal, */
/*                                      Cluster *primalApproxm, Is_it *reduced); */
#endif
