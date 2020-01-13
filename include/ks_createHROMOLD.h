
#ifndef KS_CREATEHROM_H_
#define KS_CREATEHROM_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "ks_podv2.h"
#include "ks_gaussNewtonSolve.h"
#include <time.h>
#include "system_initialization.h"
#include "qRFactorization.h"
#include "ks_offlineMatrix.h"
#include "ks_copyUtype.h"
#include "ks_greedyAlgorithm.h"

void runReducedOrderModel(Multiverse multiquation, Cluster *primal,
			  Cluster *primalApprox, Is_it *reduced);
void runHyperReducedOrderModel(Multiverse multiquation, Cluster *primal,
			       Cluster *primalApprox, Is_it *reduced);
void createReducedOrderModel(Multiverse multiquation, Cluster *primal,
       			     Cluster *primalApprox, Is_it *reduced); 
void createHyperReducedOrderModel(Multiverse multiquation, Cluster *primal,   
                                  Cluster * primalApprox, Is_it *reduced); 
/* void createHyperReducedOrderModel(Universe equation, Galaxy *primal, int dss, */
/*                                   int numBasisFunctions, int ns, */
/*                                   int nResidualSnapshots, */
/*                                   int nJR, int *iterationCount, Mat Rc, Mat B); */
/* void reducedState2FullState(Universe equation, Galaxy *primal, */
/*                             Galaxy *primalReduced, Mat ROBState); */
/* void minElements(Galaxy *primal, PetscInt *nodeSet, Is_it *reduced); */
/* /\* void findReducedApproxMesh(Galaxy *primal, Mesh *reduced, int **tempM, int *sum); *\/ */
/* void findApproxMesh(Galaxy *primal, Galaxy *primalApprox, Is_it *reduced); */
#endif
