
#ifndef YK_CREATEHROM_H_
#define YK_CREATEHROM_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "yk_gaussNewtonSolve.h"
#include <time.h>
#include "system_initialization.h"
#include "qRFactorization.h"
#include "ks_offlineMatrix.h"
#include "ks_copyUtype.h"
#include "ks_greedyAlgorithm.h"
#include "yk_solverCFD.h"
#include "yk_performPOD.h"
#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_EqnSetHook.h"
#include "xf_Memory.h"
#include "xf_Residual.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Output.h"
#include "xf_Solver.h"
#include "xf_SolverUnsteady.h"
#include "xf_Arg.h"
/* #include "xf_AdaptChaos_Common.h"  */

void yk_precomputeCoefMat(Mat Z, Is_it *reduced);

void yk_createReducedOrderModel_ST(yk_PrimalSolver *ykflow,
                                   Multiverse *multiquation, Cluster *primal,
                                   Is_it *reduced);

void yk_runReducedOrderModel_ST(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Cluster *primalApprox,
			     Is_it *reduced);

void yk_destroyReducedOrderModel_ST(yk_PrimalSolver *ykflow,
                                    Multiverse *multiquation, Cluster *primal,
                                    Cluster *primalApprox, Is_it *reduced);
void yk_createReducedOrderModel(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Cluster *primalApprox, Is_it *reduced);

void yk_runReducedOrderModel(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Cluster *primalApprox,
			     Is_it *reduced);
void yk_createHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow,
                                        Multiverse *multiquation,
                                        Cluster *primal, Cluster *primalApprox,
                                        Is_it *reduced);
void yk_runHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Cluster *primalApprox,
			     Is_it *reduced);


void yk_destroyHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow,
					 Multiverse *multiquation, Cluster *primal,
					 Cluster *primalApprox, Is_it *reduced);

void yk_createHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				     Multiverse *multiquation, Cluster *primal,
				     Cluster *primalApprox, Is_it *reduced);

void yk_runHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				  Multiverse *multiquation, Cluster *primal,
				  Cluster *primalApprox,Is_it *reduced);

/* void minElements(Multiverse *multiquation, Cluster *primal, int *nodeSet, /\*  *\/ */
/*                  Mat *spaceIdentity, Is_it *reduced); */
void minElements(Multiverse *multiquation, Cluster *primal, int *nodeSet, int *timeSet,
                 Mat *spaceIdentity, Mat *ZZ, Is_it *reduced);
void findApproxMesh(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *primal, Cluster *primalApprox, Is_it *reduced);

void yk_findBoundarySeeds1D(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, int nSampNodes, int *numSeeds,
			    int *nodeSet);

void yk_findAdjacentElems1D(yk_PrimalSolver *ykflow, Cluster *primal, int elem,
                            int *elemSet);

void createInitialConditions_ST(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Mat *st_rOBTemp, Vec *estimatediniST,
				Is_it *reduced);
#endif
