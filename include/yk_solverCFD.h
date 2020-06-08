
#ifndef YK_SOLVERCFD_H_
#define YK_SOLVERCFD_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "petsc.h"
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
#include "xf_SolverTools.h"
#include "xf_SolverUnsteady.h"
#include "xf_Arg.h"
/* #include "xf_AdaptChaos_Common.h"   */
#include "xf_State.h"
#include "xf_ErrEst.h"
#include "xf_Adapt.h"
#include "struct_def.h"
//#include "tools.h"
#include "xf_SolverMultiStepStage.h"
#include "xf_Basis.h"

typedef struct yk_PrimalSolver yk_PrimalSolver;

//-----------------------------------------------------------------------------
// Pointers to function
//-----------------------------------------------------------------------------
typedef void (*ptrBSeeds)(yk_PrimalSolver*, Multiverse *, Cluster *, int, int*,
			  int *);
typedef void (*ptrFunction)(yk_PrimalSolver *, Multiverse *, Cluster *, Vec,
			    Is_it *, int);
typedef void (*ptrdFunctiondu)(yk_PrimalSolver *, Multiverse *, Cluster *, Mat,
			       Is_it *,int);
typedef void (*ptrResidual)(yk_PrimalSolver*, Multiverse *, Cluster*, Vec,
			    Mat, Is_it *,int);
typedef void (*ptrFomName)(yk_PrimalSolver *, Universe, Is_it *, int, char *);

typedef void (*ptrAdjacentElem)(yk_PrimalSolver*, Cluster *, PetscInt **,
				PetscInt **,Is_it *);
typedef void (*ptrSpatialAveOutput)(yk_PrimalSolver *, Multiverse *, Galaxy *,
				    double *,int);
typedef void (*ptrSpatialdJdU)(yk_PrimalSolver *, Multiverse *, Galaxy *, int,
			       Vec, Is_it *);
typedef void (*ptrMultInvMass)(yk_PrimalSolver *, Galaxy *, Vec, Is_it*);
typedef void (*ptrInjection)(yk_PrimalSolver *, Multiverse *, Galaxy *,
			     Galaxy *, Is_it *);
typedef void (*ptrDelete)(yk_PrimalSolver *);

typedef void (*ptrAnyFunction)(yk_PrimalSolver *, Multiverse *, Galaxy *, Vec,
			       Is_it *);

typedef int (*ptrGetMassCoef)(yk_PrimalSolver *, Galaxy *, Mat *, Is_it *);

typedef void (*ptrTargetElem)(yk_PrimalSolver *, Multiverse *, Cluster *,
			      PetscInt *, Is_it *);

typedef void (*ptrResidualFile)(yk_PrimalSolver *, Multiverse *, Cluster *,
				Cluster *, Is_it*);

typedef struct yk_PrimalSolver{
  int numElemCol;
  void *solver;
  char path[xf_MAXSTRLEN];
  ptrFunction Function;
  ptrdFunctiondu dFunctiondu;
  ptrResidual Residual;
  ptrTargetElem minElements;
  ptrBSeeds boundarySeeds;
  ptrAdjacentElem adjacentElems;
  ptrSpatialAveOutput aveSpatialOutput;
  ptrSpatialdJdU spatialdJdU;
  ptrMultInvMass multiplyByInvMass;
  ptrInjection injectH2h;
  ptrDelete delete;
  ptrAnyFunction anyFunction;
  ptrFomName FomName;
  ptrGetMassCoef findMassCoefBDF;
  ptrResidualFile residualSnapshotFile;
}yk_PrimalSolver;

yk_PrimalSolver *new_Ykflow();

void delete_Ykflow(yk_PrimalSolver *ykflow);

void yk_ykflow_totalResidual(yk_PrimalSolver *primalSolverObj,
                             Multiverse *multiquation, Cluster *object,
			     Vec residual, Mat Jacobian, Is_it *reduced,
			     int timeSolveNum);

void yk_fomBuffer(yk_PrimalSolver *ykflow, Universe equation, Is_it *reduced,
		  int p, char *nameOfDir);

void yk_findBoundarySeeds1D(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, int nSampNodes, int *numSeeds,
			    int *nodeSet);

/* void yk_findAdjacentElems1D(yk_PrimalSolver *ykflow, Cluster *primal, int elem, */
/*                             int *elemSet); */

void yk_ykflow_function(yk_PrimalSolver *ykflow, Multiverse* multiquation,
                        Cluster *primal, Vec vecObj, Is_it *reduced,
                        int timeNode);

void yk_ykflow_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
                           Cluster *primal, Mat matObj, Is_it *reduced,
                           int timeNode);

void spatialOutputAverage(yk_PrimalSolver * ykflow, Multiverse *multiquation,
			  Galaxy *primal, double *value, int node);

void ks_dJdU(yk_PrimalSolver * ykflow, Multiverse *multiquation,
	     Galaxy *primal, int time_i, Vec dJdU, Is_it *reduced);

void inject(yk_PrimalSolver *ykflow, Multiverse *multiquation,
               Galaxy *primal_h, Galaxy * primal_H, Is_it *reduced);

void yk_any1D_function(yk_PrimalSolver *ykflow, Multiverse *multiquation,
                       Galaxy * primal, Vec vecObj, Is_it *reduced);

#endif
