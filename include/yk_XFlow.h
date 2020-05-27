
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
#include <math.h>
#include "xf_Output.h"
#include "xf_Solver.h"
#include "xf_SolverUnsteady.h"
#include "xf_Arg.h"
/* include "xf_AdaptChaos_Common.h" */
#include "xf_State.h"
#include "xf_ErrEst.h"
#include "xf_Adapt.h"
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "struct_def.h"
#include <time.h>
#include "read_file.h"
#include "system_initialization.h"
#include "ks_initialConditions.h"
#include "ks_jbaraverage.h"
#include "injection.h"
#include "petsc.h"
#include "ks_errorEstimate.h"
#include "ks_errorFind.h"
#include "yk_createHROM.h"
#include "slepcsvd.h"
#include <mpi.h>
#include "ks_adjoint_Unsteady.h"
#include "yk_solverCFD.h"
#include "xf_SolverMultiStepStage.h"
#include "yk_leastSquaresShadowing.h"
#include "yk_calculateError.h"
#include "rdtsc.h"

//-----------------------------------------------------------------------------
// Structure used for storing for xflow specific fucntions and variables
//-----------------------------------------------------------------------------
typedef struct yk_Xflow{
  yk_PrimalSolver *ykflow;
  xf_KeyValue *KeyValueArg;
  xf_All *All;
  /* xf_UnsteadyData *UnsteadyData; */
  xf_VectorGroup *UG;
  xf_VectorGroup *UG_h;
  xf_VectorGroup *UG_H;
  xf_TimeHistData *TimeHistData;
  xf_TimeHistData *TimeHistData_h;
}yk_Xflow;

//-----------------------------------------------------------------------------
// Definitions defined here for function pointers
//-----------------------------------------------------------------------------
yk_PrimalSolver* new_Xflow();

void yk_findElementOfTarget(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, PetscInt *spatialSet,
			    Is_it *reduced);

void yk_findBoundarySeeds(yk_PrimalSolver *ykflow, Multiverse *multiquaiton,
			  Cluster *primal, int nSampNodes, int *numSeeds,
			  int *nodeSet);

void yk_fomBuffer(yk_PrimalSolver *ykflow, Universe equation, Is_it *reduced, int p,
                  char *nameOfDir);

void yk_xflow_totalResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, Vec residual, Mat Jacobian,
			    Is_it *reduced, int timeSolveNum);

void yk_xflow_function(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Cluster *primal, Vec vecObj, Is_it *reduced,
		       int timeNode);

void yk_xflow_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiqation,
			  Cluster *primal, Mat matObj, Is_it *reduced,
			  int timeNode);

void yk_findAdjacentElems(yk_PrimalSolver *ykflow, Cluster *primal, PetscInt **mshOs, PetscInt **temporalSet,
 			  Is_it *reduced);

void yk_uploadStateCalcSpatialAverageOutput(yk_PrimalSolver *ykflow,
					    Multiverse *multiquation,
                                            Galaxy *primal, double *aveSpace,
                                            int timeNode);
void yk_uploadStateCalcSpatialdJdU(yk_PrimalSolver *ykflow,
				   Multiverse *multiquation, Galaxy *primal,
                                   int timeNode, Vec vecObj, Is_it *reduced);

int yk_getOffDiagonalJacobians(yk_PrimalSolver *ykflow, Galaxy *primal,
			       Mat *Jacobian,
                                xf_Vector *U, xf_Vector *R,
				xf_JacobianMatrix *R_U,
				Is_it *reduced);

int yk_getOffDiagonalJacobiansGroup(yk_PrimalSolver *ykflow, Galaxy *primal,
				     Mat *Jacobian, Is_it *reduced);

void yk_MultInvMassMatrixArray(yk_PrimalSolver *ykflow, Galaxy * adjObj,
                               Vec adjoint, Is_it *reduced);

void yk_inject(yk_PrimalSolver *ykflow, Multiverse *multiquation,
	       Galaxy *primal_H, Galaxy * primal_h, Is_it *reduced);

void yk_anyArrayFunction(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Galaxy *primal, Vec vecObj, Is_it *reduced);

void delete_Xflow(yk_PrimalSolver *ykflow);
