

#ifndef YK_GAUSSNEWTONSOLVE_H_
#define YK_GAUSSNEWTONSOLVE_H_

#include <stdio.h>
#include <stdlib.h>
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
/* #include "xf_AdaptChaos_Common.h"                                               */
#include "xf_State.h"
#include "xf_ErrEst.h"
#include "xf_Adapt.h"
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "ks_copyUtype.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "ks_podv2.h"
#include "qRFactorization.h"
#include "yk_solverCFD.h"

void yk_gaussNewtonSolve_ST(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primalApprox, Is_it *reduced, int tSolve,
			    int *iterationCount);

void yk_gaussNewtonSolve(yk_PrimalSolver *solver, Multiverse *multiquation,
			 Cluster *primalApprox, Is_it *reduced, int tSolve,
			 int *iterationCount);

#endif
