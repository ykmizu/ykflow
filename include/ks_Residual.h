



#ifndef KS_RESIDUAL_
#define KS_RESIDUAL_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <stdarg.h>   
#include "struct_def.h" 
#include "ks_res_advection.h"
#include "ks_res_burgers.h"
#include "ks_res_diffusion.h"
#include "ks_res_fourth_IPDG.h"
#include "ks_res_tangent.h"
#include "ks_res_lagrange.h"
#include "petsc.h"
#include "yk_solverCFD.h"

void ks_totalResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster* object, Vec residual, Is_it *reduced,
		      int timeSolveNum);
void ks_spatialResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			Cluster *object, Vec residual, Is_it *reduced);
void ks_spatialStateResidual(Universe equation, Galaxy *object, Vec residual, 
			     Is_it *reduced);
void yk_ykflow_totalResidual(yk_PrimalSolver *primalSolverObj,
			     Multiverse *multiquation, Cluster *object,
			     Vec residual, Mat Jacobian, Is_it *reduced,
			     int timeSolveNum);
#endif
