
#ifndef KS_FUNCTION_H_
#define KS_FUNCTION_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "petsc.h"
#include "xiL.h"
#include "u_dg.h"
#include "ks_mass.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "tools.h"
#include "yk_solverCFD.h" 

void ks_function(yk_PrimalSolver *object, Multiverse *multiquation,
		 Cluster *primal, Is_it *reduced, Vec f, int timeNode);

void ks_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiquaiton,
		    Cluster *primal, Is_it *reduced, Mat dRdU, int timeNode);

#endif
