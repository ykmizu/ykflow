
#ifndef KS_GAUSSNEWTONSOLVE_H_
#define KS_GAUSSNEWTONSOLVE_H_

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
#include "ks_podv2.h"
#include "qRFactorization.h"


void ks_gaussNewtonSolve(Multiverse multiquation, Cluster *primalReduced,      
                         Is_it *reduced, int innertimesolver,
			 int *iterationCount);

#endif

  
