
#ifndef KS_RES_LAGRANGE_H_
#define KS_RES_LAGRANGE_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "petsc.h"
#include "ks_dRdU.h"
#include "petscvec.h"
#include "petscmat.h"
#include "u_dg.h"
#include "ks_read_solution_file.h"     
#include "ks_mass.h"
#include "ks_function.h"
#include "tools.h"
#include "basisfunc.h"
#include "ks_jbaraverage.h"

void ks_res_lagrange(yk_PrimalSolver *ykflow, Multiverse *multiqation,
		     Cluster *primal, Vec R, Is_it *reduced);
  
#endif
