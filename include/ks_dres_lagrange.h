
#ifndef KS_DRES_LAGRANGE_H_
#define KS_DRES_LAGRANGE_H_

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
#include "yk_solverCFD.h"

void ks_dres_lagrange(yk_PrimalSolver *ykflow, Multiverse *multiquaiton,
		      Cluster *Psi1, Mat dRdU, Is_it *reduced);
#endif
