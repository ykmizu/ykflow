
#ifndef KS_DRES_TANGENT_H_
#define KS_DRES_TANGENT_H_

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
#include "tools.h"
#include "yk_solverCFD.h"

void ks_dres_tangent(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		     Cluster *Psi2, Mat dRdU, Is_it *reduced);

#endif
        
