
#ifndef KS_ERRORESTIMATE_H_
#define KS_ERRORESTIMATE_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "petsc.h" 
#include "injection.h"
#include "ks_GMRES.h"
#include "ks_adjoint_MATVEC.h"
#include "petsc.h"
#include <stdarg.h>
#include "yk_createHROM.h"
/* double ks_unsteadyError(Universe eqn, Utype* Ufine, int numT); */

double ks_leastSquaresShadow(Multiverse multiquation, Cluster *primal,
			     char *argv[],  Is_it *reduced);

#endif
