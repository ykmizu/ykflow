
#ifndef KS_ADJOINT_MATVEC_H_
#define KS_ADJOINT_MATVEC_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "struct_def.h"
#include "ks_mass.h"
#include "ks_read_solution_file.h"
#include "petsc.h"
#include "ks_copyUtype.h"
#include "system_initialization.h"
#include "tools.h"
#include "ks_ApplyTimeScheme.h"
#include "ks_function.h"
#include "ks_jbaraverage.h"
#include "yk_solverCFD.h"

void ks_adjoint_MATVEC(yk_PrimalSolver *ykflow, Multiverse * multiquation,
		       Cluster *primal, double *R, double *Rnew,
		       Is_it *Reduced);

#endif
