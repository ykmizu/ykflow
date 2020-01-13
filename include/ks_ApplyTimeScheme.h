



#ifndef KS_APPLYTIMESCHEME_
#define KS_APPLYTIMESCHEME_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
//#include "ks_totalResidual.h"
//#include "ks_dRdU.h"
//#include "ks_copyUtype.h"
#include "ks_implicitTimeSolver.h"
#include "basisfunc.h"
#include "ks_mass.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "yk_solverCFD.h"

void ks_ApplyTimeScheme(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			Cluster *primal, Is_it *reduced);

#endif
