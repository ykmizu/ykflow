

#ifndef INJECTION_H_
#define INJECTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct_def.h"
#include "u_dg.h"
#include "xiL.h"
#include "ks_read_solution_file.h"
#include "yk_solverCFD.h"

void inject(yk_PrimalSolver *ykflow, Multiverse *multiquation, Galaxy *primal,
            Galaxy *primalFine, Is_it *reduced);
  
//void injectAll(Universe eqn, Galaxy *U, Galaxy *Ufine); 
#endif 