


#ifndef KS_JBARAVERAGE_H_
#define KS_JBARAVERAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "u_dg.h"
#include "struct_def.h"
#include "ks_read_solution_file.h"
#include "basisfunc.h"
#include <complex.h>
#include "yk_solverCFD.h"

void ks_jbaraverage(yk_PrimalSolver *ykflow, Multiverse *multiquaiton,
		    Galaxy *primal, Is_it *reduced);


void spatialOutputAverage(yk_PrimalSolver * ykflow, Multiverse *multiquation,   
                          Galaxy *primal, double *value, int node); 

void ks_dJdU(yk_PrimalSolver *ykflow, Multiverse *multiquation, Galaxy *primal,
             int time_i, Vec dJdU, Is_it *reduced);

void ks_calculatedjdu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Galaxy *primal, Vec djdu, int timeNode, Is_it *reduced);
#endif
