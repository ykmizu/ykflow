
#ifndef KS_ERRORFIND_H_
#define KS_ERRORFIND_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "ks_Residual.h"
#include "petsc.h"
#include "tools.h"
#include "ks_ApplyTimeScheme.h"
#include "ks_read_solution_file.h"
#include "injection.h"
#include "ks_jbaraverage.h"
#include "ks_adjoint_Unsteady.h"

double yk_sensitivityLSScalc(Multiverse multiquation, Cluster *primal, 
                         int numTimeSegs);
void ks_errorD(Multiverse multiquation, Cluster *primal_h,
		 Cluster *primal_H, Is_it *reduced);

void yk_tradErrorEst(Multiverse multiquation, Cluster *primal_h, 
				    Cluster *primal_H, Is_it *reduced);
#endif
