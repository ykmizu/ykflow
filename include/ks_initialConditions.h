

#ifndef KS_INITIALCONDITIONS_H_
#define KS_INITIALCONDITIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "ks_ApplyTimeScheme.h"
#include "time.h"
#include "rdtsc.h"

void ks_findAttractor(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *primal_H, Is_it *reduced);

#endif
