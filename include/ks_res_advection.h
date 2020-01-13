

#ifndef KS_RES_ADVECTION_
#define KS_RES_ADVECTION_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "u_dg.h"
#include "basisfunc.h"
#include "tools.h"

void ks_res_advection(Universe eqn, Galaxy *U, Vec R, Is_it *reduced);

#endif
