#ifndef KS_DRES_DIFFUSION_H_
#define KS_DRES_DIFFUSION_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "u_dg.h"
#include "basisfunc.h"   
#include "petsc.h"
#include "tools.h"

void ks_dres_diffusion(Universe eqn, Galaxy *U, Mat dRdU, Is_it *reduced);

#endif
