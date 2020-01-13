


#ifndef KS_RES_DIFFUSION_H_
#define KS_RES_DIFFUSION_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "u_dg.h"
#include "basisfunc.h"
#include "tools.h"

void ks_res_diffusion(Universe eqn, Galaxy *U, Vec R, Is_it *reduced);

#endif
