#ifndef KS_DRDU_H_
#define KS_DRDU_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <stdarg.h>
#include "struct_def.h"
#include "ks_Residual.h"
#include "petsc.h"
#include "getBDFCoef.h"
#include "ks_mass.h"
#include "ks_dres_advection.h"
#include "ks_dres_burgers.h"
#include "ks_dres_diffusion.h"
#include "ks_dres_fourth_IPDG.h"
#include "ks_dres_tangent.h"
#include "ks_dres_lagrange.h"
#include "ks_podv2.h"
#include "petscmat.h"


void ks_totaldRdU(yk_PrimalSolver *ykflow, Multiverse *multiquation,       
                  Cluster *object, Mat dRdU, Is_it *reduced, 
                  int innertimesolver);

void ks_spatialdRdU(yk_PrimalSolver *ykflow, Multiverse *multiquation,  
                    Cluster *object, Mat dRdU, Is_it *reduced);

  
void ks_spatialStatedRdU(Universe equation, Galaxy *object, Mat dRdU,
			 Is_it *reduced);

void ks_spatialdRdUCheck(Multiverse multiquation, Cluster *object,      
			 Is_it *reduced);
  
#endif
