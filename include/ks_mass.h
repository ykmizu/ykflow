
#ifndef KS_MASS_H_
#define KS_MASS_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "basisfunc.h"
#include "petsc.h"
#include "petscsys.h" 
#include "tools.h"
void ks_mass(Universe equation, Galaxy* primal, Mat Mij, Is_it *reduced);

void ks_massBlock(Universe equation, Galaxy *primal, Mat *MijBlock);

void ks_massReduced(Universe equation, Galaxy* primal, Mat Mij, Is_it*reduced);

void yk_pMass(Universe equation, Galaxy *primal_h, Galaxy* primal_H, Mat *Mij);

void yk_project_I(Multiverse *multiquation, Galaxy *primal_h,              
                  Galaxy *primal_H, Mat *I_proj);
void yk_massH1(Universe equation, Galaxy *primal_c, Galaxy *primal_r, Mat *Mij);

void yk_projectH1_I(Multiverse *multiquation, Galaxy *primal_h,
		    Galaxy *primal_H, Mat *I_proj);

void yk_Identity(Mat *Imat, int size);
#endif
