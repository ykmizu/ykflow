
#ifndef KS_OFFLINEMATRIX_H_
#define KS_OFFLINEMATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "petsc.h"
#include "ks_podv2.h"
#include "qRFactorization.h" 
#include <petscksp.h>

void ks_offlineMatrix(Universe equation, Galaxy *primal, Galaxy *primalApprox,
                       Is_it *reduced);
  
void ks_offlineMatrixLSS(Universe equation, Galaxy *object, Is_it *reduced);

void restrictedIdentity(Universe equation, Galaxy *object,Galaxy *objectApprox,
			PetscInt *nodeSet, Is_it *reduced);
  
void restrictedMatrices(Universe equation, Galaxy *object,                     
                        Galaxy *objectApprox, Mat *rOBResidualHat,
			Mat *rOBJacobianHat,                
                        PetscInt *nodeSet, Is_it *reduced);
#endif
