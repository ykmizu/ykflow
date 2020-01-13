#ifndef _SYSTEM_INITIALIZATION_H_
#define _SYSTEM_INITIALIZATION_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "quadrature.h"
#include "basisfunc.h"
#include <string.h>
#include "ks_copyUtype.h"
//#include "ks_createHROM.h"
#include "yk_solverCFD.h"

#include <stdarg.h>  

void copyUniverse(Universe equation, Universe *equationCopy);
void createSystem(Universe eq, Galaxy *U);
//void createSystem(Universe equation, Galaxy *state, int numElems, int numNode);
void initializeStates(Universe eq, Galaxy *U);
void createMesh(Mesh *U);
void assignBasisQuad(Galaxy *U);
void destroySystem(Universe eq, Galaxy *U);
void destroyStates(Universe eq, Galaxy *U);
void destroyMesh(Mesh *space);
void deleteBasisQuad(Galaxy *U);
void copySystem(Universe equation, Galaxy *state, Galaxy *stateCopy);
void copyMesh(Mesh *space, Mesh *spaceCopy);
void createAdjointCluster(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Cluster *primal, Cluster *Psi, Is_it *reduced);
void destroyAdjointCluster(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			   Cluster *primal, Cluster *Psi, Is_it * reduced);
 
/* void createCluster(Galaxy *primal, Cluster *object); */
/* void findApproxSystem(Universe equation, Galaxy *state, Galaxy *stateApprox, */
/* 		      Is_it reduced); */


#endif
