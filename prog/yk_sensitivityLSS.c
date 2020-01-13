/* 
Sensitivity LSS code 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct_def.h"
#include <time.h>
#include "read_file.h"
#include "system_initialization.h"                                             
#include "ks_initialConditions.h"                                              
#include "ks_jbaraverage.h"                                                    
#include "ks_sensitivity.h"                                                   
#include "injection.h"                                                         
#include "petsc.h"                                                             
#include "ks_errorEstimate.h"                                                  
#include "ks_errorFind.h"                                                     
#include "ks_createHROM.h"   
#include "slepcsvd.h"
#include "mpi.h"
#include "ks_adjoint_Unsteady.h"


int main(int argc, char *argv[]){
  int i;
  double actualJbar = 0;
  double estimateJbar = 0;
  int nTimeSegs;
  int *parameters = (int *) malloc (7*sizeof(int));
  Multiverse multiquation_C;
  Multiverse multiquation_c; //since we're doing sensitivity calculation;
  Cluster *primal_C = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primal_c = (Cluster *) malloc (sizeof(Cluster));
  Is_it *reduced = (Is_it *) malloc (sizeof(Is_it));
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  primal_C->self = (Galaxy *) malloc (sizeof(Galaxy));
  primal_c->self = (Galaxy *) malloc (sizeof(Galaxy));
  reduced->reducedSolution = 0;
  for (i=2; i<6; i++)
    parameters[i] = 0;
  //-------------------------------------------------------------------------- 
  // Find parameter values                                                     
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  multiquation_C.equation = read_eqn(argc, argv);
  multiquation_c.equation = read_eqn(argc, argv);
  read_input(argc, argv, multiquation_C.equation, primal_c->self,
	     primal_C->self, parameters);   
  multiquation_c.equation.c[multiquation_C.equation.paramIndex] =
    multiquation_C.equation.paramValue;
  nTimeSegs = parameters[0];                                                   
  reduced->dss = parameters[1];          
  reduced->typeLSS = 0; //Sensitivity calculations
  //---------------------------------------------------------------------------
  //PETSC 
  //---------------------------------------------------------------------------
  SlepcInitialize(&argc, &argv, NULL, NULL);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  //---------------------------------------------------------------------------
  // Set the original C solution with the first c value to compare to
  //---------------------------------------------------------------------------
  strcpy(primal_C->clusterId, "state");
  strcpy(primal_C->self->utypeId, "state");
  createSystem(multiquation_C.equation, primal_C->self);
  ks_copyUtype(primal_C->self, primal_c->self);
  sprintf(multiquation_c.equation.nameEqn, "%sPerturbed",
  	  multiquation_c.equation.nameEqn); 
  strcpy(primal_c->clusterId, "state");
  strcpy(primal_c->self->utypeId, "state");
  createSystem(multiquation_c.equation, primal_c->self);
  //---------------------------------------------------------------------------
  // Initial Conditions
  //---------------------------------------------------------------------------
  ks_findAttractor(multiquation_C, primal_C, reduced);
  ks_readSolution(multiquation_C.equation, primal_C->self, 0);
  ks_copySolutions(multiquation_c.equation, primal_C->self, primal_c->self);
  ks_printSolution(multiquation_c.equation, primal_c->self, 0);
  //---------------------------------------------------------------------------
  // Discontinuous Galerkin Solver
  //---------------------------------------------------------------------------
  printf("----------------------------------------------------------\n");
  printf("|                 s_0 State Solve                        |\n");
  printf("----------------------------------------------------------\n");
  ks_ApplyTimeScheme(multiquation_C, primal_C, reduced);
  ks_jbaraverage(multiquation_C.equation, primal_C->self, reduced);
  printf("----------------------------------------------------------\n");
  printf("|                 s_1 State Solve                        |\n");
  printf("----------------------------------------------------------\n");
  ks_ApplyTimeScheme(multiquation_c, primal_c, reduced);
  ks_jbaraverage(multiquation_c.equation, primal_c->self, reduced);
  actualJbar = primal_C->self->j_bar-primal_c->self->j_bar;
  //---------------------------------------------------------------------------
  // Least Squares Shadowing for Sensitivity Analysis to double check implme
  //---------------------------------------------------------------------------
  multiquation_c.equationLSS.numStates = primal_c->self->space.node.count;
  strcpy(multiquation_c.equationLSS.nameEqn, "meks");
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  ks_leastSquaresShadow(multiquation_c, primal_c, nTimeSegs, argv, 0, reduced);
  //---------------------------------------------------------------------------
  // Sensitivity calculation/ output error prediction
  //---------------------------------------------------------------------------
  estimateJbar = yk_sensitivityLSScalc(multiquation_c, primal_c, nTimeSegs);
  estimateJbar *=
    multiquation_C.equation.c[multiquation_C.equation.paramIndex]-
    multiquation_c.equation.c[multiquation_C.equation.paramIndex];
  //---------------------------------------------------------------------------
  // Results
  //---------------------------------------------------------------------------
  printf("Actual Error = %0.16e\n", actualJbar);
  printf("Estiamted Error = %0.16e\n", estimateJbar);
  //---------------------------------------------------------------------------
  // Destroy and Free Everything
  //---------------------------------------------------------------------------
  VecDestroy(&primal_C->self->djdu);
  VecDestroy(&primal_c->self->djdu);
  destroySystem(multiquation_C.equation, primal_C->self);
  destroySystem(multiquation_c.equation, primal_c->self);
  free(multiquation_C.equation.c);
  free(multiquation_c.equation.c);
  free(parameters);
  free(reduced);
  free(primal_C->self);
  free(primal_c->self);
  free(primal_C);
  free(primal_c);
  SlepcFinalize();
}
