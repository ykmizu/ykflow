
#include "ks_initialConditions.h"

/* void setUtypeNil(Universe eqn, Galaxy *U){ */
/*   int i, j;   //initialization for iteration */
/*   for (i=0; i<eqn.numStates; i++){ */
/*     for (j=0; j<U->space.node.count; j++){ */
/*       U->solution[i].array[j] = 0; */
/*     } */
/*   } */
/* } */

void ks_findAttractor(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *primal_H, Is_it *reduced){
  int i, j;                 //initialization for iteration
  char removeFile[1000];
  primal_H->self->time.t_f = primal_H->self->burnT_f;
  double totalTime = primal_H->self->time.t_f-primal_H->self->time.t_0;
  double totalSpace = primal_H->self->space.x_f-primal_H->self->space.x_0;
  double halfValue = floor(totalSpace/2.0);
  double randomV;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  primal_H->self->time.count = totalTime/primal_H->self->time.dt; 
  //---------------------------------------------------------------------------
  // Implementation of the Initial Condition
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Random number generator
  //---------------------------------------------------------------------------
  srand(rdtsc());
  //randomV = 0;	      
  randomV = (double )rand()/(double)RAND_MAX*0.2-0.1;
  // Initial conditions assigned along with perturbed initial conditions
  for (i=0; i<multiquation->equation.numStates; i++)
    for (j=0; j<primal_H->self->space.node.count; j++)
      primal_H->self->solution[i].array[j] = 0 + randomV;
  primal_H->self->solution[0].array[1] = 1 + randomV; //x = 0.5
  ks_printSolution(multiquation->equation, primal_H->self, 0);
  //---------------------------------------------------------------------------
  // Run a burn Time and then print the new initial solution
  //---------------------------------------------------------------------------
  ks_ApplyTimeScheme(ykflow, multiquation, primal_H, reduced); //Sho
  for (i=0; i<primal_H->self->time.count+1; i++){
    sprintf(removeFile, "%s_%d.dat", primal_H->self->id, i);
    remove(removeFile);
  }
  // print the initial conditions for the next phase of the simulation
  ks_printSolution(multiquation->equation, primal_H->self, 0); //Print initial 
  //---------------------------------------------------------------------------
  // Need to return the value of count to its original value
  //---------------------------------------------------------------------------
  primal_H->self->time.t_f = primal_H->self->time.globalT_f;
  primal_H->self->time.count =
    (primal_H->self->time.globalT_f-
     primal_H->self->time.globalT_0)/primal_H->self->time.dt;
}
