
#include "ks_sensitivity.h"

void ks_sensitivity(Multiverse multiquation, Cluster *primal){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j;                             //Initialization for iterations
  int systemSize = primal->self->space.node.count*
    multiquation.equation.numStates;
  int *index = (int *) malloc (systemSize*sizeof(int));
  double derrorEst;                     //Initialize to 0: sum value
  char outfile[10000];
  PetscScalar val = .01;               //Perturbed the coarse solution
  Vec Uoutput;                          //del Residual Value (dRdU * Delta U)
  Vec dJdU;
  Vec psi1data;
  Vec perturb;
  Mat dRdU;                             //Jacobian Matrix at the initial time  
  Galaxy* Psi = (Galaxy *) malloc (sizeof(Galaxy));  //Adjoint Solution
  Cluster* primalU = (Cluster *) malloc (sizeof(Cluster));
  Cluster* primalDelU = (Cluster *) malloc (sizeof(Cluster)); //Perturbed S
  Is_it *reduced = (Is_it*) malloc (sizeof(Is_it));
  primalU->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalDelU->self = (Galaxy *) malloc (sizeof(Galaxy));
  //---------------------------------------------------------------------------
  // Initialization and Mallocing 
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &perturb); 
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &Uoutput);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &dJdU);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &psi1data);
  MatCreateSeqBAIJ(PETSC_COMM_SELF, primal->self->basis.p+1, systemSize,
		   systemSize, 3, NULL, &dRdU);
  for (i=0; i<systemSize; i++)
    index[i] = i;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Non-perturbed solution
  ks_copyUtype(primal->self, primalU->self);
  primalU->self->time.t_f = primal->self->time.dt*5;
  primalU->self->time.globalT_f = primal->self->time.dt*5;
  sprintf(primalU->clusterId, "state");
  sprintf(primalU->self->utypeId, "state");
  createSystem(multiquation.equation, primalU->self);
  sprintf(primalU->self->id, "%s_sensitivity", primalU->self->id);
  //Perturbed Solution
  ks_copyUtype(primal->self, primalDelU->self);
  primalDelU->self->time.t_f = primal->self->time.dt*5;
  primalDelU->self->time.globalT_f = primal->self->time.dt*5;
  sprintf(primalDelU->clusterId, "state");
  sprintf(primalDelU->self->utypeId, "state");
  createSystem(multiquation.equation, primalDelU->self);
  sprintf(primalDelU->self->id, "%s_perturbed", primalDelU->self->id);  //spec
  //---------------------------------------------------------------------------
  // Perturbed adjoint Solution
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, Psi); //Can do this since interpolation o
  sprintf(Psi->utypeId, "psi_unsteady");   //only for one segment
  createSystem(multiquation.equation, Psi);   //Create the system which 
  //---------------------------------------------------------------------------
  // Non-perturbed solution
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation.equation, primal->self, 0);
  ks_copySolutions(multiquation.equation, primal->self, primalU->self);
  ks_printSolution(multiquation.equation, primalU->self, 0);
  ks_ApplyTimeScheme(multiquation, primalU, reduced);
  ks_jbaraverage(multiquation.equation, primalU->self, reduced);
  //---------------------------------------------------------------------------
  // Perturbed injected solution
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation.equation, primal->self, 0);     
  ks_copySolutions(multiquation.equation, primal->self, primalDelU->self);
  VecSet(perturb, val);
  for (i=0; i<multiquation.equation.numStates; i++)
    for (j=0; j<primal->self->space.node.count; j++)
      primalDelU->self->solution[i].array[j] += val;
  ks_printSolution(multiquation.equation, primalDelU->self, 0);
  ks_ApplyTimeScheme(multiquation, primalDelU, reduced);
  ks_jbaraverage(multiquation.equation, primalDelU->self, reduced);
  //---------------------------------------------------------------------------
  // Tradtional adjoint calculation
  //---------------------------------------------------------------------------
  ks_adjoint_Unsteady(multiquation, primalDelU);  
  //---------------------------------------------------------------------------
  // Need to find dRdU
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation.equation, primalU->self, 0);
  ks_totaldRdU(multiquation, primalU, dRdU, reduced, 1);
  MatMult(dRdU, perturb, Uoutput); 
  //---------------------------------------------------------------------------
  // Restrieve the first adjoint values for error estimation
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation.equation, Psi, 0);
  array2Vec(multiquation.equation, Psi, psi1data);
  //---------------------------------------------------------------------------
  // error = Adjoint* \delR
  //---------------------------------------------------------------------------
  VecDot(psi1data, Uoutput, &derrorEst);
  printf("----------------------------------------------------------\n");
  printf("|    Sensitivity test results for traditional adjoints   |\n");
  printf("----------------------------------------------------------\n");
  printf("Actual %.16f\n", primalDelU->self->j_bar-primalU->self->j_bar);
  printf("Estimated %.16f\n", -derrorEst);
  //---------------------------------------------------------------------------
  // Delete the Ufine files and the Adjoint files /Don't need them
  //---------------------------------------------------------------------------
  for (i=0; i<primal->self->time.count+1; i++){
    ks_removeSolution(multiquation.equation, primalU->self, i);
    ks_removeSolution(multiquation.equation, primalDelU->self, i);
    ks_removeSolution(multiquation.equation, Psi, i);
  }
  sprintf(outfile, "%s_jbar.dat", primalDelU->self->id);  //special cas  
  remove(outfile);
  sprintf(outfile, "%s_jbar.dat", primalU->self->id);
  remove(outfile);
  //---------------------------------------------------------------------------
  // Free eveerything
  //---------------------------------------------------------------------------
  VecDestroy(&perturb);
  VecDestroy(&dJdU);
  VecDestroy(&psi1data);
  VecDestroy(&Uoutput);
  VecDestroy(&primalDelU->self->djdu);
  VecDestroy(&primalU->self->djdu);
  MatDestroy(&dRdU);
  destroySystem(multiquation.equation, primalU->self);
  destroySystem(multiquation.equation, primalDelU->self);
  destroySystem(multiquation.equation, Psi);
  free(reduced);
  free(primalDelU->self);
  free(primalDelU);
  free(primalU->self);
  free(primalU);
  free(Psi);
  free(index);
}


