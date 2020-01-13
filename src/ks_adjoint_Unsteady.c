
#include "ks_adjoint_Unsteady.h"

void ks_adjoint_Unsteady(Multiverse multiquation, Cluster *primal){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j; //initialization for iteration
  int innertimesolver = primal->self->time_method;
  int iterNum;
  int numPreviousTime;
  int systemSize = multiquation.equation.numStates*
    primal->self->space.node.count;
  Is_it *reduced = (Is_it *) malloc (sizeof(Is_it));  
  Vec djdu;
  Vec djdu_time;
  Vec adjoint;
  Vec tempAdjoint = NULL;
  Mat Mij;
  Mat dRdU;
  KSP ksp;
  PC pc;
  PetscInt maxits = 1000;                                         
  PetscReal rtol = pow(10,-12);                                    
  PetscReal abstol = pow(10,-16);                          
  PetscReal dtol= 1000;   
  Galaxy *Psi = (Galaxy *) malloc (sizeof(Galaxy));
  double **c = (double **)malloc((primal->self->time_method)*sizeof(double *));
  //  c[0] = (double *) malloc ((primal->self->time_method+1)*
  //			    primal->self->time_method*sizeof(double));
  for (i=0; i<primal->self->time_method; i++){
    c[i] = getBDFCoef(i+1);
  }
  //---------------------------------------------------------------------------
  // Initialization and Mallocing
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  MatCreateSeqAIJ(PETSC_COMM_SELF, systemSize, systemSize,           
                  primal->self->basis.p+1, PETSC_NULL, &Mij); 
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &djdu);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &djdu_time);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &tempAdjoint);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &adjoint);
  KSPCreate(PETSC_COMM_SELF, &ksp);  
  //---------------------------------------------------------------------------
  // Set up dJ/dU, you already have this info but put into Vec form
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, Psi);
  sprintf(Psi->utypeId, "psi_unsteady");
  createSystem(multiquation.equation, Psi);
  for (i = primal->self->time.count+1; i-->0;){
    MatCreateSeqBAIJ(PETSC_COMM_SELF, primal->self->basis.p+1, systemSize,  
                     systemSize, 3, NULL, &dRdU);   
    //-------------------------------------------------------------------------
    // Calculate the Jacobian dR/dU for current time i (Left hand side A)
    //-------------------------------------------------------------------------
    ks_readSolution(multiquation.equation, primal->self, i);
    ks_totaldRdU(multiquation, primal, dRdU, reduced, innertimesolver);
    //-------------------------------------------------------------------------
    //Now find the right hand side of this adjoint calculation (right hand B)
    //-------------------------------------------------------------------------
    ks_dJdU(multiquation.equation, primal->self, djdu_time, i);
    if (i < primal->self->time.count-primal->self->time_method)
      numPreviousTime = primal->self->time_method;
    else
      numPreviousTime = primal->self->time.count - i;
    // Extract the adjoint solution found previously (backwards substitution)
    for (j=1; j<numPreviousTime+1; j++){
      ks_mass(multiquation.equation, primal->self, Mij, reduced); //Don't
      ks_readSolution(multiquation.equation, Psi, i+j); 
      array2Vec(multiquation.equation, Psi, adjoint);
      if (i <primal->self->time_method -1 && j<primal->self->time_method-i){
	MatScale(Mij,  c[i+j-1][j]/primal->self->time.dt);
      }else{
	MatScale(Mij, c[primal->self->time_method-1][j]/primal->self->time.dt);
      }
      MatMult(Mij, adjoint, djdu);
      VecAXPY(djdu_time, 1, djdu);
    }
    // Copy djdu value so that you don't lose it when doing KSPSOLVE
    VecScale(djdu_time, -1);
    //-------------------------------------------------------------------------
    // Perform Ax = b to find the unsteady traditional adjoint values
    //-------------------------------------------------------------------------
    MatTranspose(dRdU, MAT_INPLACE_MATRIX, &dRdU);   
    KSPSetOperators(ksp, dRdU, dRdU);                                
    KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);              
    KSPSetType(ksp, KSPGMRES);                                
    KSPSetFromOptions(ksp);                                          
    KSPGetPC(ksp, &pc); //Check signs of A and b.                  
    PCSetType(pc, PCBJACOBI); //KSPGetTolerance. <--- check this to see if 
    KSPSolve(ksp, djdu_time, djdu_time);   
    KSPGetIterationNumber(ksp, &iterNum);
    VecCopy(djdu_time, adjoint);
    vec2Array(multiquation.equation, Psi, adjoint);
    ks_printSolution(multiquation.equation, Psi, i);
    //-------------------------------------------------------------------------
    // Determine which time solver to use based on which element in 
    //-------------------------------------------------------------------------
    if (i == 2) //Use BDF1 at thbeginn
      innertimesolver = 1;
    else if (i == 3 && primal->self->time_method == 3)
      innertimesolver = 2;
         
    MatDestroy(&dRdU);
  }
  //---------------------------------------------------------------------------
  // Destroy and Clear out Memory
  //---------------------------------------------------------------------------
  destroySystem(multiquation.equation, Psi);
  KSPDestroy(&ksp); 
  MatDestroy(&Mij);
  VecDestroy(&djdu);
  VecDestroy(&djdu_time);
  VecDestroy(&tempAdjoint);
  VecDestroy(&adjoint);
  free(reduced);
  free(Psi);
  for (i=0; i<primal->self->time_method; i++)
    free(c[i]);
  free(c);
}
