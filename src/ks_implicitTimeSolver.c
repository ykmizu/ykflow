#include "ks_implicitTimeSolver.h"

void ks_implicitTimeSolve(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Cluster *primal, Is_it *reduced, int innert){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int newtonFlag = 0;
  int newtCounter = 0;
  int systemSize;
  int iterationCount = 0;
  int elementSize;
  PetscInt maxits = 1000;
  PetscScalar l2norm;
  PetscReal rtol = pow(10,-16);                                                
  PetscReal abstol = pow(10,-16);                                            
  PetscReal dtol= 1000;                                                      
  PetscInt iterNum;
  Vec residual = NULL;
  Vec objectVec = NULL;
  Mat dRdU;
  KSP ksp = NULL;
  PC pc;
  Universe _eqnOfInterest;
  int i;
  clock_t start, end;                                                          
  double cpu_time_used; 
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  int _eqnOfOriginal = multiquation->equation.numStates;
  if (reduced->reducedSolution == 1)
    _eqnOfInterest = multiquation->equationReduced;
  else if (reduced->reducedSolution == 0 &&
	   strncmp(primal->self->utypeId, "psi", 3)==0)
    _eqnOfInterest = multiquation->equationLSS;
  else
    _eqnOfInterest  = multiquation->equation;
  systemSize = primal->self->space.node.count*_eqnOfInterest.numStates;
  if (reduced->reducedSolution == 1)
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, systemSize, NULL, &dRdU);
  else if (reduced->reducedSolution == 0 &&
  	   strncmp(primal->self->utypeId, "psi", 3)==0)
    MatCreateSeqBAIJ(PETSC_COMM_SELF, _eqnOfOriginal*primal->state->basis.nodes
  		     ,systemSize, systemSize, 4, NULL, &dRdU); //check
  else
    MatCreateSeqBAIJ(PETSC_COMM_SELF, _eqnOfOriginal*primal->self->basis.nodes
  		     , systemSize, systemSize, 4, NULL, &dRdU);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &residual);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &objectVec);
  KSPCreate(PETSC_COMM_SELF, &ksp);
  while (newtonFlag == 0){
    start=clock();
    ks_totalResidual(ykflow, multiquation, primal, residual, reduced, innert);
    ks_totaldRdU(ykflow, multiquation, primal, dRdU, reduced, innert);    
    //-------------------------------------------------------------------------
    // Calculate the l2norm to check the progress of the newton iteration
    //-------------------------------------------------------------------------
    VecDot(residual, residual, &l2norm);
    newtCounter++; //Make sure it's converging and everything you  know
    if (newtCounter > 10){ //error management                               
      printf("NOT CONVERGING: IT'S TAKING A LONG TIME TO CONVERGE, YUKI!\n");
      exit(EXIT_FAILURE);         
    }
    iterationCount++;
    if (sqrt(l2norm) < pow(10, -9)) //Check 
      newtonFlag = 1;               //If so, the set the flag to break loop
    //-------------------------------------------------------------------------
    // Perform Newton Solve
    //-------------------------------------------------------------------------
    KSPSetOperators(ksp, dRdU, dRdU);
    KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);                         
    KSPSetType(ksp, KSPGMRES);                                    
    KSPSetFromOptions(ksp);    
    KSPGetPC(ksp, &pc); //Check signs of A and b.                              
    PCSetType(pc, PCBJACOBI); //KSPGetTolerance. <--- check this to see if
    VecScale(residual, -1); //Check to certain magnitude.                      
    KSPSolve(ksp, residual, residual);
    KSPGetIterationNumber(ksp, &iterNum);
    //-------------------------------------------------------------------------
    // Update the Current Solution                                          
    //-------------------------------------------------------------------------
    array2Vec(_eqnOfInterest, primal->self, objectVec);
    VecAXPY(objectVec, 1, residual);
    vec2Array(_eqnOfInterest, primal->self, objectVec);
    end=clock();
    cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC; 
    printf("n = %d, ε = %0.16e, t(s) =  %0.13f\n", iterationCount, sqrt(l2norm), cpu_time_used);
  }
  MatDestroy(&dRdU);
  VecDestroy(&residual);
  VecDestroy(&objectVec);
  KSPDestroy(&ksp);
}                                                                           
  
