
#include "ks_errorEstimate.h"

/* double ks_unsteadyError(Universe eqn, Utype* U, int numT){ */
/*   double delJ; */
/*   ks_adjoint_calcv2(eqn, U, numT); */
/*   delJ = ks_errorD(eqn, U, numT);     */
/*   return delJ; */
/* } */

double ks_leastSquaresShadow(Multiverse multiquation, Cluster *primal,
			     char *argv[], Is_it *reduced){
  int i;                                //initialization for iteration
  int solverSize = 0;
  double delJ = 0;
  double* R0;
  double* R;
  double* Rcheck;
  int systemSize;
  Universe _eqnOfInterest;
  Galaxy * _solOfInterest;
  Mat MijBlock;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    _eqnOfInterest = multiquation.equationReduced;
    _solOfInterest = primal->reduced;
  }else{
    _eqnOfInterest = multiquation.equation;
    _solOfInterest = primal->self;
  }
  systemSize = _eqnOfInterest.numStates*_solOfInterest->space.node.count;
  solverSize = systemSize*(2*reduced->nTimeSegs-1);
  R0 = (double *) malloc (solverSize*sizeof(double));
  R = (double *) malloc (solverSize*sizeof(double));
  Rcheck = (double *) malloc (solverSize*sizeof(double));
  //---------------------------------------------------------------------------
  // Make an initial guess for the GMRES solver
  //---------------------------------------------------------------------------
  if (reduced->restart)
    yk_readLSSInitialCon(argv[2], R0);
  else
    for (i=0; i<solverSize; i++) //Make initial guess for ubold for lss
      R0[i] = 0;
  //---------------------------------------------------------------------------
  // Mass Matrix 
  //---------------------------------------------------------------------------
  ks_massBlock(multiquation.equation, primal->self, &MijBlock);
  MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
  inverseMatLU(MijBlock, &primal->self->Mij);
  //---------------------------------------------------------------------------
  // Calculation of offline matrices
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    ks_offlineMatrixLSS(multiquation.equation, primal->self, reduced);
    VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs,
		 &primal->reduced->djdu);
    MatMultTranspose(reduced->rOBState, primal->stateFull->djdu,
		     primal->reduced->djdu);
  }
  //---------------------------------------------------------------------------
  // Perform GMRES to solve Ax = b
  //---------------------------------------------------------------------------
  ks_GMRES(multiquation, primal, R0, reduced);
  //---------------------------------------------------------------------------
  // Verify that the adjoint has been found now
  //---------------------------------------------------------------------------
  primal->self->beta = 1;
  ks_adjoint_MATVEC(multiquation, primal, R0, R, reduced);
  printf("----------------------------------------------------------------\n");
  printf("|                      GMRES Error Checking                    |\n");
  printf("----------------------------------------------------------------\n");
  for (i=0; i<solverSize; i++)
    printf("Linear System Solve ε = %0.16e\n", R[i]);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    MatDestroy(&reduced->A);
    MatDestroy(&reduced->B);
    VecDestroy(&primal->reduced->djdu);  
  }
  MatDestroy(&MijBlock);
  MatDestroy(&primal->self->Mij);
  free(R0);
  free(R);
  free(Rcheck);
  return 0;
}
