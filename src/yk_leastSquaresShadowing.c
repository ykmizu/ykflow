
#include "yk_leastSquaresShadowing.h"

double yk_leastSquaresShadowing(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				char *argv[], Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;                                // initialization for iteration
  int solverSize;                       // size of the R vector in LSS
  int systemSize;                       // size of the system (FOM, ROM, HROM)
  double *R0;                           // initial guess for LSS adjoint sys
  double *R;                            // final solution for LSS adjoint sys
  Universe _eqnOfInterest;              // reduced or full order solution eqn
  Galaxy *_solOfInterest;               // differentiate from FOM and HROM red
  Mat MijBlock;                         // Mass Block Matrix
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    _eqnOfInterest = multiquation->equationReduced;
    _solOfInterest = primal->reduced;
  }else{
    _eqnOfInterest = multiquation->equation;
    _solOfInterest = primal->self;
  }
  systemSize = _eqnOfInterest.numStates*_solOfInterest->space.node.count;
  solverSize = systemSize*(2*reduced->nTimeSegs-1);
  R0 = (double *) malloc (solverSize*sizeof(double));
  R = (double *) malloc (solverSize*sizeof(double));
  //---------------------------------------------------------------------------
  // Make an initial guess for the GMRES solver
  //---------------------------------------------------------------------------
  if (reduced->restart)
    yk_readLSSInitialCon(argv[2], R0);
  else
    for (i=0; i<solverSize; i++)
      R0[i] = 0;
  //---------------------------------------------------------------------------
  // Mass Matrix
  //---------------------------------------------------------------------------
  if (ykflow->numElemCol == 3){
    ks_massBlock(multiquation->equation, primal->self, &MijBlock);
    MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
    inverseMatLU(MijBlock, &primal->self->Mij);
  }
  //---------------------------------------------------------------------------
  // Calculation of Offline Matrices
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)
    ks_offlineMatrixLSS(multiquation->equation, primal->self, reduced);
  ks_GMRES(ykflow, multiquation, primal, R0, reduced);
  //---------------------------------------------------------------------------
  // Now find the actual adjoints
  //---------------------------------------------------------------------------
  primal->self->beta = 1;
  ks_adjoint_MATVEC(ykflow, multiquation, primal, R0, R, reduced);
  for (i=0; i<solverSize; i++)
    printf("%g\n", R[i]);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  if (ykflow->numElemCol == 3){
    MatDestroy(&MijBlock);
    MatDestroy(&primal->self->Mij);
  }
  free(R0);
  free(R);
  return 0;
}
