
#include "ks_dres_tangent.h"

void ks_dres_tangent(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		     Cluster *Psi2, Mat dRdU, Is_it *reduced){
  //---------------------------------------------------------------------------
  // Creaet the reduced priaml solution to be used to find the approx solut
  //---------------------------------------------------------------------------
  int timeNode = Psi2->self->time.node;
  //---------------------------------------------------------------------------
  // Calculation of the spatial residual                                 
  //---------------------------------------------------------------------------
  ks_dfunctiondu(ykflow, multiquation, Psi2, reduced, dRdU, timeNode);   
  MatScale(dRdU, -1);
  //---------------------------------------------------------------------------
  // Release and Destroy everything
  //---------------------------------------------------------------------------
}


/* void ks_dres_tangent(Universe equation, Utype *Psi2, Mat dRdUs, */
/*                      Mesh *reducedMesh, ...){ */
/*   int i,j; //initialization for iteration */
/*   int systemSize = Psi2->space.node.count*equation.numStates; */
/*   int nBasisFuncs = Psi2->space.node.count*equation.numStates;   */
/*   Utype *primal = (Utype*) malloc (sizeof(Utype)); //Iono have this here */
/*   int *numStatesIndex = (int *) malloc (equation.numStates*sizeof(int));  */
/*   va_list lls_list;  */
/*   Mat rOBState; */
/*   Mat B; */
/*   Vec uSolution; */
/*   Vec primalApproxVec; */
/*   Universe equationOriginal; */
/*   Utype *primalApprox; */
/*   Mesh *reducedMeshOriginal; */
/*   PetscInt rOBrow; */
/*   PetscInt rOBcol; */
/*   Mat dfdu; */
/*   Mat Zmult; */
/*   va_start(lls_list, reducedMesh); */
/*   equationOriginal = va_arg(lls_list, Universe); */
/*   primalApprox = va_arg(lls_list, Utype*); */
/*   reducedMeshOriginal = va_arg(lls_list, Mesh*);  */
/*   rOBState = va_arg(lls_list, Mat); */
/*   Zmult = va_arg(lls_list, Mat); */
/*   B = va_arg(lls_list, Mat); */
/*   va_end(lls_list); */
/*   MatGetSize(rOBState, &rOBrow, &rOBcol);   */
/*   for (i=0; i<equation.numStates; i++) */
/*     numStatesIndex[i] =i; */
/*   //---------------------------------------------------------------------------- */
/*   //Create the reduced primal solution to be used to find the approx solution */
/*   //---------------------------------------------------------------------------- */
/*   ks_copyUtype(equation, Psi2, primal); */
/*   strcpy(primal->utypeId, "state"); */
/*   primal->solution = (Array *) malloc (sizeof(Array)); */
/*   primal->space.node.count = nBasisFuncs; */
/*   initArray(primal->solution, nBasisFuncs);  */
/*   sprintf(primal->id, "%s_%s_%d_%d_hyper_reduced", equationOriginal.nameEqn, */
/*           primal->utypeId, primal->basis.p, primal->time_method);  */
/*   //---------------------------------------------------------------------------- */
/*   // Implementation */
/*   //---------------------------------------------------------------------------- */
/*   VecCreateSeq(PETSC_COMM_SELF, equation.numStates, &uSolution); */
/*   VecCreateSeq(PETSC_COMM_SELF, rOBrow, &primalApproxVec); */
/*   //---------------------------------------------------------------------------- */
/*   // Need to retrieve the U solution for the particular time node */
/*   //---------------------------------------------------------------------------- */
/*   ks_readSolution(equationOriginal, primal, Psi2->time.node); */
/*   VecSetValues(uSolution, equation.numStates, */
/*                numStatesIndex, primal->solution[0].array, INSERT_VALUES); */
/*   MatMult(rOBState, uSolution, primalApproxVec); */
/*   vec2Array(equationOriginal, primalApprox, primalApproxVec);   */
/*   //---------------------------------------------------------------------------- */
/*   // Find the dfdu */
/*   //---------------------------------------------------------------------------- */
/*   ks_dfunctiondu(equationOriginal, primalApprox, reducedMeshOriginal, */
/*                  rOBState, &dfdu, systemSize, Zmult, B); */
/*   MatScale(dfdu, -1); */
/*   MatCopy(dfdu, dRdUs, SAME_NONZERO_PATTERN); */



/*   VecDestroy(&uSolution); */
/*   VecDestroy(&primalApproxVec); */
/*   MatDestroy(&dfdu); */
/*   delArray(primal->solution); */
/*   free(primal->solution); */
/*   free(primal); */
/*   free(numStatesIndex); */
/* } */
  
