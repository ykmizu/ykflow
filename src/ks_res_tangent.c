#include "ks_res_tangent.h"
//-----------------------------------------------------------------------------
// Calculates the residual R  for tangent equation forleastsquaresshadowing
//                                                                         
// -   dv   ∂f                                                             
// R = -- - -- v = 0  I think dfdu is a constant scalar term...            
//     dt   ∂u                                                             
//                                                                         
// Yukiko Shimizu                                                          
// June 21, 2016                                                  
// Error Estimation for Choatic Systems                           
//-----------------------------------------------------------------------------

void ks_res_tangent(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster* Psi2, Vec residual, Is_it *reduced){
  int systemSize;          //total nodes in this particular system
  int elementSize;
  Vec v_petsc;
  Mat dfdu = NULL;
  Vec dfds;
  Universe _eqnOfInterest;
  //---------------------------------------------------------------------------
  // Mallocing things and crafting things and initializing things
  //---------------------------------------------------------------------------
  elementSize = multiquation->equation.numStates*
    Psi2->state->basis.nodes;
  if (reduced->reducedSolution == 1)                                     
    _eqnOfInterest = multiquation->equationReduced;       
  else                                                                    
    _eqnOfInterest = multiquation->equationLSS;  //equaitonLSS ?
  systemSize = _eqnOfInterest.numStates*Psi2->self->space.node.count;

  VecCreateSeq(PETSC_COMM_SELF, systemSize, &v_petsc);
  if (reduced->reducedSolution == 1)                                       
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, systemSize, NULL, &dfdu); 
  else       
    MatCreateSeqBAIJ(PETSC_COMM_SELF, elementSize, systemSize, systemSize, 4,
		     NULL, &dfdu); 

  //---------------------------------------------------------------------------
  // C0nvert the current solution for v to Vector Petsc form
  //---------------------------------------------------------------------------
  array2Vec(_eqnOfInterest, Psi2->self, v_petsc);  
  //---------------------------------------------------------------------------
  // Calculation of the spatial residual
  //---------------------------------------------------------------------------
  ks_dres_tangent(ykflow, multiquation, Psi2, dfdu, reduced); //find dfdu
  MatMult(dfdu, v_petsc, residual);
  /* //--------------------------------------------------------------------------- */
  /* // Additional term for sensitivity calculation */
  /* //--------------------------------------------------------------------------- */
  /* if (reduced->typeLSS == 0){ */
  /*   ks_function(ykflow, Psi2, reduced, 1, dfds); */
  /*   VecAXPY(residual, -Psi2->self->beta, dfds); */
  /* } */
  
  //---------------------------------------------------------------------------
  // Release and Destroy everything
  //---------------------------------------------------------------------------
  MatDestroy(&dfdu);                                                       
  VecDestroy(&v_petsc);
}
