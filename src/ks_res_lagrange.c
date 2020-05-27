
#include "ks_res_lagrange.h"

void ks_res_lagrange(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		     Cluster* Psi1, Vec residual, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int systemSize;
  int elementSize;
  Vec w_petsc;
  Vec v_petsc;
  Vec djdu;
  Vec f;
  Vec vp;
  Mat dfdu = NULL;
  Universe _equation;
  Universe _eqnOfInterest;
  Cluster *Psi2 = (Cluster *) malloc (sizeof(Cluster));
  Galaxy *_stateOfInterest;
  Psi2->self = (Galaxy *) malloc (sizeof(Galaxy));
  //-----------------------------------------------------------------------
  // Mallocing thing and crafting things and initializing things
  //-----------------------------------------------------------------------
  elementSize = multiquation->equation.numStates*
    Psi1->state->basis.nodes;
  _equation = multiquation->equation;
  if (reduced->reducedSolution == 1){
    _eqnOfInterest = multiquation->equationReduced;
    _stateOfInterest = Psi1->stateFull;
  }else{
    _eqnOfInterest = multiquation->equationLSS;
    _stateOfInterest = Psi1->state;
  }
  systemSize = _eqnOfInterest.numStates*Psi1->self->space.node.count;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &w_petsc);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &v_petsc);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &djdu);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &f);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &vp);
  if (reduced->reducedSolution == 1)
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, systemSize, NULL, &dfdu);
  else
    MatCreateSeqBAIJ(PETSC_COMM_SELF, elementSize, systemSize, systemSize, 4,
                     NULL, &dfdu);
    /* MatCreateSeqBAIJ(PETSC_COMM_SELF, Psi1->state->basis.p+1, systemSize, */
    /*                  systemSize, 3, NULL, &dfdu);    */
  //-----------------------------------------------------------------------
  // Convert the current solution for w to Vector Petsc form
  //-----------------------------------------------------------------------
  array2Vec(_eqnOfInterest, Psi1->self, Psi1->self->space,  w_petsc);
  //-----------------------------------------------------------------------
  // Initialize the v solution via the Psi1 Utype
  //-----------------------------------------------------------------------
  ks_copyUtype(Psi1->self, Psi2->self);
  strcpy(Psi2->clusterId, Psi1->clusterId);
  Psi2->self->basis.p = 0;
  Psi2->self->basis.nodes = 1;
  Psi2->self->space.elem.count = 1;
  Psi2->self->space.node.count = 1;
  Psi2->clusterId[3] = '2';
  Psi2->self->utypeId[3] = '2';
  createSystem(_eqnOfInterest, Psi2->self);
  //---------------------------------------------------------------------------
  // Given Psi1, find the corresponding state solution depending on size ofPsi1
  //---------------------------------------------------------------------------
  ks_dres_lagrange(ykflow, multiquation, Psi1, dfdu, reduced);
  MatMult(dfdu, w_petsc, residual); //calculates dfdu^T*w
  //---------------------------------------------------------------------------
  // Calculate the Projection*tangent value: second term
  //---------------------------------------------------------------------------
  ks_readSolution(_eqnOfInterest, Psi2->self, Psi1->self->time.node);
  array2Vec(_eqnOfInterest, Psi2->self, Psi2->self->space, vp);
  ks_function(ykflow, multiquation, Psi1, reduced, f, Psi1->self->time.node);
  projection(_eqnOfInterest, f, vp);
  VecAXPY(residual, -1, vp);
  //---------------------------------------------------------------------------
  // Third term here and add to the other terms
  //---------------------------------------------------------------------------
  ks_calculatedjdu(ykflow, multiquation, _stateOfInterest, djdu,
		   Psi1->self->time.node, reduced);
  VecAXPY(residual, -Psi1->self->beta/Psi1->self->time.globalT_f, djdu);
  //---------------------------------------------------------------------------
  // Destroy everything!
  //---------------------------------------------------------------------------
  VecDestroy(&w_petsc);
  VecDestroy(&v_petsc);
  VecDestroy(&djdu);
  VecDestroy(&f);
  VecDestroy(&vp);
  MatDestroy(&dfdu);
  destroySystem(_eqnOfInterest, Psi2->self);
  free(Psi2->self);
  free(Psi2);
}
