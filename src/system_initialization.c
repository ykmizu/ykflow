#include "system_initialization.h"

void copyUniverse(Universe equation, Universe *equationCopy){
  int i; //initialization for iteration
  equationCopy->c = (double *) malloc (equation.numParams*sizeof(double));
  strcpy(equationCopy->nameEqn, equation.nameEqn);
  equationCopy->numStates = equation.numStates;
  equationCopy->numParams = equation.numParams;
  for (i=0; i<equation.numParams; i++)
    equationCopy->c[i] = equation.c[i];
}

void createSystem(Universe equation, Galaxy *state){
  state->quad.n = 2*state->basis.p+1;      //Quadrature points


  state->time.count = (state->time.t_f-state->time.t_0)/state->time.dt; //num
  createMesh(&state->space);
  assignBasisQuad(state); //Establish the basis functions and quadrature point
  initializeStates(equation, state);
  sprintf(state->id, "%s_%s_%d_%d", equation.nameEqn, state->utypeId,
          state->basis.p, state->time_method);
}

void copySystem(Universe equation, Galaxy *state, Galaxy *stateCopy){
  stateCopy->quad.n = state->quad.n;
  stateCopy->time.count = state->time.count;
  copyMesh(&state->space, &stateCopy->space);
  assignBasisQuad(stateCopy);
  initializeStates(equation, stateCopy);
  strcpy(stateCopy->id, state->id);
}

/* void findApproxSystem(Universe equation, Galaxy *state, Galaxy *stateApprox, */
/* 		      Is_it reduced){ */
/*   stateApprox->quad.n = state->quad.n; */
/*   stateApprox->time.count = state->time.count; */
/*   findApproxMesh(state, stateApprox, reduced);  */
/*   assignBasisQuad(stateApprox); */
/*   initializeStates(equation, stateApprox); */
/*   strcpy(stateApprox->id, state->id); */
/* } */

/* void createCluster(Multiverse multiquation, Galaxy *primal, Cluster *object){  */
/*   ks_copyUtype(primal, object->self); */
/*   copySystem(primal, object->self); */
/*   ks_copyUtye (primal, object->stateApprox); */
/*   findApproxSystem(primal, Psi); */
/* } */
void createAdjointCluster(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Cluster *primal, Cluster *Psi, Is_it *reduced){
  PetscInt size;
  Psi->self = (Galaxy *) malloc (sizeof(Galaxy));
  Psi->state = (Galaxy *) malloc (sizeof(Galaxy));
  //-----------------------------------------------------------------
  // Initialization
  //----------------------------------------------------------------
  Psi->self->beta = primal->self->beta;
  //--------------------------------------------------------------------
  // Implementation
  //-------------------------------------------------------------------
  if (ykflow->numElemCol == 3){
    MatCreateSeqDense(PETSC_COMM_SELF, primal->self->basis.nodes,
		      primal->self->basis.nodes, NULL, &Psi->state->Mij);
    MatCopy(primal->self->Mij, Psi->state->Mij, SAME_NONZERO_PATTERN);
  }
  ks_copyUtype(primal->self, Psi->state);
  copySystem(multiquation->equation, primal->self, Psi->state);
  if (reduced->reducedSolution == 1){  //Reduced Solution
    //    VecGetSize(primal->reduced->djdu, &size);
    //VecCreateSeq(PETSC_COMM_SELF, size, &Psi->self->djdu);
    //VecCopy(primal->reduced->djdu, Psi->self->djdu);
    Psi->reduced = (Galaxy *) malloc (sizeof(Galaxy));
    Psi->stateFull = (Galaxy *) malloc (sizeof(Galaxy));
    ks_copyUtype(primal->reduced, Psi->reduced);
    ks_copyUtype(primal->stateFull, Psi->stateFull);
    copySystem(multiquation->equationReduced, primal->reduced, Psi->reduced);
    copySystem(multiquation->equation, primal->stateFull, Psi->stateFull);
  }else{  //Full Solution here
    //VecGetSize(primal->self->djdu, &size);
    //VecCreateSeq(PETSC_COMM_SELF, size, &Psi->self->djdu);
    //VecCopy(primal->self->djdu, Psi->self->djdu);
  }
}

void destroyAdjointCluster(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			   Cluster *primal, Cluster *Psi, Is_it *reduced){
  destroySystem(multiquation->equation, Psi->state);
  if (ykflow->numElemCol == 3)
    MatDestroy(&Psi->state->Mij);
  free(Psi->state);
  free(Psi->self);
  if (reduced->reducedSolution == 1){
    destroySystem(multiquation->equationReduced, Psi->reduced);
    destroySystem(multiquation->equation, Psi->stateFull);
    free(Psi->reduced);
    free(Psi->stateFull);
  }
}


/* void createAdjointCluster(Multiverse multiquation, Cluster *primal, */
/* 			  Cluster *Psi){ */
/*   PetscInt size; */
/*   VecGetSize(primal->reduced->djdu, &size); */
/*   Psi->self = (Galaxy *) malloc (sizeof(Galaxy));                             */
/*   Psi->reduced = (Galaxy *) malloc (sizeof(Galaxy)); */
/*   Psi->state = (Galaxy *) malloc (sizeof(Galaxy));                   */
/*   Psi->stateFull = (Galaxy *) malloc (sizeof(Galaxy));    */
/*   Psi->stateFull->beta = primal->stateFull->beta; */
/*   MatCreateSeqDense(PETSC_COMM_SELF, primal->self->basis.p+1, */
/* 		    primal->self->basis.p+1, NULL, &Psi->state->Mij); */
/*   VecCreateSeq(PETSC_COMM_SELF, size, &Psi->self->djdu);  */
/*   ks_copyUtype(primal->self, Psi->state); */
/*   ks_copyUtype(primal->reduced, Psi->reduced); */
/*   ks_copyUtype(primal->stateFull, Psi->stateFull); */
/*   VecCopy(primal->reduced->djdu, Psi->self->djdu);  */
/*   MatCopy(primal->self->Mij, Psi->state->Mij, SAME_NONZERO_PATTERN); */
/*   copySystem(multiquation.equation, primal->self, Psi->state); */
/*   copySystem(multiquation.equationReduced, primal->reduced, Psi->reduced); */
/*   copySystem(multiquation.equation, primal->stateFull, Psi->stateFull); */
/* } */

void initializeStates(Universe equation, Galaxy *state){
  int i;               //initialization for iteration
  state->solution = (Array *) malloc (equation.numStates*sizeof(Array));
  for (i=0; i<equation.numStates; i++)
    initArray(&state->solution[i], state->space.node.count);
}

void copyMesh(Mesh *space, Mesh *spaceCopy){
  int i;
  int fullElementCount;
  int fullNodeCount;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  fullElementCount = space->elem.count;
  fullNodeCount = space->node.count;
  initArrayInt(&spaceCopy->elem, fullElementCount);
  initArrayInt(&spaceCopy->node, fullNodeCount);
  for (i=0; i<space->elem.count; i++)
    spaceCopy->elem.array[i] = space->elem.array[i];
  for (i=0; i<space->node.count; i++)
    spaceCopy->node.array[i] = space->node.array[i];
}



void createMesh(Mesh *space){
  int i, j;                       //initialization for iteration
  int fullElementCount;
  int fullNodeCount;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  fullElementCount = space->elem.count;
  fullNodeCount = space->node.count;
  initArrayInt(&space->elem, fullElementCount);
  initArrayInt(&space->node, fullNodeCount);

  for (i=0; i<space->elem.count; i++)
    space->elem.array[i] = i;
  for (i=0; i<space->node.count; i++)
    space->node.array[i] = i;

}



void assignBasisQuad(Galaxy *state){
  int i, j;                           //initialization for iteration
  initArray(&state->quad.x_i, state->quad.n); //Quadrature points: primal sol
  initArray(&state->quad.w, state->quad.n);   //Quadrature weights: primal sol
  initArray2D(&state->basis.phi, state->basis.p+1, state->quad.n);  //ϕ
  initArray2D(&state->basis.dphi, state->basis.p+1, state->quad.n); //dϕ/dɛ
  quadrature(state->quad);                //Primal points

  for (i=0; i<state->basis.p+1; i++){
    for (j=0; j<state->quad.n; j++){
      state->basis.phi.array2[i][j] =
        phi(state->basis.p, i, state->quad.x_i.array[j]);
      state->basis.dphi.array2[i][j] =
        dphi(state->basis.p, i, state->quad.x_i.array[j]);
    }
  }
}

void destroySystem(Universe equation, Galaxy *state){
  destroyMesh(&state->space);
  destroyStates(equation, state);
  deleteBasisQuad(state);
}

void destroyStates(Universe equation, Galaxy *state){
  int i;
  for (i=0; i<equation.numStates; i++){
    delArray(&state->solution[i]);
  }
  free(state->solution);
}

void destroyMesh(Mesh *space){
  //delArray(&state->space.x);
  //space->elem.array = (int *)malloc(fullElementCount*sizeof(int));
  delArrayInt(&space->elem);
  delArrayInt(&space->node);
  //free(space->elem.array);
  //free(space->node.array);
  //delArrayInt(&state->space.elem);
  //delArrayInt(&state->space.node);
}

void deleteBasisQuad(Galaxy *state){
  delArray(&state->quad.x_i);
  delArray(&state->quad.w);
  delArray2D(&state->basis.phi);
  delArray2D(&state->basis.dphi);
}
