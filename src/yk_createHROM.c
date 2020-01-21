



#include "yk_createHROM.h"



void yk_createReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				   Multiverse *multiquation, Cluster *primal,
				   Cluster *primalApprox, Is_it *reduced){

  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  char approximateID[1000];             // Name given to the approximate state
  char reducedID[1000];                 // Name given to the reduced state
  Mat snapshot;                         // Snapshot Matrix before POD
  Mat *st_basis_temp;
  Mat *st_snapshot;
  int snapshotCount = primal->self->time.count;
  /* int snapshotCount = primal->self->time.count_wi; */
  /* int snapshotCount = ceil(((primal->self->time.t_f-primal->self->time.t_0)*10)/ */
  /* 			   (primal->self->time.dt*10)); */
  Universe _equationFull = multiquation->equation;
  /* Universe *_equationReduced = &multiquation->equationReduced; */
  Mat *s_mu;
  Vec init;
  /* Mat *tryrOB; */
  /* reduced->ST_rOBState_i = tryrOB; */
  PetscInt rows, nBasis;
  /* PetscMalloc1(reduced->nBasisFuncs, &st_basis_temp); */

  reduced->ST_rOBState_i = (Mat *) malloc (snapshotCount*sizeof(Mat));
  /* PetscMalloc1(primal->self->time.count, &tryrOB); */
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  strcpy(approximateID, "ROM_approximate");
  strcpy(reducedID, "ROM_reduced");
  strcpy(primalApprox->clusterId, primal->clusterId);

  //---------------------------------------------------------------------------
  // Create reduced model calculation for state
  //---------------------------------------------------------------------------
  s_mu = yk_createSnapshotState(ykflow, multiquation, primal, &snapshot,
				reduced->dss, reduced);
  yk_properOrthogonalDecompose(&snapshot, primal->self->systemSize,
			       primal->self->index, &reduced->nBasisFuncs,
			       reduced->eBasisSpace, &reduced->rOBState);
  //---------------------------------------------------------------------------
  //Calculate the space time basis and make it be the reduced->ST_rOBState
  //---------------------------------------------------------------------------
  MatGetSize(reduced->rOBState, &rows, &nBasis);
  PetscMalloc1(nBasis, &st_snapshot);
  PetscMalloc1(nBasis, &st_basis_temp);
  yk_createTemporalSnapshotState(ykflow, multiquation, primal, st_snapshot,
				 s_mu,
				 reduced->dss, reduced);
  yk_properOrthogonalDecomposeGroup(st_snapshot, nBasis, snapshotCount,
				    primal->self->index, &reduced->nBasisTime,
				    reduced->eBasisTime, st_basis_temp);

  //---------------------------------------------------------------------------
  // Callcualte the ST_HOSVD Tailored space-time basis
  //---------------------------------------------------------------------------
  yk_createSpaceTimeBasis(&reduced->rOBState, st_basis_temp,
			  &reduced->ST_rOBState, reduced->ST_rOBState_i);
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->self->systemSize = primal->self->systemSize;
  primalApprox->self->index = (PetscInt *) malloc
    (primalApprox->self->systemSize*sizeof(PetscInt));
  for (i=0; i<primalApprox->self->systemSize; i++)
    primalApprox->self->index[i] = i;
  createSystem(multiquation->equation, primalApprox->self);
  reduced->reducedMesh = primalApprox->self->space;
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  /* sprintf(eqnReduced, "%s%s", jobFile, "Reduced"); */
  /* strcpy(multiquation->equationReduced.nameEqn, eqnReduced); */
  multiquation->equationReduced.numStates = nBasis;
  //---------------------------------------------------------------------------
  // Create the reduced states
  //---------------------------------------------------------------------------
  //THis part is questionable
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(multiquation->equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  //CreateInitial conditions based on s_mu
  //---------------------------------------------------------------------------
  createInitialConditions_ST(ykflow, primalApprox, s_mu, &init, reduced);
  vec2Array(multiquation->equationReduced, primalApprox->reduced, init);
  ks_printSolution(multiquation->equationReduced, primalApprox->reduced, 0);

  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  for (i=0; i<nBasis; i++){
    MatDestroy(&st_snapshot[i]);
    MatDestroy(&st_basis_temp[i]);
  }

  PetscFree(st_snapshot);
  PetscFree(st_basis_temp);
  VecDestroy(&init);
  MatDestroy(&snapshot);
  for (i=0; i<reduced->numParamSet; i++)
    MatDestroy(&s_mu[i]);
  PetscFree(s_mu);
}

void yk_runReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int innertimesolver = 1;              //time solver
  Vec primal_0;
  Vec primalReducedVec;
  Universe _equationFull = multiquation->equation;
  PetscInt row, basis;
  int leftTimeNode;
  Universe _equationReduced = multiquation->equationReduced;
  MatGetSize(reduced->ST_rOBState, &row, &basis);
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->nSnapshotsRJ = 0;
  VecCreate(PETSC_COMM_SELF, &primal_0);//Set up the singular vector
  VecSetSizes(primal_0, primal->self->systemSize, primal->self->systemSize);
  VecSetType(primal_0, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, basis, basis);
  VecSetType(primalReducedVec, VECSEQ);

  fclose(fopen("residualSnapshot.dat", "w"));
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------
  leftTimeNode = primal->self->time.t_0/primal->self->time.dt;
  if (leftTimeNode == 0){
    ks_readSolution(_equationFull, primal->self, leftTimeNode);
    ks_copySolutions(_equationFull, primal->self, primalApprox->self);
  } else{
    ks_readSolution(_equationFull, primalApprox->self, leftTimeNode);
  }
  ks_printSolution(_equationFull, primalApprox->self, leftTimeNode);
  //---------------------------------------------------------------------------
  // Gauss NEwton to minimize 1/2 norm(Residua, 2)^2
  //---------------------------------------------------------------------------
  yk_gaussNewtonSolve_ST(ykflow, multiquation, primalApprox, reduced,
  		      innertimesolver, &reduced->nSnapshotsRJ);
  VecDestroy(&primal_0);
  VecDestroy(&primalReducedVec);
}

void yk_storeBasisInfoWindow(int wnum_i, Is_it *reduced){
  if (reduced->nTimeSegs>1){
    reduced->win_i[wnum_i]->ST_rOBState = reduced->ST_rOBState;
    reduced->win_i[wnum_i]->ST_rOBState_i = reduced->ST_rOBState_i;
  }
}

void yk_destroyReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				    Multiverse *multiquation, Cluster *primal,
				    Cluster *primalApprox, Is_it *reduced){
  int i; //initialization for iteration
  int snapshotCount = primal->self->time.count;
  Universe _equationReduced = multiquation->equationReduced;
  MatDestroy(&reduced->rOBState);
  for (i=0; i<snapshotCount; i++)
    MatDestroy(&reduced->ST_rOBState_i[i]);
  free(reduced->ST_rOBState_i);
  MatDestroy(&reduced->ST_rOBState);
  free(primalApprox->self->index);
  destroySystem(_equationReduced, primalApprox->reduced);
  destroySystem(multiquation->equation, primalApprox->self);
  free(primalApprox->reduced);
  free(primalApprox->self);
}


//-----------------------------------------------------------------------------
// Initialize and set up everything for the potential to run a ROM
//-----------------------------------------------------------------------------
void yk_createReducedOrderModel(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  char approximateID[1000];             // Name given to the approximate state
  char reducedID[1000];                 // Name given to the reduced state
  Mat snapshot;                         // Snapshot Matrix before POD
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  strcpy(approximateID, "ROM_approximate");
  strcpy(reducedID, "ROM_reduced");
  strcpy(primalApprox->clusterId, primal->clusterId);

  //---------------------------------------------------------------------------
  // Create reduced model calculation for states
  //---------------------------------------------------------------------------
  yk_createSnapshotState(ykflow, multiquation, primal,&snapshot, reduced->dss,
			 reduced);
  yk_properOrthogonalDecompose(&snapshot, primal->self->systemSize,
			       primal->self->index, reduced->nBasisFuncs,
			       reduced->eBasisSpace, &reduced->rOBState);
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->self->systemSize = primal->self->systemSize;
  primalApprox->self->index = (int *) malloc
    (primalApprox->self->systemSize*sizeof(int));
  for (i=0; i<primalApprox->self->systemSize; i++)
    primalApprox->self->index[i] = i;
  createSystem(multiquation->equation, primalApprox->self);
  reduced->reducedMesh = primalApprox->self->space;
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(_equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&snapshot);
}

//-----------------------------------------------------------------------------
// Execute to find the reduced order model
//-----------------------------------------------------------------------------
void yk_runReducedOrderModel(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Cluster *primalApprox,
			     Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int innertimesolver = 1;              //time solver
  Vec primal_0;
  Vec primalReducedVec;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->nSnapshotsRJ = 0;
  VecCreate(PETSC_COMM_SELF, &primal_0);//Set up the singular vector
  VecSetSizes(primal_0, primal->self->systemSize, primal->self->systemSize);
  VecSetType(primal_0, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, reduced->nBasisFuncs, reduced->nBasisFuncs);
  VecSetType(primalReducedVec, VECSEQ);

  fclose(fopen("residualSnapshot.dat", "w"));
  fclose(fopen("jacobianSnapshot.dat", "w"));
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  ks_copySolutions(_equationFull, primal->self, primalApprox->self);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // Set r0 initial solution of reduced solution to 0
  //---------------------------------------------------------------------------
  /* array2Vec(_equationFull, primal->self, primal_0); */
  /* MatMultTranspose(reduced->rOBState, primal_0, primalReducedVec); */
  VecSet(primalReducedVec, 0);
  vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);
  //---------------------------------------------------------------------------
  // Gauss NEwton to minimize 1/2 norm(Residua, 2)^2
  //---------------------------------------------------------------------------
  yk_gaussNewtonSolve(ykflow, multiquation, primalApprox, reduced,
		      innertimesolver, &reduced->nSnapshotsRJ);
  VecDestroy(&primal_0);
  VecDestroy(&primalReducedVec);

}

void yk_createHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				     Multiverse *multiquation, Cluster *primal,
				     Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int seedCount;                        // Number of Boundary Nodes
  int *nodeSet = (int *) malloc (reduced->nSampleNodes*sizeof(int));
  char approximateID[1000];
  char reducedID[1000];
  char realID[1000];
  char *residualSnapshots = "residualSnapshot.dat";
  char *jacobianSnapshots = "jacobianSnapshot.dat";
  Mat sResMat, sJacMat;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 1;
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->stateFull = (Galaxy *) malloc (sizeof(Galaxy));
  strcpy(approximateID, "HROM_approximate");
  strcpy(reducedID, "HROM_reduced");
  strcpy(realID, "HROM_full");
  strcpy(primalApprox->clusterId, primal->clusterId);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for residuals
  //---------------------------------------------------------------------------
  ks_readMatrix(residualSnapshots, primal->self->systemSize, &sResMat);
  yk_properOrthogonalDecompose(&sResMat, primal->self->systemSize,
			       primal->self->index, reduced->nBasisFuncsRJ,
			       reduced->eBasisSpace, &reduced->rOBResidual);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for jacobians
  //---------------------------------------------------------------------------
  ks_readMatrix(jacobianSnapshots, primal->self->systemSize, &sJacMat);
  yk_properOrthogonalDecompose(&sJacMat, primal->self->systemSize,
  			       primal->self->index, reduced->nBasisFuncsRJ,
  			       reduced->eBasisSpace, &reduced->rOBJacobian);
  printf("POD DONE\n");
  //---------------------------------------------------------------------------
  // Greedy Algorithm and create of the reduced mesh
  //---------------------------------------------------------------------------
  ykflow->boundarySeeds(ykflow, multiquation, primal, reduced->nSampleNodes,
			&seedCount, nodeSet);
  printf("here\n");
  ks_greedyAlgorithm(seedCount, nodeSet, reduced);
  printf("mayube here\n");
  minElements(multiquation, primal, nodeSet, reduced);
  for (i=0; i<reduced->reducedMesh.elem.count; i++)
    printf("%d\n", reduced->reducedMesh.elem.array[i]);
  //---------------------------------------------------------------------------
  // Create the hyper approximate states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->stateFull->systemSize = primal->self->systemSize;
  primalApprox->stateFull->index = (int *) malloc
    (primalApprox->stateFull->systemSize*sizeof(int));
  for (i=0; i<primalApprox->stateFull->systemSize; i++)
    primalApprox->stateFull->index[i] = i;
  findApproxMesh(ykflow, multiquation, primal, primalApprox, reduced);
  assignBasisQuad(primalApprox->self);
  initializeStates(_equationFull, primalApprox->self);
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(_equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Create the hyper real states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->stateFull);
  createSystem(_equationFull, primalApprox->stateFull);
  sprintf(primalApprox->stateFull->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->stateFull->basis.p,
          primalApprox->stateFull->time_method, realID);
  //---------------------------------------------------------------------------
  // Create the corresponding identity matrices from greedy algorithm
  //---------------------------------------------------------------------------
  restrictedIdentity(_equationFull, primal->self, primalApprox->self, nodeSet,
		     reduced);
  MatMatMult(reduced->Zp, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
             &reduced->rOBStateBar);
  MatMatMult(reduced->Zu, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
             &reduced->rOBStateRed);
  //---------------------------------------------------------------------------
  // Destroy and Free things
  //---------------------------------------------------------------------------
  free(nodeSet);
  MatDestroy(&sResMat);
  MatDestroy(&sJacMat);
}

void yk_runHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				  Multiverse *multiquation, Cluster *primal,
				  Cluster *primalApprox,Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int systemApproxSize;                 //Size of the primalApprox vector
  int innertimesolver = 1;              //Only for 1D... gotta figure out this
  Vec primalApprox_0;
  Vec primal_0;
  Vec primalReducedVec;
  Vec primalRealVec;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialization and Mallocing
  //---------------------------------------------------------------------------
  reduced->hrom = 1;
  systemApproxSize =
    _equationFull.numStates*primalApprox->self->space.node.count;
  VecCreate(PETSC_COMM_SELF, &primal_0);
  VecSetSizes(primal_0, primal->self->systemSize, primal->self->systemSize);
  VecSetType(primal_0, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalRealVec);
  VecSetSizes(primalRealVec, primal->self->systemSize,
	      primal->self->systemSize);
  VecSetType(primalRealVec, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, reduced->nBasisFuncs, reduced->nBasisFuncs);
  VecSetType(primalReducedVec, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalApprox_0);
  VecSetSizes(primalApprox_0, systemApproxSize, systemApproxSize);
  VecSetType(primalApprox_0, VECSEQ);
  //---------------------------------------------------------------------------
  // Find A and B for the calculation of residual and jacobian approx
  //---------------------------------------------------------------------------
  ks_offlineMatrix(_equationFull, primal->self, primalApprox->self,
		   reduced);
  //---------------------------------------------------------------------------
  // Print out the initial conditions for completeness (primalApprox);
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  array2Vec(_equationFull, primal->self, primal_0);
  MatMult(reduced->Zp, primal_0, primalApprox_0);
  vec2Array(_equationFull, primalApprox->self, primalApprox_0);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // Find the r0 initial solution of reduced solution
  //---------------------------------------------------------------------------


  VecSet(primalReducedVec, 0);
  vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);


  /* array2Vec(_equationFull, primal->self, primal_0); */
  /* MatMultTranspose(reduced->rOBState, primal_0, primalReducedVec); */
  /* vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /* ks_printSolution(_equationReduced, primalApprox->reduced, 0); */

  //---------------------------------------------------------------------------
  // Read Solution initial conditions
  //---------------------------------------------------------------------------
  ks_copySolutions(_equationFull, primal->self, primalApprox->stateFull);
  ks_printSolution(_equationFull, primalApprox->stateFull, 0);
  //---------------------------------------------------------------------------
  // Gauss Newton Solver
  //---------------------------------------------------------------------------
  yk_gaussNewtonSolve(ykflow, multiquation, primalApprox, reduced,
		      innertimesolver, &reduced->nSnapshotsRJ);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  VecDestroy(&primal_0);
  VecDestroy(&primalRealVec);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalApprox_0);
  MatDestroy(&reduced->A);
  MatDestroy(&reduced->B);
}
void minElements(Multiverse *multiquation, Cluster *primal, int *nodeSet,
		 Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             // initialization for iteration
  int elem;                             // elem of interest given nSampleNodes
  int count = 0;
  int totNumElems = 0;                  // saves total number of elements to c
  int numBasis = primal->self->basis.nodes;
  int elemSize = numBasis*multiquation->equation.numStates;

  int *tempNodeUSet = (int *) malloc     // array info on which element to keep
    (primal->self->space.elem.count*sizeof(int));
  //---------------------------------------------------------------------------
  // Implementation given the greedyalgorithm giving me nodes to keep
  //---------------------------------------------------------------------------
  for (i=0; i<primal->self->space.elem.count; i++)
    tempNodeUSet[i] = 0;
  for (i=0; i<reduced->nSampleNodes; i++){ //tells me which element included
    elem = floor(nodeSet[i]/(elemSize)); //nodeSet has only samplenodes
    tempNodeUSet[elem] = 1;             // include this element
  }
  // Count the total of elements in the reduced mesh
  for (i=0; i<primal->self->space.elem.count; i++)
    totNumElems += tempNodeUSet[i];
  //Need to fill in has table thingie in xflow

  //---------------------------------------------------------------------------
  // Create the actual temporary mesh
  //---------------------------------------------------------------------------
  //Take this out later, this is weird
  initArrayInt(&reduced->reducedMesh.elem, totNumElems);
  initArrayInt(&reduced->reducedMesh.node, numBasis*totNumElems);
  //Don't need the dx/x thing. Want to be indepedent on the physical mesh trait
  //Make sure this is possible
  for (i=0; i<primal->self->space.elem.count; i++){
    if (tempNodeUSet[i] == 1){ //If this element exists in reduced mesh
      reduced->reducedMesh.elem.array[count] = i;
      for (j=0; j<numBasis; j++)
	reduced->reducedMesh.node.array[count*numBasis+j] = i*numBasis+j;
      count++;
    }
  }
  free(tempNodeUSet);
}

void findApproxMesh(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *primal, Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             // initialization for iteration
  int totNumElems = 0;
  int count = 0;
  int elemCount = primal->self->space.elem.count;
  int numBasis = primal->self->basis.nodes;
  int *elemSet = (int *) malloc (elemCount*sizeof(int));
  //---------------------------------------------------------------------------
  // Count the number of total elements in the approximate mesh given reduced m
  //---------------------------------------------------------------------------
  for (i=0; i<elemCount; i++)
    elemSet[i] = 0; //initialize it to empty 0
  for (i=0; i<reduced->reducedMesh.elem.count; i++)
    ykflow->adjacentElems(ykflow, primal, reduced->reducedMesh.elem.array[i],
			  elemSet);
  //Count the number of total elements in the approximate mesh
  for (i=0; i<elemCount; i++)
    totNumElems += elemSet[i];
  //---------------------------------------------------------------------------
  // Now create that approximate mesh
  //---------------------------------------------------------------------------
  initArrayInt(&primalApprox->self->space.elem, totNumElems);
  initArrayInt(&primalApprox->self->space.node, numBasis*totNumElems);
  //Don't need the dx/x thing. Want to be indepdent off the physical mesh
  for (i=0; i<primal->self->space.elem.count; i++){
    if (elemSet[i] == 1){ //If it exists
      primalApprox->self->space.elem.array[count] = i;
      for (j=0; j<numBasis; j++)
        primalApprox->self->space.node.array[count*numBasis+j] = i*numBasis+j;
      count++;
    }
  }
  free(elemSet);
}


void yk_findBoundarySeeds1D(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, int nSampNodes, int *numSeeds,
			    int *nodeSet){
  *numSeeds = 2;
  nodeSet[0] = 0;
  nodeSet[multiquation->equation.numStates] =  primal->self->systemSize-1;
}

void yk_findAdjacentElems1D(yk_PrimalSolver *ykflow, Cluster *primal, int elem,
			    int *elemSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  //1D
  elemSet[elem]= 1;
  if (elem >0 )
    elemSet[elem-1] = 1;
  if (elem < primal->self->space.elem.count-1)
    elemSet[elem+1] = 1;
}


void createInitialConditions_ST(yk_PrimalSolver *ykflow, Cluster *primal,
				Mat * s_mu, Vec *estimatedinitST,
				Is_it *reduced){

  PetscInt i,j; //initialization for iteration
  Mat sumS_MU;
  Vec vecsnapshot;
  PetscInt *index = (PetscInt *) malloc (primal->self->systemSize*
					 sizeof(PetscInt));
  PetscScalar *colvalues = (PetscScalar *) malloc (primal->self->systemSize*
						   sizeof(PetscScalar));
  PetscInt row, nBasis;
  MatGetSize(reduced->ST_rOBState, &row, &nBasis);
  PetscInt snapshotCount = primal->self->time.count;
  /* PetscInt snapshotCount = primal->self->time.count_wi; */
  //---------------------------------------------------------------------------
  // initialization
  //---------------------------------------------------------------------------
  MatCreate(PETSC_COMM_SELF,&sumS_MU);
  MatSetSizes(sumS_MU,primal->self->systemSize,snapshotCount,
	      primal->self->systemSize,snapshotCount);
  MatSetType(sumS_MU,MATDENSE);
  MatSetUp(sumS_MU);


  VecCreate(PETSC_COMM_SELF, &vecsnapshot);
  VecSetSizes(vecsnapshot, snapshotCount*primal->self->systemSize,
	      snapshotCount*primal->self->systemSize);
  VecSetType(vecsnapshot, VECSEQ);

  VecCreate(PETSC_COMM_SELF, estimatedinitST);
  VecSetSizes(*estimatedinitST, nBasis, nBasis);
  VecSetType(*estimatedinitST, VECSEQ);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->numParamSet; i++)
    MatAXPY(sumS_MU, 1.00, s_mu[i], SAME_NONZERO_PATTERN);

  MatAssemblyBegin(sumS_MU,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(sumS_MU,MAT_FINAL_ASSEMBLY);

  MatScale(sumS_MU, 1.0/reduced->numParamSet);

  //Create the vector of the snapshots for the average of the parameter stuff
  for (i=0; i<snapshotCount; i++){
    for (j=0; j<primal->self->systemSize; j++)
      index[j] = i*primal->self->systemSize + j;
    MatGetValues(sumS_MU, primal->self->systemSize, primal->self->index,
		 1, &i, colvalues);
    VecSetValues(vecsnapshot, primal->self->systemSize, index, colvalues,
		 INSERT_VALUES);
  }
  VecAssemblyBegin(vecsnapshot);
  VecAssemblyEnd(vecsnapshot);

  //Find the reduced states for the average snapshots
  //This is the initial condition. Print to file to be used.
  MatMultTranspose(reduced->ST_rOBState, vecsnapshot, *estimatedinitST);

  MatDestroy(&sumS_MU);
  VecDestroy(&vecsnapshot);
  free(index);
  free(colvalues);
}
