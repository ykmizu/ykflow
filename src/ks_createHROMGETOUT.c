
#include "ks_createHROM.h"
                                                                               
void createReducedOrderModel(Multiverse multiquation, Cluster *primal,
			     Cluster *primalApprox, Is_it *reduced){
  int i;                            //initialization for iteration
  int systemSize;
  int innertimesolver = 1;
  char approximateID[1000];
  char reducedID[1000];
  Vec primal_0;                        //initial conditions for full state
  Vec primalApprox_0;
  Vec primalReducedVec;
  Vec primalApproxVec;
  Mat snapshot;                        //snapshot matrix for primal solution
  Universe _equationFull= multiquation.equation;
  Universe _equationReduced = multiquation.equationReduced;
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->nSnapshotsRJ = 0;
  strcpy(approximateID, "ROM_approximate");
  strcpy(reducedID, "ROM_reduced");
  strcpy(primalApprox->clusterId, primal->clusterId); 
  systemSize  = primal->self->space.node.count*_equationFull.numStates;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primal_0);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalApprox_0);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &primalReducedVec);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalApproxVec);
  fclose(fopen("residualSnapshot.dat", "w"));
  fclose(fopen("jacobianSnapshot.dat", "w")); 
  //---------------------------------------------------------------------------
  //initializae the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self); // number: Ufine.
  createSystem(multiquation.equation, primalApprox->self);
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
	  primalApprox->clusterId, primalApprox->self->basis.p,
	  primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  //initialize the reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);  //doesn't have createSys
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->space.dx = 0; 
  createSystem(_equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
	  primalApprox->clusterId, primalApprox->self->basis.p,
	  primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for states
  //---------------------------------------------------------------------------
  createSnapshotState(_equationFull, primal->self, &snapshot, reduced->dss); 
  properOrthogonalDecompose(snapshot, systemSize, &reduced->nBasisFuncs,
			    &reduced->rOBState);
  //---------------------------------------------------------------------------
  // For complete sake, find the r0 initial solution of reduced solution
  //---------------------------------------------------------------------------
  array2Vec(_equationFull, primal->self, primal_0);
  MatMultTranspose(reduced->rOBState, primal_0, primalReducedVec);
  vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  ks_copySolutions(_equationFull, primal->self, primalApprox->self);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // Gauss Netwon to minimize 1/2norm(Residual, 2)^2
  //---------------------------------------------------------------------------
  reduced->reducedMesh = primalApprox->self->space;
  for (i=1; i<primalApprox->self->time.count+1; i++){
    primalApprox->self->time.node = i;
    primalApprox->reduced->time.node = i;   //the time reduced solution to gue
    ks_setSolutionZero(_equationReduced, primalApprox->reduced);
    ks_gaussNewtonSolve(multiquation, primalApprox, reduced, innertimesolver,
			&reduced->nSnapshotsRJ);
    if (i == 1 && primalApprox->self->time_method == 2) //Use BDF1
      innertimesolver = 2;
    else if (i == 1 && primalApprox->self->time_method == 3)
      innertimesolver = 2;
    else if (i == 2 && primalApprox->self->time_method == 3)
      innertimesolver = 3; //After two iterations you can now go and use BDF3
    //-------------------------------------------------------------------------
    // Print out the reduced states found at everytime time iteration
    //-------------------------------------------------------------------------
    ks_printSolution(_equationReduced, primalApprox->reduced, i);
    ks_printSolution(_equationFull, primalApprox->self, i);
  }
  //---------------------------------------------------------------------------
  // Free and Destroy and everything
  //---------------------------------------------------------------------------
  VecDestroy(&primal_0);
  VecDestroy(&primalApprox_0);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalApproxVec);
  MatDestroy(&snapshot);
}




void createHyperReducedOrderModel(Multiverse multiquation, Cluster *primal,
				  Cluster *primalApprox, Is_it *reduced){
  int i;                            //initialization for iteration
  int iterationCount = 0;
  int systemSize;
  int systemApproxSize;
  int reducedMeshSize;
  int seedCount = 2; //boundary condition values 2 of them in 2D greedy
  int innertimesolver = 1;             //signifying BDF1 BDF2 etc
  char approximateID[1000];
  char reducedID[1000];
  char realID[1000];
  char *residualSnapshots = "residualSnapshot.dat";
  char *jacobianSnapshots = "jacobianSnapshot.dat";
  PetscInt *nodeSet =
    (PetscInt*) malloc (reduced->nSampleNodes*sizeof(PetscInt));
  Vec primalApprox_0;
  Vec primal_0;
  Vec primalReducedVec;
  Vec primalRealVec;
  Vec primalApproxVec;
  Mat sResMat, sJacMat;
  Universe _equationFull= multiquation.equation;                             
  Universe _equationReduced = multiquation.equationReduced;  
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy)); //Reduced solutio
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->stateFull = (Galaxy *) malloc (sizeof(Galaxy));
  //---------------------------------------------------------------------------
  // Initialization and Mallocing
  //---------------------------------------------------------------------------
  reduced->hrom = 1;
  strcpy(approximateID, "HROM_approximate");                            
  strcpy(reducedID, "HROM_reduced");                                    
  strcpy(realID, "HROM_full");
  strcpy(primalApprox->clusterId, primal->clusterId);
  systemSize  = primal->self->space.node.count*_equationFull.numStates; //FOM 
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primal_0);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalRealVec);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalApproxVec);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &primalReducedVec);
  //---------------------------------------------------------------------------
  //initialize the hyper reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);  //doesn't have createSys
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->space.dx = 0;
  createSystem(_equationReduced, primalApprox->reduced);  
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primal->self->basis.p,
	  primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  //initialize the hyper real states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->stateFull);
  createSystem(_equationFull, primalApprox->stateFull);
  sprintf(primalApprox->stateFull->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->stateFull->basis.p,
	  primalApprox->stateFull->time_method, realID);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for residuals
  //---------------------------------------------------------------------------
  ks_readMatrix(residualSnapshots, systemSize, reduced->nSnapshotsRJ,&sResMat);
  properOrthogonalDecompose(sResMat, systemSize, &reduced->nBasisFuncsRJ,
			    &reduced->rOBResidual);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for jacobians
  //---------------------------------------------------------------------------
  ks_readMatrix(jacobianSnapshots, systemSize, reduced->nSnapshotsRJ,&sJacMat);
  properOrthogonalDecompose(sJacMat, systemSize, &reduced->nBasisFuncsRJ,
			    &reduced->rOBJacobian);
  //---------------------------------------------------------------------------
  // Greedy Algorithim
  //---------------------------------------------------------------------------
  nodeSet[0] = 0;
  nodeSet[1] = systemSize-1;
  ks_greedyAlgorithm(seedCount, nodeSet, reduced);
  //---------------------------------------------------------------------------
  // Creation of the reduced mesh
  //---------------------------------------------------------------------------
  minElements(primal->self, nodeSet, reduced);
  reducedMeshSize = reduced->reducedMesh.node.count*_equationFull.numStates;
  //---------------------------------------------------------------------------
  // initialize the hyper approximate states   
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  findApproxMesh(primal->self, primalApprox->self, reduced);
  assignBasisQuad(primalApprox->self);
  initializeStates(_equationFull, primalApprox->self);
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
	  primalApprox->self->time_method, approximateID);
  systemApproxSize =
    _equationFull.numStates*primalApprox->self->space.node.count;
  VecCreateSeq(PETSC_COMM_SELF, systemApproxSize, &primalApprox_0);
  //---------------------------------------------------------------------------
  // Find A and B for the calculation of residual and jacobian approx 4 p
  //---------------------------------------------------------------------------
  ks_offlineMatrix(_equationFull, primal->self, primalApprox->self, nodeSet,
		   reduced);
  //---------------------------------------------------------------------------
  // Print out the inidital conditions for completeness (primalApprox)
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  array2Vec(_equationFull, primal->self, primal_0);
  MatMult(reduced->Zp, primal_0, primalApprox_0);
  vec2Array(_equationFull, primalApprox->self, primalApprox_0);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // For complete sake, find the r0 initial solution of reduced solution
  //---------------------------------------------------------------------------
  MatMultTranspose(reduced->rOBStateBar, primalApprox_0, primalReducedVec);
  vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);
  //---------------------------------------------------------------------------
  // Real Solution initial conditions (primal->approx->self in ROM)
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);                            
  ks_copySolutions(_equationFull, primal->self, primalApprox->stateFull);
  ks_printSolution(_equationFull, primalApprox->stateFull, 0);  
  //---------------------------------------------------------------------------
  //Gauss Newton Solver for hyper reduced
  //---------------------------------------------------------------------------
  for (i=1; i<primalApprox->self->time.count+1; i++){
    primalApprox->self->time.node = i;
    primalApprox->reduced->time.node = i;
    ks_setSolutionZero(_equationReduced, primalApprox->reduced); 
    ks_gaussNewtonSolve(multiquation, primalApprox, reduced, innertimesolver,
			&iterationCount); 
    if (i == 1 && primalApprox->self->time_method == 2) //Use BDF1 at the be
      innertimesolver = 2;
    else if (i == 1 && primalApprox->self->time_method == 3)
      innertimesolver = 2;
    else if (i == 2 && primalApprox->self->time_method == 3)
      innertimesolver = 3; //After two iterations you can now go and use BDF3
    //-------------------------------------------------------------------------
    // Print out the reduced states found at everytime time iteration
    //-------------------------------------------------------------------------
    ks_printSolution(_equationReduced, primalApprox->reduced, i);
    ks_printSolution(_equationFull, primalApprox->self, i);
   //--------------------------------------------------------------------------
    // Real solution discovery with the regular V (for verification)
    //-------------------------------------------------------------------------
    array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec);
    MatMult(reduced->rOBState, primalReducedVec, primalRealVec);
    vec2Array(_equationFull, primalApprox->stateFull, primalRealVec);
    ks_printSolution(_equationFull, primalApprox->stateFull, i);
  }  
  VecDestroy(&primal_0);
  VecDestroy(&primalApprox_0);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalRealVec);
  VecDestroy(&primalApproxVec);
  MatDestroy(&sResMat);
  MatDestroy(&sJacMat);
  MatDestroy(&reduced->A);
  MatDestroy(&reduced->B);
  free(nodeSet);
}


void minElements(Galaxy *primal, PetscInt *nodeSet, Is_it *reduced){
  int i, j;
  int elem;
  int systemSize = primal->space.node.count;
  int *tempNodeUSet = (int *) malloc ((systemSize)*sizeof(int));
  int *nodeMain = (int *) malloc (systemSize*sizeof(int));
  int countReducedMesh = 0;
  int sum = 0;
  int sumMain = 0;
  for (i=0; i<systemSize; i++){
    tempNodeUSet[i] = 0;
    nodeMain[i] = 0;
  }
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->nSampleNodes; i++){
    elem = floor(nodeSet[i]/(primal->basis.p+1));
    tempNodeUSet[elem] = 1;
    nodeMain[elem] = 1;
    if (elem > 0)
      tempNodeUSet[elem-1] = 1;
    if (elem<primal->space.elem.array[primal->space.elem.count-1])
      tempNodeUSet[elem+1] = 1;
  }
  for (i=0; i<systemSize; i++){
    sum += tempNodeUSet[i];
    sumMain += nodeMain[i];
  }
  reduced->reducedMesh.dx = primal->space.dx;
  initArrayInt(&reduced->reducedMesh.elem, sumMain);
  initArrayInt(&reduced->reducedMesh.node, (primal->basis.p+1)*sumMain);
  initArray(&reduced->reducedMesh.x, (primal->basis.p+1)*sumMain);
  reduced->reducedMesh.x_0 = primal->space.x_0;
  reduced->reducedMesh.x_f = primal->space.x_f;
  for (i=0; i<systemSize; i++){
    if (nodeMain[i] == 1){
      reduced->reducedMesh.elem.array[countReducedMesh] = i;
      for (j=0; j<primal->basis.p+1; j++){
        reduced->reducedMesh.node.array
	  [countReducedMesh*(primal->basis.p+1)+j]= i*(primal->basis.p+1)+j;
        reduced->reducedMesh.x.array[countReducedMesh*(primal->basis.p+1)+j] =
          i*primal->space.dx +(xiL(primal->basis.p, j)+1)/2.0*primal->space.dx;
      }
      countReducedMesh++;
    }
  }
  free(tempNodeUSet);
  free(nodeMain);
}
  
  
/* void findReducedApproxMesh(Galaxy *primal, Mesh *reduced, int **tempM, int *sum){ */
/*   int i; */
/*   int sumD = *sum; */
/*   int *nodeSet = (int *) malloc (primal->space.node.count*sizeof(int)); */
/*   sumD = 0; */
/*   int count = 0; */
/*   for (i=0; i<primal->space.node.count; i++) */
/*     nodeSet[i] = 0; //Zero everything out here */
/*   for (i=0; i<reduced->elem.count; i++){ */
/*     nodeSet[reduced->elem.array[i]*3] = 1; */
/*     nodeSet[reduced->elem.array[i]*3 +1] = 1; */
/*     nodeSet[reduced->elem.array[i]*3 +2] = 1; */
/*     if (reduced->elem.array[i]>0){ */
/*       nodeSet[(reduced->elem.array[i]-1)*3] = 1; */
/*       nodeSet[(reduced->elem.array[i]-1)*3+1] = 1; */
/*       nodeSet[(reduced->elem.array[i]-1)*3+2] = 1;          */
/*     } */
/*     if (reduced->elem.array[i]< */
/*         primal->space.elem.array[primal->space.elem.count-1]){ */
/*       nodeSet[(reduced->elem.array[i]+1)*3] = 1; */
/*       nodeSet[(reduced->elem.array[i]+1)*3+1] = 1; */
/*       nodeSet[(reduced->elem.array[i]+1)*3+2] = 1; */
/*     } */
/*   } */
 
/*   //nodeSet[120] = 0; */
/*   /\* for (i=0; i<primal->space.node.count; i++) *\/ */
/*   /\*   nodeSet[i] = 1; *\/ */
/*   for (i=0; i<primal->space.node.count; i++) */
/*     sumD += nodeSet[i]; */
/*   *sum = sumD; */
/*   *tempM = (int *) malloc (*sum*sizeof(int)); */
/*   for (i=0; i<primal->space.node.count; i++){ */
/*     if (nodeSet[i] == 1){ */
/*       (*tempM)[count] = i; */
/*       count++; */
/*     } */
/*   } */
/* } */

void findApproxMesh(Galaxy *primal, Galaxy *primalApprox, Is_it *reduced){
  int i, j;
  int sum = 0;
  int count = 0;
  int *nodeSet = (int *) malloc (primal->space.elem.count*sizeof(int));
  //---------------------------------------------------------------------------
  // Count number of total elements in the approx mesh
  //---------------------------------------------------------------------------
  for (i=0; i<primal->space.elem.count; i++)
    nodeSet[i] = 0;
  for (i=0; i<reduced->reducedMesh.elem.count; i++){
    nodeSet[reduced->reducedMesh.elem.array[i]] = 1;
    if (reduced->reducedMesh.elem.array[i] > 0)
      nodeSet[reduced->reducedMesh.elem.array[i]-1] = 1;
    if (reduced->reducedMesh.elem.array[i] <
	primal->space.elem.array[primal->space.elem.count-1])
      nodeSet[reduced->reducedMesh.elem.array[i]+1] = 1;
  }
  //---------------------------------------------------------------------------
  // Now create that approximate mesh
  //---------------------------------------------------------------------------
  for (i=0; i<primal->space.elem.count; i++)
    sum += nodeSet[i]; //Number of elemenets
  initArrayInt(&primalApprox->space.elem, sum);
  initArrayInt(&primalApprox->space.node, sum*(primal->basis.p+1));
  initArray(&primalApprox->space.x, sum*(primal->basis.p+1));
  for (i=0; i<primal->space.elem.count; i++){
    if (nodeSet[i] == 1){
      primalApprox->space.elem.array[count] = i;
      count ++;
    }
  }
  for (i=0; i<sum; i++){
    for (j=0; j<primal->basis.p+1; j++){
      primalApprox->space.node.array[i*(primal->basis.p+1)+j] =
        primalApprox->space.elem.array[i]*(primal->basis.p+1)+j;
      primalApprox->space.x.array[i*(primal->basis.p+1)+j]=
        primalApprox->space.elem.array[i]*primal->space.dx +
        (xiL(primal->basis.p, j)+1)/2.0*primal->space.dx;
    }
  }
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  free(nodeSet);
}
