#include "yk_gaussNewtonSolve.h"
#include <time.h>

void yk_gaussNewtonSolve_ST(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primalApprox, Is_it *reduced, int tSolve,
			    int *iterationCount){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             //initialization for iteration
  int t_i;
  int alfa = 1, flag = 0, st_systemSize;
  int basisNodes =primalApprox->self->basis.nodes;
  int eqnStates = multiquation->equation.numStates;
  int elemSize = eqnStates*basisNodes;
  int systemSize = primalApprox->self->systemSize;
  int gaussCount = 0;
  int innertimesolver = 1;
  int left_t_node, right_t_node;
  double dt = primalApprox->self->time.dt;
  double cpu_time_used;
  double *snapshotsRJ;
  Universe _eqn = multiquation->equation;
  Universe _eqnRd = multiquation->equationReduced;
  clock_t start, end;
  FILE *residualSnapshotFile;
  PetscInt size_Os, size_Xs;
  PetscInt rowRes, colRes;
  PetscInt rowB;
  PetscInt bckct_i, local_i;
  PetscInt we = 0;
  PetscInt basis;
  PetscInt *rowIndex, *index_Xs, *nSTBasisIndex;
  PetscScalar *residual_i, *cc;
  PetscScalar normp=1.0;
  PetscReal val;
  Vec p;
  Vec state_0, state_0_final;
  Vec residualTemp;
  Vec primalReducedVec;
  Vec primalApproxVec_i, primalApproxVec_i_final;
  Vec st_residual, r_st_residual;
  Vec *ptr_residual = NULL;
  Mat st_rOBState = reduced->ST_rOBState;
  Mat dRdU=NULL, dRdUTemp = NULL;
  Mat C, r_dRdU;
  Mat coefMat, coefMat1, coefMat2;
  Mat ST_rOBState = NULL, *ptr_dRdU= NULL;
  Mat *ST_rOBState_i;
  Mesh reducedSpace = primalApprox->reduced->space;
  Mesh approxSpace = primalApprox->self->space;
  Mesh mesh_Os;
  PetscInt mm, nn;
  PetscScalar test;
  char residualSnapshot[1000];
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatGetSize(reduced->ST_rOBState, &rowB, &basis);
  left_t_node = primalApprox->self->time.t_0/dt;
  right_t_node = primalApprox->self->time.t_f/dt;

  if (reduced->hrom==0){
    ST_rOBState = reduced->ST_rOBState;
    ST_rOBState_i = reduced->ST_rOBState_i;
    size_Os = systemSize;
    size_Xs = systemSize;
    mesh_Os = primalApprox->self->space;
    st_systemSize = size_Xs*primalApprox->self->time.count;
    snapshotsRJ = (double *) malloc (st_systemSize*sizeof(double));
  }else if (reduced->hrom==1){
    MatGetSize(reduced->ST_rOBResidual, &rowRes, &colRes);
    ST_rOBState = reduced->ST_rOBStateBar;
    ST_rOBState_i = reduced->ST_rOBStateBar_i;
    size_Os = reduced->reducedMesh_Os.elem.count*elemSize;
    size_Xs = reduced->reducedMesh.elem.count*elemSize;
    mesh_Os = reduced->reducedMesh_Os;
    st_systemSize = size_Xs*(reduced->nSampleTime-1);
  }
  rowIndex = (PetscInt *) malloc (st_systemSize*sizeof(PetscInt));
  nSTBasisIndex = (PetscInt *) malloc (basis*sizeof(PetscInt));
  cc = (PetscScalar *) malloc (size_Xs*basis*sizeof(PetscScalar));
  /* residual_i = (PetscScalar *) malloc (size_Xs*sizeof(PetscScalar)); */
  index_Xs = (PetscInt *) malloc(size_Xs*sizeof(PetscInt));

  for (i=0; i<size_Xs; i++) {index_Xs[i] = i;}

  for (i=0; i<basis; i++)  {nSTBasisIndex[i] = i;}

  for (i=0; i<st_systemSize; i++) {rowIndex[i] = i;}

  VecCreate(PETSC_COMM_SELF, &p);
  VecSetSizes(p, basis, basis);
  VecSetType(p, VECSEQ);

  printf("%d\n", size_Xs);
  printf("guass newton\n");
  VecCreate(PETSC_COMM_SELF, &residualTemp);
  VecSetSizes(residualTemp, size_Xs, size_Xs);
  VecSetType(residualTemp, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, basis, basis);
  VecSetType(primalReducedVec, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &state_0_final);
  VecSetSizes(state_0_final, systemSize, systemSize);
  VecSetType(state_0_final, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &state_0);
  VecSetSizes(state_0, size_Os, size_Os);
  VecSetType(state_0, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &st_residual);
  VecSetSizes(st_residual, PETSC_DECIDE, st_systemSize);
  VecSetBlockSize(st_residual, size_Xs);
  VecSetFromOptions(st_residual);

  if (reduced->hrom==1){
    VecCreate(PETSC_COMM_SELF, &r_st_residual);
    VecSetSizes(r_st_residual, colRes, colRes);
    VecSetType(r_st_residual, VECSEQ);
  }

  VecCreate(PETSC_COMM_SELF, &primalApproxVec_i);
  VecSetSizes(primalApproxVec_i, size_Os, size_Os);
  VecSetType(primalApproxVec_i, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &primalApproxVec_i_final);
  VecSetSizes(primalApproxVec_i_final, systemSize, systemSize);
  VecSetType(primalApproxVec_i_final, VECSEQ);

  MatCreate(PETSC_COMM_SELF, &dRdU);
  MatSetSizes(dRdU, PETSC_DECIDE, basis, st_systemSize, basis);
  MatSetBlockSizes(dRdU, size_Xs, basis);
  MatSetType(dRdU, MATSEQDENSE);
  MatSetFromOptions(dRdU);
  MatSetUp(dRdU);

  MatCreate(PETSC_COMM_SELF, &dRdUTemp);
  MatSetSizes(dRdUTemp, size_Xs, size_Os, size_Xs, size_Os);
  MatSetType(dRdUTemp, MATSEQBAIJ);
  MatSeqBAIJSetPreallocation(dRdUTemp, elemSize ,ykflow->numElemCol, NULL);

  MatCreate(PETSC_COMM_SELF, &C);
  MatSetSizes(C, size_Xs, basis, size_Xs, basis);
  MatSetType(C, MATSEQDENSE);
  MatSetUp(C);
  //---------------------------------------------------------------------------
  // Need to build up the Residual and the Jaociban vector and matrix
  //---------------------------------------------------------------------------
  if (reduced->hrom == 0){
    sprintf(residualSnapshot, "%s_%d.dat", "residualSnapshot",
	    primalApprox->self->time.window_i);
    residualSnapshotFile = fopen(residualSnapshot, "w+");

  }
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Rertrieve the initial conditions for the full order model solution
  ks_readSolution(_eqnRd, primalApprox->reduced, primalApprox->self->time.window_i);

  array2Vec(_eqnRd, primalApprox->reduced, reducedSpace, primalReducedVec);
  VecView(primalReducedVec, PETSC_VIEWER_STDOUT_SELF);
  //initial conditions (I don't think I need this
  if (reduced->hrom==0){
    ks_readSolution(_eqn, primalApprox->self, left_t_node);
    array2Vec(_eqn, primalApprox->self, mesh_Os, state_0);
  }else if (reduced->hrom == 1){
    ks_readSolution(_eqn, primalApprox->stateFull, left_t_node);
    array2Vec(_eqn, primalApprox->stateFull, mesh_Os, state_0);
  }
  MatCreate(PETSC_COMM_SELF, &coefMat);
  MatSetSizes(coefMat, size_Xs, size_Os, size_Xs, size_Os);
  MatSetType(coefMat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(coefMat, elemSize*elemSize, NULL);
  ykflow->findMassCoefBDF(ykflow, primalApprox->self, &coefMat, reduced);
  MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat1);
  MatScale(coefMat1, -2/0.1);
  MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat2);
  MatScale(coefMat2, 0.5/0.1);
  ptr_dRdU = &dRdU;
  ptr_residual = &st_residual;

  do{
    MatZeroEntries(dRdU);
    //-------------------------------------------------------------------------
    // Precalculate all the states at each time for the given reduced states
    //-------------------------------------------------------------------------
    for (t_i = 0; t_i<primalApprox->self->time.time_Os.count; t_i++){
      i = left_t_node+primalApprox->self->time.time_Os.array[t_i]+1;
      local_i = primalApprox->self->time.time_Os.array[t_i]+1;
      primalApprox->self->time.node = i;
      MatMultAdd(ST_rOBState_i[local_i-1], primalReducedVec, state_0,
                 primalApproxVec_i);
      vec2Array(_eqn, primalApprox->self, mesh_Os, primalApproxVec_i);
      VecNorm(primalApproxVec_i, NORM_2, &test);
      ks_printSolution(_eqn, primalApprox->self, i);
    }
    //-------------------------------------------------------------------------
    // Loop through each time to obtain the residual and Jacobian at tstep i
    //-------------------------------------------------------------------------
    for (t_i = 0; t_i<primalApprox->self->time.tNode_i.count; t_i++){
      i = left_t_node+primalApprox->self->time.tNode_i.array[t_i]+1;
      local_i = primalApprox->self->time.tNode_i.array[t_i]+1;
      ks_readSolution(_eqn, primalApprox->self, i);
      //-----------------------------------------------------------------------
      //Calculate the Residual and the Jacobian here
      //-----------------------------------------------------------------------
      primalApprox->self->time.node = i;
      VecZeroEntries(residualTemp);
      ykflow->Residual(ykflow, multiquation, primalApprox->self, residualTemp,
      		       dRdUTemp, reduced, innertimesolver);
      VecNorm(residualTemp, NORM_2, &val);
      //-----------------------------------------------------------------------
      //Add the residual for time i to a large ass space time vector
      //-----------------------------------------------------------------------
      VecGetArray(residualTemp, &residual_i);

      bckct_i = t_i;
      VecSetValuesBlocked(st_residual, 1, &bckct_i, residual_i,INSERT_VALUES);
      VecRestoreArray(residualTemp, &residual_i);
      //-----------------------------------------------------------------------
      // Jacobain manipulation
      //-----------------------------------------------------------------------
      MatMatMultBAIJ(dRdUTemp, ST_rOBState_i[local_i-1], C);
      MatGetValues(C, size_Xs, index_Xs, basis, nSTBasisIndex, cc);
      MatSetValuesBlocked(dRdU, 1, &bckct_i,1, &we, cc,ADD_VALUES);

      if (local_i>1){
    	MatMatMultBAIJ(coefMat1, ST_rOBState_i[local_i-2], C);
    	MatGetValues(C, size_Xs, index_Xs, basis, nSTBasisIndex, cc);
    	MatSetValuesBlocked(dRdU, 1, &bckct_i,1, &we, cc,ADD_VALUES);
      }

      if (local_i>2){ //BDF 2 section
      	MatMatMultBAIJ(coefMat2, ST_rOBState_i[local_i-3], C);
      	MatGetValues(C, size_Xs, index_Xs, basis, nSTBasisIndex, cc);
      	MatSetValuesBlocked(dRdU, 1, &bckct_i, 1, &we, cc, ADD_VALUES);
      }
    }
    MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(st_residual);
    VecAssemblyEnd(st_residual);
    //-------------------------------------------------------------------------
    // if hrom = 0 (just reduced order model)
    //-------------------------------------------------------------------------
    if (reduced->hrom ==0){
      VecGetValues(st_residual, st_systemSize, rowIndex, snapshotsRJ);
      fprintf(residualSnapshotFile, "%0.20e", snapshotsRJ[0]);
      for (j=1; j<st_systemSize; j++)
    	fprintf(residualSnapshotFile, " %0.20e", snapshotsRJ[j]);
      fprintf(residualSnapshotFile, "\n");
    }
    if (normp >pow(10,-6)){
      VecScale(st_residual, -1);
      //---------------------------------------------------------------------
      // Least Squares Solve for the Direction p
      //-----------------------------------------------------------------------
      if (reduced->hrom==1){
	MatMatMult(reduced->A, dRdU, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &r_dRdU);
	MatMult(reduced->A, st_residual, r_st_residual);
	linearLeastSquares(r_dRdU, r_st_residual, p);
      }else if (reduced->hrom==0)
	linearLeastSquares(dRdU, st_residual, p);
      //-----------------------------------------------------------------------
      // Perform line search or set \alpha to 1
      //-----------------------------------------------------------------------
      gaussCount++;
      VecNorm(p, 2, &normp);
      printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp);
      /* flag = normp <pow(10, -6) ? 1:  0; */
      VecAXPY(primalReducedVec, alfa, p);
      vec2Array(_eqnRd, primalApprox->reduced,  reducedSpace, primalReducedVec);


    }else
      flag = 1;

    if (reduced->hrom==1){
      MatDestroy(&r_dRdU);
    }
  }while (flag == 0);

  //---------------------------------------------------------------------------
  //Print the reduced solution based on the time window number
  //---------------------------------------------------------------------------
  ks_printSolution(_eqnRd, primalApprox->reduced,
		   primalApprox->self->time.window_i);

  if (reduced->hrom==0){
    ks_readSolution(_eqn, primalApprox->self, left_t_node);
    array2Vec(_eqn, primalApprox->self, primalApprox->self->space, state_0_final);
  }else if (reduced->hrom == 1){
    ks_readSolution(_eqn, primalApprox->stateFull, left_t_node);
    array2Vec(_eqn, primalApprox->stateFull, primalApprox->stateFull->space, state_0_final);
  }
  //Should print out the approximated entire final solution too

  //---------------------------------------------------------------------------
  // Printf
  //---------------------------------------------------------------------------
  for (i=left_t_node+1; i<right_t_node+1; i++){
    local_i = i-left_t_node;
    MatMultAdd(reduced->ST_rOBState_i[local_i-1], primalReducedVec,
  	       state_0_final, primalApproxVec_i_final);
    if (reduced->hrom==1){
      vec2Array(_eqn, primalApprox->stateFull, primalApprox->stateFull->space,
  		primalApproxVec_i_final);
      ks_printSolution(_eqn, primalApprox->stateFull, i);
    }else if (reduced->hrom==0){
      vec2Array(_eqn, primalApprox->self, approxSpace,primalApproxVec_i_final);
      ks_printSolution(_eqn, primalApprox->self, i);
    }
  }
  //-------------------------------------------------------------------------
  // Print out the reduced states found at everytime time iteration
  //-------------------------------------------------------------------------
  VecDestroy(&p);
  VecDestroy(&residualTemp);
  VecDestroy(&primalReducedVec);
  VecDestroy(&state_0);
  VecDestroy(&state_0_final);
  VecDestroy(&st_residual);
  VecDestroy(&primalApproxVec_i);
  VecDestroy(&primalApproxVec_i_final);
  if (reduced->hrom ==1)
    VecDestroy(&r_st_residual);

  MatDestroy(&dRdU);
  MatDestroy(&dRdUTemp);
  MatDestroy(&C);
  MatDestroy(&coefMat);
  MatDestroy(&coefMat1);
  MatDestroy(&coefMat2);
  if (reduced->hrom == 0){
    fclose(residualSnapshotFile);
    free(snapshotsRJ);
  }
  free(index_Xs);
  free(nSTBasisIndex);
  free(rowIndex);
  /* free(residual_i); */
  free(cc);
}

void yk_gaussNewtonSolve(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Cluster *primalApprox, Is_it *reduced, int tSolve,
			 int *iterationCount){
/*   int i, j; */
/*    clock_t start, end; */
/*   int flag = 0; */
/*   int innertimesolver = 1; */
/*   int meshSize; */
/*   int gaussCount = 0; */
/*   int systemSize; */
/*   int numNodes; */
/*   int *index; */
/*   int nBasisFuncs = reduced->nBasisFuncs; */
/*   double cpu_time_used; */
/*   double *snapshotsRJ; */
/*   double alfa; */
/*   FILE *residualSnapshotFile; */
/*   FILE *jacobianSnapshotFile; */
/*   PetscScalar normp; */
/*   Vec p; */
/*   Vec residual = NULL, residualTemp = NULL, residualRed = NULL; */
/*   Vec primalReducedVec, primalApproxVec, snapJacobian; */
/*   Mat dRdU, dRdUTemp = NULL, dRdUTempHat = NULL, dRdUbasis = NULL; */
/*   Vec primalRealVec; */
/*   Vec state_0; */
/*   Mesh _meshOfInterest; */
/*   Mat rOBState = NULL; */
/*   Universe _equation = multiquation->equation; */
/*   Universe _equationReduced = multiquation->equationReduced; */
/*   //--------------------------------------------------------------------------- */
/*   // Mallocing things occur here/Initialize things here */
/*   //--------------------------------------------------------------------------- */
/*   systemSize = _equation.numStates*primalApprox->self->space.node.count; */
/*   index = (int *) malloc (systemSize*sizeof(int)); */
/*   for (i=0; i<systemSize; i++) */
/*     index[i] = i; */
/*   if (reduced->hrom == 1){ */
/*     numNodes = reduced->nSampleNodes; */
/*     rOBState = reduced->rOBStateBar; */
/*     _meshOfInterest = reduced->reducedMesh; */
/*     meshSize = _equation.numStates*reduced->reducedMesh.node.count; */
/*     VecCreate(PETSC_COMM_SELF, &residualRed); */
/*     VecSetSizes(residualRed, numNodes, numNodes); */
/*     VecSetType(residualRed, VECSEQ); */
/*       //primalRealVec */
/*     //VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalRealVec); */
/*     VecCreate(PETSC_COMM_SELF, &primalRealVec); */
/*     VecSetSizes(primalRealVec, primalApprox->stateFull->systemSize, */
/* 		primalApprox->stateFull->systemSize); */
/*     VecSetType(primalRealVec, VECSEQ); */
/*     //dRdUbasis */
/*     MatCreate(PETSC_COMM_SELF, &dRdUbasis); */
/*     MatSetSizes(dRdUbasis, meshSize, nBasisFuncs, meshSize, nBasisFuncs); */
/*     MatSetType(dRdUbasis, MATSEQDENSE); */
/*     MatSetUp(dRdUbasis); */
/*     //dRdUTempHat */
/*     MatCreate(PETSC_COMM_SELF, &dRdUTempHat); */
/*     MatSetSizes(dRdUTempHat, numNodes, reduced->nBasisFuncs, numNodes, */
/* 		reduced->nBasisFuncs); */
/*     MatSetType(dRdUTempHat, MATSEQDENSE); */
/*     MatSetUp(dRdUTempHat); */
/*     //dRdU */
/*   }else{ */
/*     numNodes = primalApprox->self->systemSize; */
/*     residualSnapshotFile = fopen("residualSnapshot.dat", "a"); */
/*     jacobianSnapshotFile = fopen("jacobianSnapshot.dat", "a"); */
/*     snapshotsRJ = (double *) malloc (numNodes*sizeof(double)); */
/*     rOBState = reduced->rOBState; */
/*     meshSize = _equation.numStates*primalApprox->self->space.node.count; */
/*     _meshOfInterest = primalApprox->self->space; */
/*     //snapJacobian */
/*     VecCreate(PETSC_COMM_SELF, &snapJacobian); */
/*     VecSetSizes(snapJacobian, numNodes, numNodes); */
/*     VecSetType(snapJacobian, VECSEQ); */
/*   } */


/*   VecCreate(PETSC_COMM_SELF, &state_0); */
/*   VecSetSizes(state_0, systemSize, systemSize); */
/*   VecSetType(state_0, VECSEQ); */


/*   VecCreate(PETSC_COMM_SELF, &residual); */
/*   VecSetSizes(residual, numNodes, numNodes); */
/*   VecSetType(residual, VECSEQ); */
/*   //residualTemp */
/*   VecCreate(PETSC_COMM_SELF, &residualTemp); */
/*   VecSetSizes(residualTemp, meshSize, meshSize); */
/*   VecSetType(residualTemp, VECSEQ); */
/*   //p direction */
/*   VecCreate(PETSC_COMM_SELF, &p); */
/*   VecSetSizes(p, nBasisFuncs, nBasisFuncs); */
/*   VecSetType(p, VECSEQ); */
/*   //primalReducedVec */
/*   VecCreate(PETSC_COMM_SELF, &primalReducedVec); */
/*   VecSetSizes(primalReducedVec, nBasisFuncs, nBasisFuncs); */
/*   VecSetType(primalReducedVec, VECSEQ); */
/*   //primalApproxVec */
/*   VecCreate(PETSC_COMM_SELF, &primalApproxVec); */
/*   VecSetSizes(primalApproxVec, systemSize, systemSize); */
/*   VecSetType(primalApproxVec, VECSEQ); */
/*   //dRdU */
/*   MatCreate(PETSC_COMM_SELF, &dRdU); */
/*   MatSetSizes(dRdU, numNodes, nBasisFuncs, numNodes, nBasisFuncs); */
/*   MatSetType(dRdU, MATSEQDENSE); */
/*   MatSetUp(dRdU); */
/*   //dRdUTemp */
/*   MatCreate(PETSC_COMM_SELF, &dRdUTemp); */
/*   MatSetSizes(dRdUTemp, meshSize, systemSize, meshSize, systemSize); */
/*   MatSetType(dRdUTemp, MATSEQBAIJ); */
/*   MatSeqBAIJSetPreallocation(dRdUTemp, primalApprox->self->basis.nodes* */
/* 			     multiquation->equation.numStates, */
/* 			     ykflow->numElemCol, NULL);  //Need to think about */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   //Rertrieve the initial conditions for the full order model solution */
/*   ks_readSolution(_equation, primalApprox->self, 0); */
/*   array2Vec(_equation, primalApprox->self, primalApprox->self->space, state_0); */
/*   //Run for each time step */
/*   for (i = 1; i<primalApprox->self->time.count+1; i++){ */
/*   //  for(i=1; i<2; i++){ */
/*     flag = 0; */
/*     gaussCount=0; */

/*     ks_readSolution(_equationReduced, primalApprox->reduced, i-1); */
/*     //    ks_setSolutionZero(_equationReduced, primalApprox->reduced); */
/*     array2Vec(_equationReduced, primalApprox->reduced, primalApprox->reduced->space, primalReducedVec); */
/*     primalApprox->self->time.node = i; */
/*     do{ */
/*       //----------------------------------------------------------------------- */
/*       // Extract the guesses made from the previous solution */
/*       //----------------------------------------------------------------------- */
/*       //Calculate the approximate solution given initial guess */
/*       MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec); */

/*       //MatMult(rOBState, primalReducedVec, primalApproxVec); */
/*       vec2Array(_equation, primalApprox->self, primalApprox->self->space, */
/* 		primalApproxVec); */
/*       //----------------------------------------------------------------------- */
/*       //Calculate the Residual and the Jacobian here */
/*       //----------------------------------------------------------------------- */
/*       ykflow->Residual(ykflow, multiquation, primalApprox, residualTemp, */
/* 		       dRdUTemp, reduced, innertimesolver); */
/*       //----------------------------------------------------------------------- */
/*       //Gauss Newton with Approximated Tensors for hyper-reduced order modeling */
/*       //----------------------------------------------------------------------- */
/*       if (reduced->hrom == 1){ */
/*         MatMult(reduced->Zmult, residualTemp, residualRed); */
/* 	MatMult(reduced->B, residualRed, residual); */
/* 	MatMatMultBAIJ(dRdUTemp, rOBState, dRdUbasis); */
/* 	MatMatMult(reduced->Zmult, dRdUbasis, MAT_REUSE_MATRIX, PETSC_DEFAULT, */
/* 		   &dRdUTempHat); */
/* 	MatMatMult(reduced->A, dRdUTempHat, MAT_REUSE_MATRIX, PETSC_DEFAULT, */
/* 		   &dRdU); */
/*       }else{ */
/* 	VecCopy(residualTemp, residual); */
/* 	MatMatMultBAIJ(dRdUTemp, rOBState, dRdU); */
/*       } */
/*       //----------------------------------------------------------------------- */
/*       // Residual Snapshots */
/*       //----------------------------------------------------------------------- */
/*       if (reduced->hrom ==0){ //If we are doing regular rom */
/* 	VecGetValues(residual, systemSize, index, snapshotsRJ); */
/* 	fprintf(residualSnapshotFile, "%0.16f", snapshotsRJ[0]); */
/* 	for (j=1; j<systemSize; j++) */
/* 	  fprintf(residualSnapshotFile, " %0.16f", snapshotsRJ[j]); */
/* 	fprintf(residualSnapshotFile, "\n"); */
/*       } */
/*       //----------------------------------------------------------------------- */
/*       // Least Squares Solve for the Direction p */
/*       //----------------------------------------------------------------------- */
/*       VecScale(residual, -1); */
/*       start = clock(); */
/*       linearLeastSquares(dRdU, residual, p); */
/*       end = clock(); */
/*       //cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC; */
/*       //      printf("%g\n", cpu_time_used); */
/*       (*iterationCount)++; */
/*       //----------------------------------------------------------------------- */
/*       // Jacobian Snapshots */
/*       //----------------------------------------------------------------------- */
/*       if (reduced->hrom == 0){ */
/* 	MatMult(dRdU, p, snapJacobian); */
/* 	VecGetValues(snapJacobian, systemSize, index, snapshotsRJ); */
/* 	fprintf(jacobianSnapshotFile, "%0.16f", snapshotsRJ[0]); */
/* 	for (j=1; j<systemSize; j++) */
/* 	  fprintf(jacobianSnapshotFile, " %0.16f", snapshotsRJ[j]); */
/* 	fprintf(jacobianSnapshotFile, "\n"); */
/*       } */
/*       //----------------------------------------------------------------------- */
/*       // Perform line search or set \alpha to 1 */
/*       //----------------------------------------------------------------------- */
/*       alfa = 1; */
/*       VecNorm(p, 2, &normp); */
/*       printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp); */
/*       gaussCount++; */
/*       if (normp <pow(10, -8)) */
/* 	flag = 1; */
/*       VecAXPY(primalReducedVec, alfa, p); */
/*       vec2Array(_equationReduced, primalApprox->reduced, */
/* 		primalApprox->reduced->space, primalReducedVec); */
/*     }while (flag ==0); */

/*     MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec); */
/*     //MatMult(rOBState, primalReducedVec, primalApproxVec); */
/*     vec2Array(_equation, primalApprox->self, primalApprox->self->space, primalApproxVec); */
/*     vec2Array(_equationReduced, primalApprox->reduced, primalApprox->self->space, primalReducedVec); */
/*     ks_printSolution(_equationReduced, primalApprox->reduced, i); */
/*     ks_printSolution(_equation, primalApprox->self, i); */
/*     if (reduced->hrom == 1){ */
/*       array2Vec(_equationReduced, primalApprox->reduced, primalApprox->reduced->space, primalReducedVec); */
/*       MatMult(reduced->rOBState, primalReducedVec, primalRealVec); */
/*       vec2Array(_equation, primalApprox->stateFull, primalApprox->stateFull->space,		primalRealVec); */
/*       ks_printSolution(_equation, primalApprox->stateFull, i); */
/*     } */
/*     if (i == 1 && primalApprox->self->time_method == 2) //Use BDF1 */
/*       innertimesolver = 2; */
/*     else if (i == 1 && primalApprox->self->time_method == 3) */
/*       innertimesolver = 2; */
/*     else if (i == 2 && primalApprox->self->time_method == 3) */
/*       innertimesolver = 3; //After two iterations you can now go and use BDF3 */
/*     //------------------------------------------------------------------------- */
/*     // Print out the reduced states found at everytime time iteration */
/*     //------------------------------------------------------------------------- */
/*   } */
/*   //--------------------------------------------------------------------------- */
/*   // Destroy all the malloc stuff and initialization stuff OMG */
/*   //--------------------------------------------------------------------------- */
/*   VecDestroy(&p); */
/*   VecDestroy(&primalReducedVec); */
/*   VecDestroy(&primalApproxVec); */
/*   VecDestroy(&residual); */
/*   VecDestroy(&state_0); */
/*   VecDestroy(&residualTemp); */
/*   MatDestroy(&dRdU); */
/*   MatDestroy(&dRdUTemp); */
/*   if (reduced->hrom == 1){ */
/*     VecDestroy(&residualRed); */
/*     VecDestroy(&primalRealVec); */
/*     MatDestroy(&dRdUTempHat); */
/*     MatDestroy(&dRdUbasis); */
/*   }else if (reduced->hrom == 0){ */
/*     fclose(residualSnapshotFile); */
/*     fclose(jacobianSnapshotFile); */
/*     VecDestroy(&snapJacobian); */
/*     free(snapshotsRJ); */
/*   } */
/*   free(index); */
}
