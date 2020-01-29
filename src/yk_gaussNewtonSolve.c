#include "yk_gaussNewtonSolve.h"
#include <time.h>


void yk_gaussNewtonSolve_ST(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primalApprox, Is_it *reduced, int tSolve,
			    int *iterationCount){
  int i, j;                          //initialization for iteration
  double *state_i = (double *) malloc      // the current full state and time i
    (primalApprox->self->systemSize*sizeof(double));
  PetscInt *spaceIndex = (PetscInt *) malloc (primalApprox->self->systemSize*
					      sizeof(PetscInt ));
  PetscInt *residualIndex = (PetscInt *) malloc (primalApprox->self->
						 systemSize*sizeof(PetscInt));

  PetscInt *index = (PetscInt *) malloc (primalApprox->self->systemSize*
					 sizeof(PetscInt ));
  for (i=0; i<primalApprox->self->systemSize; i++)
    index[i] = i;
  PetscScalar *residual_i;
  /* double *residual_i= (double *) malloc */
  /*   (primalApprox->self->systemSize*sizeof(double)); */
  PetscInt local_i;
  int gaussCount = 0;
  int innertimesolver = 1;
  Vec state_0;
  Vec residualTemp;
  Vec primalReducedVec;
  Vec primalApproxVec;
  Vec primalApproxVec_i;
  Vec st_residual;
  Mat dRdU=NULL;
  Mat dRdUTemp = NULL;
  PetscInt row, col;
  Universe _equation = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  Mat st_rOBState = reduced->ST_rOBState;
  int meshSize = _equation.numStates*primalApprox->self->space.node.count;
  Mesh  _meshOfInterest = primalApprox->self->space;
  PetscInt rowBasis, basis;
  MatGetSize(reduced->ST_rOBState, &row, &basis);
  int nSTBasis = basis;
  int systemSize = primalApprox->self->systemSize;
  int st_systemSize = systemSize*primalApprox->self->time.count;

  PetscInt *nSTBasisIndex = (PetscInt *) malloc (nSTBasis*sizeof(PetscInt));
  for (i=0; i<nSTBasis; i++)
    nSTBasisIndex[i] = i;
  int alfa;
  PetscScalar normp;
  Vec p;
  clock_t start, end;
  double cpu_time_used;

  MPI_Comm comm;
  comm = PETSC_COMM_WORLD;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  VecCreate(PETSC_COMM_SELF, &p);
  VecSetSizes(p, nSTBasis, nSTBasis);
  VecSetType(p, VECSEQ);

  MatGetSize(st_rOBState, &row, &col);
  //residualTemp
  VecCreate(PETSC_COMM_SELF, &residualTemp);
  VecSetSizes(residualTemp, systemSize, systemSize);
  VecSetType(residualTemp, VECSEQ);
  //primalReducedVec
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, nSTBasis, nSTBasis);
  VecSetType(primalReducedVec, VECSEQ);

  VecCreate(PETSC_COMM_SELF, &state_0);
  VecSetSizes(state_0, systemSize, systemSize);
  VecSetType(state_0, VECSEQ);


  VecCreate(PETSC_COMM_SELF, &primalApproxVec);
  VecSetSizes(primalApproxVec, row, row);
  VecSetType(primalApproxVec, VECSEQ);

  VecCreate(comm, &st_residual);
  VecSetSizes(st_residual, PETSC_DECIDE, row);
  VecSetBlockSize(st_residual, systemSize);
  VecSetFromOptions(st_residual);
  /* VecSetType(st_residual, VECSEQ); */

  VecCreate(PETSC_COMM_SELF, &primalApproxVec_i);
  VecSetSizes(primalApproxVec_i, systemSize, systemSize);
  VecSetType(primalApproxVec_i, VECSEQ);

 //dRdU
  MatCreate(comm, &dRdU);
  MatSetSizes(dRdU, PETSC_DECIDE, nSTBasis, st_systemSize, nSTBasis);
  MatSetBlockSizes(dRdU, systemSize, nSTBasis);
  MatSetType(dRdU, MATSEQDENSE);
  MatSetFromOptions(dRdU);
  MatSetUp(dRdU);


  MatCreate(PETSC_COMM_SELF, &dRdUTemp);
  MatSetSizes(dRdUTemp, meshSize, systemSize, meshSize, systemSize);
  MatSetType(dRdUTemp, MATSEQBAIJ);
  MatSeqBAIJSetPreallocation(dRdUTemp, primalApprox->self->basis.nodes*
                             multiquation->equation.numStates,
                             ykflow->numElemCol, NULL);

  /* MatCreate(PETSC_COMM_SELF, &dRdUTemp); */
  /* MatSetSizes(dRdUTemp, meshSize, systemSize, meshSize, systemSize); */
  /* MatSetType(dRdUTemp, MATSEQAIJ); */
  /* MatSeqAIJSetPreallocation(dRdUTemp,  primalApprox->self->basis.nodes* */
  /* 			    multiquation->equation.numStates* */
  /* 			    ykflow->numElemCol, NULL); */


  PetscInt size;
  Mat C;
  int nUnknownsElem = primalApprox->self->basis.nodes*
    multiquation->equation.numStates;
  MatCreate(PETSC_COMM_SELF, &C);
  MatSetSizes(C, systemSize, nSTBasis, systemSize, nSTBasis);
  MatSetType(C, MATSEQDENSE);
  MatSetUp(C);
  PetscInt we = 0;
  Mat coefMat1;
  Mat coefMat2;
  PetscScalar *cc = (PetscScalar *) malloc (systemSize*nSTBasis*sizeof(PetscScalar));
  //---------------------------------------------------------------------------
  // Need to build up the Residual and the Jaociban vector and matrix
  //---------------------------------------------------------------------------
  /* numNodes = primalApprox->self->systemSize; */
  /* residualSnapshotFile = fopen("residualSnapshot.dat", "a"); */
  /* snapshotsRJ = (double *) malloc (numNodes*sizeof(double)); */
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------


  int flag = 0;
  PetscInt blockcount_i;
  Mat coefMat;
  int left_t_node = primalApprox->self->time.t_0/primalApprox->self->time.dt;
  int right_t_node = primalApprox->self->time.t_f/primalApprox->self->time.dt;



  //Rertrieve the initial conditions for the full order model solution
  ks_readSolution(_equationReduced, primalApprox->reduced, left_t_node);
  array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec);

  //initial conditions (I don't think I need this
  ks_readSolution(_equation, primalApprox->self, left_t_node);
  array2Vec(_equation, primalApprox->self, state_0);
  /* /\* Mat *MatCoeff = (Mat *) malloc (2*sizeof(Mat));; //Allocated /Malloced in the *\/ following function: */
  /* printf("%d\n", nUnknownsElem); */
  /* printf("%d\n", systemSize); */
  MatCreate(PETSC_COMM_SELF, &coefMat);
  MatSetSizes(coefMat, systemSize, systemSize, systemSize,
  	      systemSize);
  /* MatSetType(coefMat, MATSEQBAIJ); */
  /* MatSeqBAIJSetPreallocation(coefMat, nUnknownsElem, 1, NULL); */
  MatSetType(coefMat, MATSEQAIJ);
  /* MatSetFromOptions(coefMat); */
  /* MatSetUp(coefMat); */
  MatSeqAIJSetPreallocation(coefMat, primalApprox->self->basis.nodes*4*12, NULL);
  /* printf("watching you\n"); */

  /* primalApprox->self->time.node = 1; */
  /* ykflow->findMassCoefBDF(ykflow, primalApprox->self, &MatCoeff, reduced); */
  /* primalApprox->self->time.node = 1; */
  ykflow->findMassCoefBDF(ykflow, primalApprox->self, &coefMat, reduced);

   MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat1);
   MatScale(coefMat1, -2/0.1);
   //MatScale(coefMat1, -1/0.1);
   MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat2);
   MatScale(coefMat2, 0.5/0.1);
   PetscReal val;
   //  MatZeroEntries(C);
   do{
     //-------------------------------------------------------------------------
     // Loop through each time to obtain the residual and Jacobian at tstep i
    //-------------------------------------------------------------------------
    //Calculate just the Pi*reducedstates part. Add in initial conditions in
    //for loop
    /* MatMult(st_rOBState, primalReducedVec, primalApproxVec); */
    /* printf("will this work hahaha\n"); */
    MatZeroEntries(dRdU);
    //for (i = 1; i<primalApprox->self->time.count+1; i++){
    for (i=left_t_node+1; i<right_t_node+1; i++){
      local_i = i-left_t_node;
      primalApprox->self->time.node = i;
      //temp temp

      //Write to file
      /* MatView(reduced->ST_rOBState_i[local_i-1], PETSC_VIEWER_STDOUT_SELF); */
      /* printf("anymore\n"); */
      /* getchar(); */
      MatMultAdd(reduced->ST_rOBState_i[local_i-1], primalReducedVec, state_0,
      		 primalApproxVec_i);

      /* VecGetValues(primalApproxVec, primalApprox->self->systemSize, residualIndex, */
      /* 		   state_i); */


      /* VecSetValues(primalApproxVec_i, primalApprox->self->systemSize, */
      /* 		   primalApprox->self->index, state_i, INSERT_VALUES); */
      /* VecAXPY(primalApproxVec_i, 1, state_0); */
      vec2Array(_equation, primalApprox->self, primalApproxVec_i);
      ks_printSolution(_equation, primalApprox->self, i);
      /* ks_readSolution(_equation, primalApprox->self, 1); */
      /* printf("working working working\n"); */
      //-----------------------------------------------------------------------
      //Calculate the Residual and the Jacobian here
      //-----------------------------------------------------------------------
      ykflow->Residual(ykflow, multiquation, primalApprox, residualTemp,
      		       dRdUTemp, reduced, innertimesolver);


      VecNorm(residualTemp, 2, &val);

      VecGetSize(residualTemp, &size);
      //-----------------------------------------------------------------------
      //Add the residual for time i to a large ass space time vector
      //-----------------------------------------------------------------------
      VecGetArray(residualTemp, &residual_i);
      blockcount_i = local_i-1;
      VecSetValuesBlocked(st_residual, 1, &blockcount_i, residual_i, INSERT_VALUES);
      //-----------------------------------------------------------------------
      // Jacobain manipulation
      //-----------------------------------------------------------------------
      start = clock();
      //Oh fuck
      MatMatMultBAIJ(dRdUTemp, reduced->ST_rOBState_i[local_i-1], C);
      /* MatMatMult(dRdUTemp, reduced->ST_rOBState_i[local_i-1], MAT_REUSE_MATRIX, */
      /* 		 PETSC_DEFAULT, &C); */
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      MatGetValues(C, systemSize,  primalApprox->self->index, nSTBasis,
		   nSTBasisIndex, cc);

      MatSetValuesBlocked(dRdU, 1, &blockcount_i,1, &we, cc,ADD_VALUES);

      if (local_i>1){

	MatMatMultBAIJ(coefMat1, reduced->ST_rOBState_i[local_i-2], C);
      	MatGetValues(C, systemSize, primalApprox->self->index, nSTBasis,
      		     nSTBasisIndex, cc);
      	MatSetValuesBlocked(dRdU, 1, &blockcount_i,1, &we, cc,ADD_VALUES);

      }


      if (local_i>2){ //BDF 2 section
      	MatMatMultBAIJ(coefMat2, reduced->ST_rOBState_i[local_i-3], C);
      	MatGetValues(C, systemSize, primalApprox->self->index,
      		     nSTBasis, nSTBasisIndex, cc);
      	MatSetValuesBlocked(dRdU, 1, &blockcount_i, 1, &we, cc, ADD_VALUES);
      }
    }
    MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY);

    /* MatViewdRdU, PETSC_VIEWER_STDOUT_SELF); */
    VecAssemblyBegin(st_residual);
    VecAssemblyEnd(st_residual);
    /*     if (reduced->hrom == 1){ */
  /*       MatMult(reduced->Zmult, residualTemp, residualRed); */
  /*     }else{ */
  /* 	VecCopy(residualTemp, residual); */
  /* 	MatMatMultBAIJ(dRdUTemp, rOBState, dRdU); */
  /*     } */

  /*   } */

  /*   //------------------------------------------------------------------------- */
  /*   // Print out the entire Space time Residual Snapshots */
  /*   //------------------------------------------------------------------------- */
  /*   if (reduced->hrom ==0){ //If we are doing regular rom */
  /*     VecGetValues(residual, systemSize, index, snapshotsRJ); */
  /*     fprintf(residualSnapshotFile, "%0.16f", snapshotsRJ[0]); */
  /*     for (j=1; j<systemSize; j++) */
  /* 	fprintf(residualSnapshotFile, " %0.16f", snapshotsRJ[j]); */
  /*     fprintf(residualSnapshotFile, "\n"); */
  /*   } */

  /*   //----------------------------------------------------------------------- */
  /*   // Least Squares Solve for the Direction p */
  /*   //----------------------------------------------------------------------- */
    VecScale(st_residual, -1);
  /*   start = clock(); */
    linearLeastSquares(dRdU, st_residual, p);
    /*   end = clock(); */
  /*   //cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC; */
  /*   //      printf("%g\n", cpu_time_used); */
  /*   (*iterationCount)++; */
  /*   //----------------------------------------------------------------------- */
  /*   // Perform line search or set \alpha to 1 */
  /*   //----------------------------------------------------------------------- */
    alfa = 1;
    gaussCount++;
    VecNorm(p, 2, &normp);
    printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp);
    /* gaussCount++;		/\*  *\/ */
    if (normp <pow(10, -10))
      flag = 1;
    VecAXPY(primalReducedVec, alfa, p);
    vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
  }while (flag ==0);


  /*   MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec); */
  /*   //MatMult(rOBState, primalReducedVec, primalApproxVec); */
  /*   vec2Array(_equation, primalApprox->self, primalApproxVec); */
  /*   vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /*   ks_printSolution(_equationReduced, primalApprox->reduced, i); */
  /*   ks_printSolution(_equation, primalApprox->self, i); */
  /*   if (reduced->hrom == 1){ */
  /*     array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /*     MatMult(reduced->rOBState, primalReducedVec, primalRealVec); */
  /*     vec2Array(_equation, primalApprox->stateFull, primalRealVec); */
  /*     ks_printSolution(_equation, primalApprox->stateFull, i); */
  /*   } */
  /*   if (i == 1 && primalApprox->self->time_method == 2) //Use BDF1 */
  /*     innertimesolver = 2; */
  /*   else if (i == 1 && primalApprox->self->time_method == 3) */
  /*     innertimesolver = 2; */
  /*   else if (i == 2 && primalApprox->self->time_method == 3) */
  /*     innertimesolver = 3; //After two iterations you can now go and use BDF3 */
  /*   //------------------------------------------------------------------------- */
  /*   // Print out the reduced states found at everytime time iteration */
  /*   //------------------------------------------------------------------------- */
  VecDestroy(&p);
  VecDestroy(&residualTemp);
  VecDestroy(&primalReducedVec);
  VecDestroy(&state_0);
  VecDestroy(&primalApproxVec);
  VecDestroy(&st_residual);
  VecDestroy(&primalApproxVec_i);
  MatDestroy(&dRdU);
  MatDestroy(&dRdUTemp);
  MatDestroy(&C);
  MatDestroy(&coefMat);
  MatDestroy(&coefMat1);
  MatDestroy(&coefMat2);
  free(nSTBasisIndex);
  free(state_i);
  free(residualIndex);
  free(spaceIndex);
  free(index);
  free(cc);
}

void yk_gaussNewtonSolve(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Cluster *primalApprox, Is_it *reduced, int tSolve,
			 int *iterationCount){
  int i, j;
   clock_t start, end;
  int flag = 0;
  int innertimesolver = 1;
  int meshSize;
  int gaussCount = 0;
  int systemSize;
  int numNodes;
  int *index;
  int nBasisFuncs = reduced->nBasisFuncs;
  double cpu_time_used;
  double *snapshotsRJ;
  double alfa;
  FILE *residualSnapshotFile;
  FILE *jacobianSnapshotFile;
  PetscScalar normp;
  Vec p;
  Vec residual = NULL, residualTemp = NULL, residualRed = NULL;
  Vec primalReducedVec, primalApproxVec, snapJacobian;
  Mat dRdU, dRdUTemp = NULL, dRdUTempHat = NULL, dRdUbasis = NULL;
  Vec primalRealVec;
  Vec state_0;
  Mesh _meshOfInterest;
  Mat rOBState = NULL;
  Universe _equation = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Mallocing things occur here/Initialize things here
  //---------------------------------------------------------------------------
  systemSize = _equation.numStates*primalApprox->self->space.node.count;
  index = (int *) malloc (systemSize*sizeof(int));
  for (i=0; i<systemSize; i++)
    index[i] = i;
  if (reduced->hrom == 1){
    numNodes = reduced->nSampleNodes;
    rOBState = reduced->rOBStateBar;
    _meshOfInterest = reduced->reducedMesh;
    meshSize = _equation.numStates*reduced->reducedMesh.node.count;
    VecCreate(PETSC_COMM_SELF, &residualRed);
    VecSetSizes(residualRed, numNodes, numNodes);
    VecSetType(residualRed, VECSEQ);
      //primalRealVec
    //VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalRealVec);
    VecCreate(PETSC_COMM_SELF, &primalRealVec);
    VecSetSizes(primalRealVec, primalApprox->stateFull->systemSize,
		primalApprox->stateFull->systemSize);
    VecSetType(primalRealVec, VECSEQ);
    //dRdUbasis
    MatCreate(PETSC_COMM_SELF, &dRdUbasis);
    MatSetSizes(dRdUbasis, meshSize, nBasisFuncs, meshSize, nBasisFuncs);
    MatSetType(dRdUbasis, MATSEQDENSE);
    MatSetUp(dRdUbasis);
    //dRdUTempHat
    MatCreate(PETSC_COMM_SELF, &dRdUTempHat);
    MatSetSizes(dRdUTempHat, numNodes, reduced->nBasisFuncs, numNodes,
		reduced->nBasisFuncs);
    MatSetType(dRdUTempHat, MATSEQDENSE);
    MatSetUp(dRdUTempHat);
    //dRdU
  }else{
    numNodes = primalApprox->self->systemSize;
    residualSnapshotFile = fopen("residualSnapshot.dat", "a");
    jacobianSnapshotFile = fopen("jacobianSnapshot.dat", "a");
    snapshotsRJ = (double *) malloc (numNodes*sizeof(double));
    rOBState = reduced->rOBState;
    meshSize = _equation.numStates*primalApprox->self->space.node.count;
    _meshOfInterest = primalApprox->self->space;
    //snapJacobian
    VecCreate(PETSC_COMM_SELF, &snapJacobian);
    VecSetSizes(snapJacobian, numNodes, numNodes);
    VecSetType(snapJacobian, VECSEQ);
  }


  VecCreate(PETSC_COMM_SELF, &state_0);
  VecSetSizes(state_0, systemSize, systemSize);
  VecSetType(state_0, VECSEQ);


  VecCreate(PETSC_COMM_SELF, &residual);
  VecSetSizes(residual, numNodes, numNodes);
  VecSetType(residual, VECSEQ);
  //residualTemp
  VecCreate(PETSC_COMM_SELF, &residualTemp);
  VecSetSizes(residualTemp, meshSize, meshSize);
  VecSetType(residualTemp, VECSEQ);
  //p direction
  VecCreate(PETSC_COMM_SELF, &p);
  VecSetSizes(p, nBasisFuncs, nBasisFuncs);
  VecSetType(p, VECSEQ);
  //primalReducedVec
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, nBasisFuncs, nBasisFuncs);
  VecSetType(primalReducedVec, VECSEQ);
  //primalApproxVec
  VecCreate(PETSC_COMM_SELF, &primalApproxVec);
  VecSetSizes(primalApproxVec, systemSize, systemSize);
  VecSetType(primalApproxVec, VECSEQ);
  //dRdU
  MatCreate(PETSC_COMM_SELF, &dRdU);
  MatSetSizes(dRdU, numNodes, nBasisFuncs, numNodes, nBasisFuncs);
  MatSetType(dRdU, MATSEQDENSE);
  MatSetUp(dRdU);
  //dRdUTemp
  MatCreate(PETSC_COMM_SELF, &dRdUTemp);
  MatSetSizes(dRdUTemp, meshSize, systemSize, meshSize, systemSize);
  MatSetType(dRdUTemp, MATSEQBAIJ);
  MatSeqBAIJSetPreallocation(dRdUTemp, primalApprox->self->basis.nodes*
			     multiquation->equation.numStates,
			     ykflow->numElemCol, NULL);  //Need to think about
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Rertrieve the initial conditions for the full order model solution
  ks_readSolution(_equation, primalApprox->self, 0);
  array2Vec(_equation, primalApprox->self, state_0);
  //Run for each time step
  for (i = 1; i<primalApprox->self->time.count+1; i++){
  //  for(i=1; i<2; i++){
    flag = 0;
    gaussCount=0;

    ks_readSolution(_equationReduced, primalApprox->reduced, i-1);
    //    ks_setSolutionZero(_equationReduced, primalApprox->reduced);
    array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec);
    primalApprox->self->time.node = i;
    do{
      //-----------------------------------------------------------------------
      // Extract the guesses made from the previous solution
      //-----------------------------------------------------------------------
      //Calculate the approximate solution given initial guess
      MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec);

      //MatMult(rOBState, primalReducedVec, primalApproxVec);
      vec2Array(_equation, primalApprox->self, primalApproxVec);
      //-----------------------------------------------------------------------
      //Calculate the Residual and the Jacobian here
      //-----------------------------------------------------------------------
      ykflow->Residual(ykflow, multiquation, primalApprox, residualTemp,
		       dRdUTemp, reduced, innertimesolver);
      //-----------------------------------------------------------------------
      //Gauss Newton with Approximated Tensors for hyper-reduced order modeling
      //-----------------------------------------------------------------------
      if (reduced->hrom == 1){
        MatMult(reduced->Zmult, residualTemp, residualRed);
	MatMult(reduced->B, residualRed, residual);
	MatMatMultBAIJ(dRdUTemp, rOBState, dRdUbasis);
	MatMatMult(reduced->Zmult, dRdUbasis, MAT_REUSE_MATRIX, PETSC_DEFAULT,
		   &dRdUTempHat);
	MatMatMult(reduced->A, dRdUTempHat, MAT_REUSE_MATRIX, PETSC_DEFAULT,
		   &dRdU);
      }else{
	VecCopy(residualTemp, residual);
	MatMatMultBAIJ(dRdUTemp, rOBState, dRdU);
      }
      //-----------------------------------------------------------------------
      // Residual Snapshots
      //-----------------------------------------------------------------------
      if (reduced->hrom ==0){ //If we are doing regular rom
	VecGetValues(residual, systemSize, index, snapshotsRJ);
	fprintf(residualSnapshotFile, "%0.16f", snapshotsRJ[0]);
	for (j=1; j<systemSize; j++)
	  fprintf(residualSnapshotFile, " %0.16f", snapshotsRJ[j]);
	fprintf(residualSnapshotFile, "\n");
      }
      //-----------------------------------------------------------------------
      // Least Squares Solve for the Direction p
      //-----------------------------------------------------------------------
      VecScale(residual, -1);
      start = clock();
      linearLeastSquares(dRdU, residual, p);
      end = clock();
      //cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC;
      //      printf("%g\n", cpu_time_used);
      (*iterationCount)++;
      //-----------------------------------------------------------------------
      // Jacobian Snapshots
      //-----------------------------------------------------------------------
      if (reduced->hrom == 0){
	MatMult(dRdU, p, snapJacobian);
	VecGetValues(snapJacobian, systemSize, index, snapshotsRJ);
	fprintf(jacobianSnapshotFile, "%0.16f", snapshotsRJ[0]);
	for (j=1; j<systemSize; j++)
	  fprintf(jacobianSnapshotFile, " %0.16f", snapshotsRJ[j]);
	fprintf(jacobianSnapshotFile, "\n");
      }
      //-----------------------------------------------------------------------
      // Perform line search or set \alpha to 1
      //-----------------------------------------------------------------------
      alfa = 1;
      VecNorm(p, 2, &normp);
      printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp);
      gaussCount++;
      if (normp <pow(10, -8))
	flag = 1;
      VecAXPY(primalReducedVec, alfa, p);
      vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
    }while (flag ==0);

    MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec);
    //MatMult(rOBState, primalReducedVec, primalApproxVec);
    vec2Array(_equation, primalApprox->self, primalApproxVec);
    vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec);
    ks_printSolution(_equationReduced, primalApprox->reduced, i);
    ks_printSolution(_equation, primalApprox->self, i);
    if (reduced->hrom == 1){
      array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec);
      MatMult(reduced->rOBState, primalReducedVec, primalRealVec);
      vec2Array(_equation, primalApprox->stateFull, primalRealVec);
      ks_printSolution(_equation, primalApprox->stateFull, i);
    }
    if (i == 1 && primalApprox->self->time_method == 2) //Use BDF1
      innertimesolver = 2;
    else if (i == 1 && primalApprox->self->time_method == 3)
      innertimesolver = 2;
    else if (i == 2 && primalApprox->self->time_method == 3)
      innertimesolver = 3; //After two iterations you can now go and use BDF3
    //-------------------------------------------------------------------------
    // Print out the reduced states found at everytime time iteration
    //-------------------------------------------------------------------------
  }
  //---------------------------------------------------------------------------
  // Destroy all the malloc stuff and initialization stuff OMG
  //---------------------------------------------------------------------------
  VecDestroy(&p);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalApproxVec);
  VecDestroy(&residual);
  VecDestroy(&state_0);
  VecDestroy(&residualTemp);
  MatDestroy(&dRdU);
  MatDestroy(&dRdUTemp);
  if (reduced->hrom == 1){
    VecDestroy(&residualRed);
    VecDestroy(&primalRealVec);
    MatDestroy(&dRdUTempHat);
    MatDestroy(&dRdUbasis);
  }else if (reduced->hrom == 0){
    fclose(residualSnapshotFile);
    fclose(jacobianSnapshotFile);
    VecDestroy(&snapJacobian);
    free(snapshotsRJ);
  }
  free(index);
}
