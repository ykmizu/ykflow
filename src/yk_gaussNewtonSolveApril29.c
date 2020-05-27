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
  int alfa;
  int systemSize = primalApprox->self->systemSize;

  int rowSize = reduced->reducedMesh.elem.count*12;
  int sampleSize = reduced->reducedMesh_Os.elem.count*12;

  int st_systemSize = systemSize*primalApprox->self->time.count;
  if (reduced->hrom==1)
    st_systemSize = rowSize*(reduced->nSampleTime-1);
  int gaussCount = 0;
  int innertimesolver = 1;
  int spaceNodes = primalApprox->self->space.node.count;
  int meshSize = multiquation->equation.numStates*spaceNodes;
  int left_t_node = primalApprox->self->time.t_0/primalApprox->self->time.dt;
  int right_t_node = primalApprox->self->time.t_f/primalApprox->self->time.dt;

  double cpu_time_used;
  int nUnknownsElem = primalApprox->self->basis.nodes*
    multiquation->equation.numStates;
  double *state_i = (double *) malloc (systemSize*sizeof(double));
  double *snapshotsRJ;
  Universe _equation = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  PetscInt local_i;
  PetscInt row, col;
  PetscInt *spaceIndex = (PetscInt *) malloc (systemSize*sizeof(PetscInt ));
  PetscInt *residualIndex = (PetscInt *) malloc (systemSize*sizeof(PetscInt));
  PetscInt *index = (PetscInt *) malloc (systemSize*sizeof(PetscInt ));
  PetscScalar *residual_i;
  Vec p;
  Vec state_0;
  Vec residualTemp;
  Vec primalReducedVec;
  Vec primalApproxVec, primalApproxVec_i, primalApproxVec_i_final;
  Vec st_residual;
  Mat st_rOBState = reduced->ST_rOBState;
  Mat dRdU=NULL, dRdUTemp = NULL;
  Mat C;
  Mat coefMat1;
  Mat coefMat2;
  PetscInt ggrow, ggcol;
  Mesh  _meshOfInterest = primalApprox->self->space;
  PetscInt rowBasis, basis;
  PetscInt *rowIndex;
  PetscInt *nSTBasisIndex;
  PetscScalar normp;
  clock_t start, end;
  FILE *residualSnapshotFile;
  PetscScalar *cc;
  PetscInt we = 0;
  int flag = 0;
  PetscInt blockcount_i;
  Mat coefMat;
  PetscReal val;
  int count = 0;
  Mat r_dRdU;
  Vec r_st_residual;
  Vec state_0_final;
  PetscInt rowRes, colRes;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->hrom == 0 )
    MatGetSize(reduced->ST_rOBState, &row, &basis);
  else
    MatGetSize(reduced->ST_rOBStateBar,  &row, &basis);

  if (reduced->hrom == 1)
    MatGetSize(reduced->rOBResidual, &rowRes, &colRes);

  snapshotsRJ = (double *) malloc (row*sizeof(double));
  rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt));
  nSTBasisIndex = (PetscInt *) malloc (basis*sizeof(PetscInt));
  cc = (PetscScalar *) malloc (systemSize*basis*sizeof(PetscScalar));

  if (reduced->hrom==0)
    residual_i = (PetscScalar *) malloc (row*sizeof(double));
  else if (reduced->hrom == 1)
    residual_i = (PetscScalar *) malloc (reduced->reducedMesh.node.count*4
					 *sizeof(double));
  for (i=0; i<basis; i++)
    nSTBasisIndex[i] = i;

  for (i=0; i<systemSize; i++)
    index[i] = i;

  for (i=0; i<row; i++)
    rowIndex[i] = i;
  PetscInt *samIndex = (PetscInt *) malloc (sampleSize*sizeof(PetscInt));

  for (i=0; i<sampleSize; i++)
    samIndex[i] = i;

  PetscInt *rowSamIndex = (PetscInt *) malloc (rowSize*sizeof(PetscInt));
  for (i=0; i<rowSize; i++)
    rowSamIndex[i] = i;
  VecCreate(PETSC_COMM_SELF, &p);
  VecSetSizes(p, basis, basis);
  VecSetType(p, VECSEQ);
  //residualTemp
  if (reduced->hrom==0){
    VecCreate(PETSC_COMM_SELF, &residualTemp);
    VecSetSizes(residualTemp, systemSize, systemSize);
    VecSetType(residualTemp, VECSEQ);
  }else if (reduced->hrom==1){
    VecCreate(PETSC_COMM_SELF, &residualTemp);
    VecSetSizes(residualTemp, rowSize,
  		rowSize);
    VecSetType(residualTemp, VECSEQ);
  }
  //primalReducedVec
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, basis, basis);
  VecSetType(primalReducedVec, VECSEQ);


  VecCreate(PETSC_COMM_SELF, &state_0_final);
  VecSetSizes(state_0_final, systemSize, systemSize);
  VecSetType(state_0_final, VECSEQ);

  if (reduced->hrom==0){
    VecCreate(PETSC_COMM_SELF, &state_0);
    VecSetSizes(state_0, systemSize, systemSize);
    VecSetType(state_0, VECSEQ);
  }else{
    VecCreate(PETSC_COMM_SELF, &state_0);
    VecSetSizes(state_0, sampleSize, sampleSize);
    VecSetType(state_0, VECSEQ);
  }


  VecCreate(PETSC_COMM_SELF, &primalApproxVec);
  VecSetSizes(primalApproxVec, row, row);
  VecSetType(primalApproxVec, VECSEQ);

  if (reduced->hrom==0){
    VecCreate(PETSC_COMM_SELF, &st_residual);
    VecSetSizes(st_residual, PETSC_DECIDE, row);
    VecSetBlockSize(st_residual, systemSize);
    VecSetFromOptions(st_residual);
  }else if (reduced->hrom==1){
    VecCreate(PETSC_COMM_SELF, &st_residual);
    VecSetSizes(st_residual, PETSC_DECIDE, rowSize*(reduced->nSampleTime-1));
    VecSetBlockSize(st_residual, rowSize);
    VecSetFromOptions(st_residual);

  }

  VecCreate(PETSC_COMM_SELF, &r_st_residual);
  VecSetSizes(r_st_residual, colRes, colRes);
  VecSetType(r_st_residual, VECSEQ);

  if (reduced->hrom==0){
    VecCreate(PETSC_COMM_SELF, &primalApproxVec_i);
    VecSetSizes(primalApproxVec_i, systemSize, systemSize);
    VecSetType(primalApproxVec_i, VECSEQ);
  }else if (reduced->hrom==1){
    VecCreate(PETSC_COMM_SELF, &primalApproxVec_i);
    VecSetSizes(primalApproxVec_i, sampleSize, sampleSize);
    VecSetType(primalApproxVec_i, VECSEQ);
  }


  VecCreate(PETSC_COMM_SELF, &primalApproxVec_i_final);
  VecSetSizes(primalApproxVec_i_final, systemSize, systemSize);
  VecSetType(primalApproxVec_i_final, VECSEQ);

 //dRdU
  MatCreate(PETSC_COMM_SELF, &dRdU);
  MatSetSizes(dRdU, PETSC_DECIDE, basis, st_systemSize, basis);
  if (reduced->hrom == 0)
    MatSetBlockSizes(dRdU, systemSize, basis);
  else if (reduced->hrom == 1)
    MatSetBlockSizes(dRdU, rowSize, basis);

  MatSetType(dRdU, MATSEQDENSE);
  MatSetFromOptions(dRdU);
  MatSetUp(dRdU);

  if (reduced->hrom==0){
    MatCreate(PETSC_COMM_SELF, &dRdUTemp);
    MatSetSizes(dRdUTemp, meshSize, systemSize, meshSize, systemSize);
    MatSetType(dRdUTemp, MATSEQBAIJ);
    MatSeqBAIJSetPreallocation(dRdUTemp, nUnknownsElem,ykflow->numElemCol, NULL);
  }else if (reduced->hrom==1){
    MatCreate(PETSC_COMM_SELF, &dRdUTemp);
    MatSetSizes(dRdUTemp, rowSize,
  		reduced->reducedMesh_Os.node.count*4,
  		rowSize,
  		reduced->reducedMesh_Os.node.count*4);
    MatSetType(dRdUTemp, MATSEQBAIJ);
    MatSeqBAIJSetPreallocation(dRdUTemp, nUnknownsElem,ykflow->numElemCol, NULL);
  }
  MatCreate(PETSC_COMM_SELF, &C);
  if (reduced->hrom==0)
    MatSetSizes(C, systemSize, basis, systemSize, basis);
  else if (reduced->hrom == 1)
    MatSetSizes(C, rowSize, basis, rowSize, basis);
  MatSetType(C, MATSEQDENSE);
  MatSetUp(C);

  //---------------------------------------------------------------------------
  // Need to build up the Residual and the Jaociban vector and matrix
  //---------------------------------------------------------------------------
  /* numNodes = primalApprox->self->systemSize; */
  residualSnapshotFile = fopen("residualSnapshot.dat", "w+");
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Rertrieve the initial conditions for the full order model solution
  ks_readSolution(_equationReduced, primalApprox->reduced, left_t_node);

  array2Vec(_equationReduced, primalApprox->reduced,
	    primalApprox->reduced->space, primalReducedVec);
  //initial conditions (I don't think I need this
  ks_readSolution(_equation, primalApprox->self, left_t_node);
  if (reduced->hrom==0)
    array2Vec(_equation, primalApprox->self, primalApprox->self->space,
	      state_0);
  else if (reduced->hrom==1)
    array2Vec(_equation, primalApprox->self, reduced->reducedMesh_Os,
	      state_0);

  if (reduced->hrom ==0){
    MatCreate(PETSC_COMM_SELF, &coefMat);
    MatSetSizes(coefMat, systemSize, systemSize, systemSize,
		systemSize);
    MatSetType(coefMat, MATSEQAIJ);
    MatSeqAIJSetPreallocation(coefMat, primalApprox->self->basis.nodes*4*12, NULL);
  }else if (reduced->hrom == 1){
    MatCreate(PETSC_COMM_SELF, &coefMat);
    MatSetSizes(coefMat, rowSize, sampleSize, rowSize, sampleSize);
    MatSetType(coefMat, MATSEQAIJ);
    MatSeqAIJSetPreallocation(coefMat, primalApprox->self->basis.nodes*4*12, NULL);
  }

  ykflow->findMassCoefBDF(ykflow, primalApprox->self, &coefMat, reduced);

  MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat1);
  MatScale(coefMat1, -2/0.1);
  //MatScale(coefMat1, -1/0.1);
  MatDuplicate(coefMat, MAT_COPY_VALUES, &coefMat2);
  MatScale(coefMat2, 0.5/0.1);
  PetscInt vvsize;


  do{
    //-------------------------------------------------------------------------
    // Loop through each time to obtain the residual and Jacobian at tstep i
    //-------------------------------------------------------------------------
    //Calculate just the Pi*reducedstates part. Add in initial conditions in
    //for loop
    MatZeroEntries(dRdU);
    //for (i = 1; i<primalApprox->self->time.count+1; i++){

    /* for (i=1; i<primalApprox->self->time.tNode_i.count+1; i++){ */
    for (t_i = 0; t_i<primalApprox->self->time.tNode_i.count; t_i++){
      i = left_t_node+primalApprox->self->time.tNode_i.array[t_i]+1;
      local_i = primalApprox->self->time.tNode_i.array[t_i]+1;
      //local_i = +1;
      primalApprox->self->time.node = i;


      /* local_i = primalApprox->self->time.tNode_i.array[i-1]+1; */
      /* primalApprox->self->time.node = left_t_node+i; */
    /* for (i=left_t_node+1; i<right_t_node+1; i++){ */

    /*   local_i = i-left_t_node; */
    /*   primalApprox->self->time.node = i; */
      //temp temp
      /* printf("i %d\n", i); */
      /* printf("local i %d\n", local_i); */
      /* printf("global i  %d\n", primalApprox->self->time.node); */



      //Write to file
      if (reduced->hrom==0)
	MatMultAdd(reduced->ST_rOBState_i[local_i-1], primalReducedVec,
		   state_0, primalApproxVec_i);
      else if (reduced->hrom==1)
	MatMultAdd(reduced->ST_rOBStateBar_i[local_i-1], primalReducedVec,
		   state_0, primalApproxVec_i);

      if (reduced->hrom==0)
	vec2Array(_equation, primalApprox->self, primalApprox->self->space,
		  primalApproxVec_i);
      else if (reduced->hrom==1)
	vec2Array(_equation, primalApprox->self, reduced->reducedMesh_Os,
                  primalApproxVec_i);
      ks_printSolution(_equation, primalApprox->self, i);
      /* ks_readSolution(_equation, primalApprox->self, 1); */
      /* printf("working working working\n"); */
      //-----------------------------------------------------------------------
      //Calculate the Residual and the Jacobian here
      //-----------------------------------------------------------------------
      ykflow->Residual(ykflow, multiquation, primalApprox, residualTemp,
      		       dRdUTemp, reduced, innertimesolver);


      /* MatView(dRdUTemp, PETSC_VIEWER_STDOUT_SELF); */
      /* printf("jacobian lets see what will happen here now \n"); */
      /* getchar(); */

      VecNorm(residualTemp, 2, &val);
      //-----------------------------------------------------------------------
      //Add the residual for time i to a large ass space time vector
      //-----------------------------------------------------------------------
      VecGetArray(residualTemp, &residual_i);

      //blockcount_i = local_i-1;
      blockcount_i = t_i;
      VecSetValuesBlocked(st_residual, 1, &blockcount_i, residual_i, INSERT_VALUES);
      //-----------------------------------------------------------------------
      // Jacobain manipulation
      //-----------------------------------------------------------------------
      start = clock();
      //Oh fuck
      if (reduced->hrom == 0)
	MatMatMultBAIJ(dRdUTemp, reduced->ST_rOBState_i[local_i-1], C);
      else if (reduced->hrom == 1)
	MatMatMultBAIJ(dRdUTemp, reduced->ST_rOBStateBar_i[local_i-1], C);
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

      if (reduced->hrom==0)
	MatGetValues(C, systemSize,  primalApprox->self->index, basis,
		     nSTBasisIndex, cc);
      else if (reduced->hrom==1)
	MatGetValues(C, rowSize, rowSamIndex, basis, nSTBasisIndex, cc);

      MatSetValuesBlocked(dRdU, 1, &blockcount_i,1, &we, cc,ADD_VALUES);

      if (local_i>1){
	if (reduced->hrom==0){
	  MatMatMultBAIJ(coefMat1, reduced->ST_rOBState_i[local_i-2], C);
	  MatGetValues(C, systemSize, primalApprox->self->index, basis,
      		     nSTBasisIndex, cc);

	}else if (reduced->hrom==1){
	  MatMatMultBAIJ(coefMat1, reduced->ST_rOBStateBar_i[local_i-2], C);
	  MatGetValues(C, rowSize, rowSamIndex, basis, nSTBasisIndex, cc);
	}
	MatSetValuesBlocked(dRdU, 1, &blockcount_i,1, &we, cc,ADD_VALUES);

      }

      if (local_i>2){ //BDF 2 section
	if (reduced->hrom == 0){
      	MatMatMultBAIJ(coefMat2, reduced->ST_rOBState_i[local_i-3], C);
      	MatGetValues(C, systemSize, primalApprox->self->index,
      		     basis, nSTBasisIndex, cc);

	}else if (reduced->hrom == 1){
          MatMatMultBAIJ(coefMat2, reduced->ST_rOBStateBar_i[local_i-3], C);
       	  MatGetValues(C, rowSize, rowSamIndex, basis, nSTBasisIndex, cc);

	}
      	MatSetValuesBlocked(dRdU, 1, &blockcount_i, 1, &we, cc, ADD_VALUES);
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
      count++;
      VecGetValues(st_residual, row, rowIndex, snapshotsRJ);
      fprintf(residualSnapshotFile, "%0.16f", snapshotsRJ[0]);
      for (j=1; j<row; j++)
	fprintf(residualSnapshotFile, " %0.16f", snapshotsRJ[j]);
      fprintf(residualSnapshotFile, "\n");
    }

    //-----------------------------------------------------------------------
    // Least Squares Solve for the Direction p
    //-----------------------------------------------------------------------
    VecScale(st_residual, -1);
  /*   start = clock(); */
    PetscInt m, n;
    if (reduced->hrom==1){
      MatGetSize(reduced->A, &ggrow, &ggcol);
      MatMatMult(reduced->A, dRdU, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &r_dRdU);
      VecGetSize(st_residual, &m);
      MatGetSize(reduced->A, &m, &n);
      MatMult(reduced->A, st_residual, r_st_residual);
      linearLeastSquares(r_dRdU, r_st_residual, p);

    }else if (reduced->hrom == 0){
      linearLeastSquares(dRdU, st_residual, p);
    }
    /*   end = clock(); */
  /*   //cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC; */
  /*   //      printf("%g\n", cpu_time_used); */
  /*   (*iterationCount)++; */
    //-----------------------------------------------------------------------
    // Perform line search or set \alpha to 1
    //-----------------------------------------------------------------------
    alfa = 1;
    gaussCount++;
    VecNorm(p, 2, &normp);
    printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp);
    /* gaussCount++;		/\*  *\/ */
    if (normp <pow(10, -10))
      flag = 1;
    VecAXPY(primalReducedVec, alfa, p);
    vec2Array(_equationReduced, primalApprox->reduced,
	      primalApprox->reduced->space, primalReducedVec);
  }while (flag ==0);

  //Print the reduced solution based on the time window number
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);

  ks_readSolution(_equation, primalApprox->self, left_t_node);
  array2Vec(_equation, primalApprox->self, primalApprox->self->space,
              state_0_final);
  //Should print out the approximated entire final solution too
  for (i=left_t_node+1; i<right_t_node+1; i++){
    local_i = i-left_t_node;
    MatMultAdd(reduced->ST_rOBState_i[local_i-1], primalReducedVec, state_0_final,
	       primalApproxVec_i_final);
    if (reduced->hrom==1){
      vec2Array(_equation, primalApprox->stateFull, primalApprox->stateFull->space, primalApproxVec_i_final);
      ks_printSolution(_equation, primalApprox->stateFull, i);
    }else if (reduced->hrom==0){
      vec2Array(_equation, primalApprox->self, primalApprox->self->space,
		primalApproxVec_i_final);
      ks_printSolution(_equation, primalApprox->self, i);
    }
  }
  /*   MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec); */
  /*   //MatMult(rOBState, primalReducedVec, primalApproxVec); */
  /*   vec2Array(_equation, primalApprox->self, primalApproxVec); */
  /*   vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /*   ks_printSolution(_equationReduced, primalApprox->reduced, i); */
  /*   if (reduced->hrom == 1){ */
  /*     array2Vec(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /*     MatMult(reduced->rOBState, primalReducedVec, primalRealVec); */
  /*     vec2Array(_equation, primalApprox->stateFull, primalRealVec); */
  /*     ks_printSolution(_equation, primalApprox->stateFull, i); */
  /*   } */
    //-------------------------------------------------------------------------
    // Print out the reduced states found at everytime time iteration
    //-------------------------------------------------------------------------
  fclose(residualSnapshotFile);

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
  free(snapshotsRJ);
  free(nSTBasisIndex);
  free(state_i);
  free(rowIndex);
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
  array2Vec(_equation, primalApprox->self, primalApprox->self->space, state_0);
  //Run for each time step
  for (i = 1; i<primalApprox->self->time.count+1; i++){
  //  for(i=1; i<2; i++){
    flag = 0;
    gaussCount=0;

    ks_readSolution(_equationReduced, primalApprox->reduced, i-1);
    //    ks_setSolutionZero(_equationReduced, primalApprox->reduced);
    array2Vec(_equationReduced, primalApprox->reduced, primalApprox->reduced->space, primalReducedVec);
    primalApprox->self->time.node = i;
    do{
      //-----------------------------------------------------------------------
      // Extract the guesses made from the previous solution
      //-----------------------------------------------------------------------
      //Calculate the approximate solution given initial guess
      MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec);

      //MatMult(rOBState, primalReducedVec, primalApproxVec);
      vec2Array(_equation, primalApprox->self, primalApprox->self->space,
		primalApproxVec);
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
      vec2Array(_equationReduced, primalApprox->reduced,
		primalApprox->reduced->space, primalReducedVec);
    }while (flag ==0);

    MatMultAdd(rOBState, primalReducedVec, state_0, primalApproxVec);
    //MatMult(rOBState, primalReducedVec, primalApproxVec);
    vec2Array(_equation, primalApprox->self, primalApprox->self->space, primalApproxVec);
    vec2Array(_equationReduced, primalApprox->reduced, primalApprox->self->space, primalReducedVec);
    ks_printSolution(_equationReduced, primalApprox->reduced, i);
    ks_printSolution(_equation, primalApprox->self, i);
    if (reduced->hrom == 1){
      array2Vec(_equationReduced, primalApprox->reduced, primalApprox->reduced->space, primalReducedVec);
      MatMult(reduced->rOBState, primalReducedVec, primalRealVec);
      vec2Array(_equation, primalApprox->stateFull, primalApprox->stateFull->space,		primalRealVec);
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
