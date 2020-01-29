#include "yk_performPOD.h"

Mat * yk_createSnapshotState(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, Mat *snapshot, int snapshotStep,
			    Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k, p;                       //initialization for iteration
  int snapshot_i;                       //state to be included in snapshot
  /* int snapshotCount = primal->self->time.count; //256 */
  /* int snapshotCount = (primal->self->time.t_f-primal->self->time.t_0)/ */
  /*   primal->self->time.dt; */
  /* int snapshotCount = ceil(((primal->self->time.t_f-primal->self->time.t_0)*10)/ */
  /*                          (primal->self->time.dt*10)); */
  int snapshotCount = primal->self->time.count;
  int totalParamSnapshot = snapshotCount*reduced->numParamSet;
  PetscInt index_i;
  PetscScalar *state_i =                //state array to add to snapshot matrix
    (PetscScalar *) malloc (primal->self->systemSize*sizeof(PetscScalar));
  PetscScalar *state_0 =                //state array to add to snapshot matrix
    (PetscScalar *) malloc (primal->self->systemSize*sizeof(PetscScalar));
  char cwd[1024];
  char fomBuffer[50];
  Mat *snapshot_mu;
  PetscInt row, col;
  PetscInt *temp = (PetscInt *) malloc (primal->self->systemSize*sizeof(PetscInt));
  for (i=0; i<primal->self->systemSize; i++)
    temp[i] = i;
  PetscMalloc1(reduced->numParamSet, &snapshot_mu);
  int timeNode0 = primal->self->time.t_0/primal->self->time.dt;
 //  Mat *snapshot_mu = (Mat *) malloc (reduced->numParamSet*sizeof(Mat));
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatCreate(PETSC_COMM_SELF, snapshot);//create snapshot matrix
  MatSetSizes(*snapshot, primal->self->systemSize, totalParamSnapshot,
	      primal->self->systemSize, totalParamSnapshot);
  MatSetType(*snapshot, MATSEQDENSE);           //set snapshot matrix type
  MatSetUp(*snapshot);

  for (i=0; i<reduced->numParamSet; i++){
    MatCreate(PETSC_COMM_SELF, &snapshot_mu[i]);//create snapshot matrix
    MatSetSizes(snapshot_mu[i], primal->self->systemSize, snapshotCount,
		primal->self->systemSize, snapshotCount);
    MatSetType(snapshot_mu[i], MATSEQDENSE);           //set snapshot matrix t
    MatSetUp(snapshot_mu[i]);
  }
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (p = 0; p< reduced->numParamSet; p++){
    //Retrieve the data set name and enter that directory
    ykflow->FomName(ykflow, multiquation->equation, reduced, p, fomBuffer);
    chdir(fomBuffer);
    getcwd(cwd, sizeof(cwd));
    //Read the state solutions for the initial time in that directory
    ks_readSolution(multiquation->equation, primal->self, timeNode0);
    array2Data(multiquation->equation, primal->self, state_0);
    //Add data to the Petsc Mat Type
    for (i=1; i<snapshotCount+1; i++){      //begin creating snapshot matrix
      ks_readSolution(multiquation->equation, primal->self, timeNode0+i);
      array2Data(multiquation->equation, primal->self, state_i);
      //Take away the initial conditions from the states
      for (j=0; j<primal->self->systemSize; j++)
      	state_i[j]-=state_0[j];
      MatGetSize(*snapshot, &row, &col);

      index_i = p*(snapshotCount)+i-1;
      //Add to the giant snapshot matrix
      MatSetValues(*snapshot, primal->self->systemSize, temp, 1,
    		   &index_i, state_i, INSERT_VALUES);
      index_i = i-1;
      MatSetValues(snapshot_mu[p], primal->self->systemSize,
		   temp,
		   1, &index_i, state_i, INSERT_VALUES);
    }
    MatAssemblyBegin(snapshot_mu[p], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(snapshot_mu[p], MAT_FINAL_ASSEMBLY);
    chdir("../");
  }

  MatAssemblyBegin(*snapshot, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*snapshot, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  // Termination
  //---------------------------------------------------------------------------
  free(state_0);
  free(state_i);
  free(temp);
  //Optional return of snapshot_mu for space time stuff
  return snapshot_mu;
}

void yk_createTemporalSnapshotState(yk_PrimalSolver *ykflow,
				    Multiverse *multiquation, Cluster *primal,
				    Mat *snapshot_phi, Mat *s_mu, int snapshotStep,
				    Is_it *reduced){
  int i, j, p; //initialization for iteration
  Vec state_0;
  Vec state_i;
  Vec phi;
  PetscScalar val;
  PetscInt row=primal->self->systemSize;
  //Number of colum.ns in that snapshot matrix basically
  int col = primal->self->time.count;
  /* int col = ceil(((primal->self->time.t_f-primal->self->time.t_0)*10)/ */
  /*                          (primal->self->time.dt*10)); */
/* PetscInt col = (primal->self->time.t_f-primal->self->time.t_0)/ */
  /*   primal->self->time.dt; */
  /* PetscInt col=primal->self->time.count; */
  char cwd[1024];
  char fomBuffer[50];
  PetscInt index_i;
  PetscInt snapshot_i;
  Vec psi_p;
  PetscInt *rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt));
  PetscInt *colIndex = (PetscInt *) malloc (col*sizeof(PetscInt));
  for (i=0; i< row; i++)
    rowIndex[i] = i;
  for (i=0; i<col; i++)
    colIndex[i] = i;
  PetscInt rbasisrow, rbasiscol;
  MatGetSize(reduced->rOBState, &row, &rbasiscol);
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  for (i=0; i<rbasiscol; i++){
    MatCreate(PETSC_COMM_SELF, &snapshot_phi[i]);//create snapshot matrix
    MatSetSizes(snapshot_phi[i], col, reduced->numParamSet, col,
		reduced->numParamSet);
    MatSetType(snapshot_phi[i], MATSEQDENSE);        //set snapshot matrix type
    MatSetUp(snapshot_phi[i]);
  }
  VecCreate(PETSC_COMM_SELF, &phi);
  VecSetSizes(phi, primal->self->systemSize, primal->self->systemSize);
  VecSetType(phi, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &psi_p);
  VecSetSizes(psi_p, col, col);
  VecSetType(psi_p, VECSEQ);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (j = 0; j<rbasiscol; j++){
    ks_MatCol2Vec(reduced->rOBState, row, rowIndex, phi, j, INSERT_VALUES);
    for (p=0; p<reduced->numParamSet; p++){
      MatMultTranspose(s_mu[p], phi, psi_p);
      ks_Vec2MatCol(snapshot_phi[j], col, colIndex, p, psi_p, INSERT_VALUES);
    }
    MatAssemblyBegin(snapshot_phi[j], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(snapshot_phi[j], MAT_FINAL_ASSEMBLY);
  }
  VecDestroy(&psi_p);
  VecDestroy(&phi);
  VecDestroy(&psi_p);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  /* ks_readSolution(multiquation->equation, primal->self, 0); */
  /* array2Vec(multiquation->equation, primal->self, state_0); */
  /* for (i=0; i<row; i++) */
  /*   rowIndex[i] = i; */
  /* for (j = 0; j<reduced->nBasisFuncs; j++){ */

  /*   ks_MatCol2Vec(reduced->rOBState, row, rowIndex, phi, j, INSERT_VALUES); */
  /*   for (p=0; p<reduced->numParamSet; p++){ */
  /*     ykflow->FomName(ykflow, multiquation->equation, reduced, p, fomBuffer); */
  /*     chdir(fomBuffer); */
  /*     getcwd(cwd, sizeof(cwd)); */
  /*     //Read the state solutions for the initial time in that directory */
  /*     //Add data to the Petsc Mat Type */
  /*     for (i=0; i<snapshotCount; i++){      //begin creating snapshot matrix */
  /* 	index_i = p*snapshotCount+i; */
  /* 	if (i==snapshotCount-1) */
  /* 	  snapshot_i = primal->self->time.count; //add the last state */
  /* 	else if (i==0) */
  /* 	  snapshot_i = 1;                //do not include the initial condition */
  /* 	else */
  /* 	  snapshot_i = snapshotStep*i; */
  /* 	//extract the state solution and add to the snapshot matrix */
  /* 	ks_readSolution(multiquation->equation, primal->self, snapshot_i); */
  /* 	array2Vec(multiquation->equation, primal->self, state_i); */
  /* 	VecAXPY(state_i, -1, state_0); */
  /* 	VecDot(state_i, phi, &val); */
  /* 	MatSetValues(*snapshot, 1, &i, 1, &p, &val, INSERT_VALUES); */
  /* 	/\* MatSetValue(snapshot, i, p, val, INSERT_VALUES); *\/ */
  /*     } */
  /*     chdir("../"); */
  /*   } */
  /* } */
  /* MatView(snapshot, PETSC_VIEWER_STDOUT_SELF); */
  /* getchar(); */
  free(rowIndex);
  free(colIndex);
}

void yk_properOrthogonalDecomposeGroup(Mat *snapshot, int numSnapshots,
				       int systemSize, PetscInt *index,
				       PetscInt *numSingularValues,
				       PetscScalar engyValues, Mat *A){
  int i;
  /* for (i=0; i<numSnapshots; i++){ */
  /*   MatCreate(PETSC_COMM_SELF, &A[i]);       //Set up the Matrix for POD */
  /*   MatSetSizes(A[i], systemSize, numSingularValues, systemSize, */
  /* 		numSingularValues); */
  /*   MatSetType(A[i], MATSEQDENSE); */
  /*   MatSetUp(A[i]); */
  /* } */
  for (i=0; i<numSnapshots; i++){
    yk_properOrthogonalDecompose(&snapshot[i], systemSize, index,
				 numSingularValues, engyValues, &A[i]);
  }
}


void yk_properOrthogonalDecompose(Mat *snapshot, int systemSize,
				  PetscInt *index,
				  PetscInt *numSingularValues,
				  PetscScalar engyValues, Mat *A){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                                //initialization for iteration
  PetscReal sigma;                      //eigenvalues
  Vec singular;                         //singular vectors rotation
  SVD svd;       //create svd object
  PetscInt nsv;
  PetscInt nconv;
  int flag = 0;
  double denomE = 0;
  double Etotal = 0;
  double EtotalPercent = 0;
  PetscInt count = 0;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  VecCreate(PETSC_COMM_SELF, &singular);//Set up the singular vector, U
  VecSetSizes(singular, systemSize, systemSize);
  VecSetType(singular, VECSEQ);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Create the Singular Value solver and set various options and solve
  //---------------------------------------------------------------------------
  SVDCreate(PETSC_COMM_SELF, &svd);
  SVDSetOperator(svd, *snapshot);
  SVDSetFromOptions(svd);
  SVDSolve(svd);
  SVDGetDimensions(svd,&nsv,NULL,NULL);
  SVDGetConverged(svd,&nconv);

  //---------------------------------------------------------------------------
  // If nBasis is undetermined, than calculate witht he energy desire thing
  //---------------------------------------------------------------------------
  if (engyValues != 0 && *numSingularValues == 0){
    for (i=0; i<nconv; i++){
      SVDGetSingularTriplet(svd, i, &sigma, NULL, NULL);
      denomE +=sigma*sigma;
    }
    while (EtotalPercent < engyValues){
      SVDGetSingularTriplet(svd, *numSingularValues, &sigma, singular, NULL);
      Etotal += (sigma*sigma);
      EtotalPercent = Etotal/denomE;
      count++;
    }
  }
  //---------------------------------------------------------------------------
  // Create Basis Matrix here
  //---------------------------------------------------------------------------
  MatCreate(PETSC_COMM_SELF, A);       //Set up the Matrix for POD
  MatSetSizes(*A, systemSize, count, systemSize, count);
  MatSetType(*A, MATSEQDENSE);
  MatSetUp(*A);
  for (i=0; i<count; i++){
    SVDGetSingularTriplet(svd, i, &sigma, singular, NULL); //extract rotation
    ks_Vec2MatCol(*A, systemSize, index, i, singular, INSERT_VALUES);
  }

  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  // Clean up
  //---------------------------------------------------------------------------
  SVDDestroy(&svd);
  VecDestroy(&singular);
}


void yk_createSpaceTimeBasis(Mat *spacebasis, Mat *timebasis, Mat *spaceTime,
			     Mat *spaceTime_i){

  /* yk_kron(&st_basis_temp, &reduced->rOBState, reduced->nBasisFuncs);   */
  Vec phi;
  Mat testest;
  PetscInt row, rowTime;
  PetscInt col, colTime;
  int i, j; //initialization for iteration
  MatGetSize(*spacebasis, &row, &col);

  PetscInt *rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt));
  for (i=0; i< row; i++)
    rowIndex[i] = i ;
  MatGetSize(*timebasis, &rowTime, &colTime);
  PetscInt *rowstIndex = (PetscInt *) malloc (row*rowTime*sizeof(PetscInt));
  double **kronpiece;
  for (i=0; i<row*rowTime; i++)
    rowstIndex[i] = i;
  PetscInt *colstIndex = (PetscInt *) malloc (colTime*sizeof(PetscInt));
  PetscInt *colIndex = (PetscInt *) malloc (col*colTime*sizeof(PetscInt));
  PetscInt *index_i = (PetscInt *) malloc (row*sizeof(PetscInt));

  PetscScalar *st_basisA = (PetscScalar *) malloc
    (row*col*colTime*sizeof(PetscScalar));
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  VecCreate(PETSC_COMM_SELF, &phi);
  VecSetSizes(phi, row, row);
  VecSetType(phi, VECSEQ);
  MatCreate(PETSC_COMM_SELF, spaceTime);       //Set up the Matrix for POD
  MatSetSizes(*spaceTime, row*rowTime, col*colTime, row*rowTime, col*colTime);
  MatSetType(*spaceTime, MATSEQDENSE);
  MatSetUp(*spaceTime);
  for (i=0; i<rowTime; i++){
    MatCreate(PETSC_COMM_SELF, &spaceTime_i[i]);
    MatSetSizes(spaceTime_i[i], row, col*colTime, row, col*colTime);
    MatSetType(spaceTime_i[i], MATSEQDENSE);
    MatSetUp(spaceTime_i[i]);
  }
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i = 0; i<col; i++){
    for (j=0; j<colTime; j++)
      colstIndex[j] = i*colTime+j;
    //Get the first basis function vector
    ks_MatCol2Vec(*spacebasis, row, rowIndex, phi, i, INSERT_VALUES);
    //Take the kronbecker product of the basis vector and the time basis
    yk_kron(timebasis[i], phi, &kronpiece);
    MatSetValues(*spaceTime, row*rowTime, rowstIndex, colTime, colstIndex,
		 *kronpiece, INSERT_VALUES);
  }
  MatAssemblyBegin(*spaceTime, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*spaceTime, MAT_FINAL_ASSEMBLY);


  /* for (i=0; i<col*colTime; i++) */
  /*   colIndex[i] = i; */

  /* for (i=0; i<rowTime; i++){ */
  /*   for (j=0; j<row; j++) */
  /*     index_i[j] = i*row+j; */
  /*   MatGetValues(*spaceTime, row, index_i, col*colTime, colIndex, st_basisA); */
  /*   /\* printf("%g\n", st_basisA[0][0]); *\/ */
  /*   MatSetValues(spaceTime_i[i], row, rowIndex, col*colTime, colIndex, */
  /* 		 st_basisA, INSERT_VALUES); */
  /*   MatAssemblyBegin(spaceTime_i[i], MAT_FINAL_ASSEMBLY); */
  /*   MatAssemblyEnd(spaceTime_i[i], MAT_FINAL_ASSEMBLY); */
  /* } */

  VecDestroy(&phi);
  free(kronpiece[0]);
  free(kronpiece);

  free(rowIndex);
  free(rowstIndex);
  free(colIndex);
  free(colstIndex);
  free(index_i);
  free(st_basisA);
}

/* void yk_createUltimateSpaceTimeBasis(yk_PrimalSolver *ykflow, */
/* 				     Multiverse *multiquation, Cluster *primal, */
/* 				     Cluster *primalApproxm, Is_it *reduced){ */
/*   Mat *st_basis_temp; */
/*   Mat *st_snapshot; */
/*   Mat *s_mu; */
/*   Mat snapshot;                         // Snapshot Matrix before POD */
/*   PetscInt rows, nBasis; */
/*   PetscInt snapshotCount; */
/*   /\* int globalSnapshotCount = primal->self->time.count; *\/ */
/*   /\* int snapshotCount = globalSnapshotCount/  *\/ */
/*   //------------------------------------//------------------------------------- */
/*   // Variables here                     // Comments section */
/*   //------------------------------------//------------------------------------- */

/*   //--------------------------------------------------------------------------- */
/*   // Initialization */
/*   //--------------------------------------------------------------------------- */
/*   reduced->ST_rOBState_i = (Mat *) malloc (snapshotCount*sizeof(Mat)); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   // For each time window, we need to build the basis based on # of time window */
/*   for (i=0; i<nSub; i++){ //Need to iterate the number of sub time windows */
/*     //------------------------------------------------------------------------- */
/*     // FOr eac subspace time window calcualte the original POD for space */
/*     //------------------------------------------------------------------------- */
/*     s_mu = yk_createSnapshotState(ykflow, multiquation, primal, &snapshot, */
/* 				  reduced->dss, reduced); */
/*     yk_properOrthogonalDecompose(&snapshot, primal->self->systemSize, */
/* 				 primal->self->index, &reduced->nBasisFuncs, */
/* 				 reduced->eBasisSpace, &reduced->rOBState); */
/*     //------------------------------------------------------------------------- */
/*     // Now need to calcualte the temmporal POD and et.c TIME */
/*     //------------------------------------------------------------------------- */
/*     //Need to find */
/*     MatGetSize(reduced->rOBState, &rows, &nBasis); */
/*     PetscMalloc1(nBasis, &st_snapshot); */
/*     PetscMalloc1(nBasis, &st_basis_temp); */

/*     yk_createTemporalSnapshotState(ykflow, multiquation, primal, st_snapshot, */
/* 				   s_mu, */
/* 				   reduced->dss, reduced); */
/*     yk_properOrthogonalDecomposeGroup(st_snapshot, nBasis, snapshotCount, */
/* 				      primal->self->index, */
/* 				      &reduced->nBasisTime, */
/* 				      reduced->eBasisTime, st_basis_temp); */
/*     //------------------------------------------------------------------------- */
/*     // Calculate the ST-HOSVD Tailored space-time basis */
/*     //------------------------------------------------------------------------- */
/*     yk_createSpaceTimeBasis(&reduced->rOBState, st_basis_temp, */
/* 			    &reduced->ST_rOBState, reduced->ST_rOBState_i); */
/*     //------------------------------------------------------------------------- */
/*     // Add the space-time basis to the whole space-time matrix the ultimate one */
/*     //------------------------------------------------------------------------- */

/*     //------------------------------------------------------------------------- */
/*     // Calculate the initial conditions for each subspace and than append */
/*     //------------------------------------------------------------------------- */
/*     createInitialConditions_ST(ykflow, primalApprox, s_mu, &init, reduced); */
/*   } */
/*   //--------------------------------------------------------------------------- */
/*   // Divide up the ST-ROBSTATE and add to the rOBState_i thing for easy solves */
/*   //--------------------------------------------------------------------------- */
/*   multiquation->equationReduced.numStates = nBasis; */
/* } */