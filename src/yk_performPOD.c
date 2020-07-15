#include "yk_performPOD.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//check check it worked

Mat * yk_createSnapshotState(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Mat *snapshot, int snapshotStep,
			     Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k, p;                       //initialization for iteration
  int snapshot_i;                       //state to be included in snapshot
  int snapshotCount = primal->self->time.count;
  int totalParamSnapshot = snapshotCount*reduced->numParamSet;
  PetscInt index_i;
  PetscScalar *state_i =                //state array to add to snapshot matrix
    (PetscScalar *) malloc (primal->self->systemSize*sizeof(PetscScalar));
  PetscScalar *state_0 =                //state array to add to snapshot matrix
    (PetscScalar *) malloc (primal->self->systemSize*sizeof(PetscScalar));
  char cwd[1024];
  char cwd2[1024];
  char fomBuffer[500];
  Mat *snapshot_mu;
  PetscInt row, col;
  PetscInt *temp = (PetscInt *) malloc (primal->self->systemSize*sizeof(PetscInt));
  for (i=0; i<primal->self->systemSize; i++){
    temp[i] = i;
  }
  char tempS[1000];
  PetscMalloc1(reduced->numParamSet, &snapshot_mu);
  int timeNode0 = round(primal->self->time.t_0/primal->self->time.dt);
 struct stat sb = {0};

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
    getcwd(cwd2, sizeof(cwd2));

    //Retrieve the data set name and enter that directory
    ykflow->FomName(ykflow, multiquation->equation, reduced, p, fomBuffer);
    printf("Get Data from --->>>> %s\n", fomBuffer);

    if (stat(fomBuffer, &sb) == 0 && S_ISDIR(sb.st_mode)) {
      chdir(fomBuffer);
      getcwd(cwd, sizeof(cwd));
    } else {
      perror(fomBuffer);
      exit(0);
    }
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
    chdir(cwd2);
  }

  MatAssemblyBegin(*snapshot, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*snapshot, MAT_FINAL_ASSEMBLY);
  /* strcpy(primal->self->id, tempS); */
  //---------------------------------------------------------------------------
  // Termination
  //---------------------------------------------------------------------------
  free(state_0);
  free(state_i);
  free(temp);
  //Optional return of snapshot_mu for space time stuff
  return snapshot_mu;
}

/* void yk_createTemporalSnapshotState(yk_PrimalSolver *ykflow, */
/* 				    Multiverse *multiquation, Cluster *primal, */
/* 				    Mat spaceBasis, Mat *snapshot_phi, */
/* 				    Mat *s_mu, int snapshotStep, */
/* 				    Is_it *reduced){ */
/*   int i, j, p; //initialization for iteration */
/*   Vec state_0; */
/*   Vec state_i; */
/*   Vec phi; */
/*   PetscScalar val; */
/*   PetscInt row=primal->self->systemSize; */
/*   //Number of colum.ns in that snapshot matrix basically */
/*   int col = primal->self->time.count; */
/*   PetscInt m,n; */
/*   /\* int col = ceil(((primal->self->time.t_f-primal->self->time.t_0)*10)/ *\/ */
/*   /\*                          (primal->self->time.dt*10)); *\/ */
/* /\* PetscInt col = (primal->self->time.t_f-primal->self->time.t_0)/ *\/ */
/*   /\*   primal->self->time.dt; *\/ */
/*   /\* PetscInt col=primal->self->time.count; *\/ */
/*   char cwd[1024]; */
/*   char fomBuffer[50]; */
/*   PetscInt index_i; */
/*   PetscInt snapshot_i; */
/*   Vec psi_p; */
/*   PetscInt *rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt)); */
/*   PetscInt *colIndex; */
/*   for (i=0; i< row; i++) */
/*     rowIndex[i] = i; */
/*   PetscInt rbasisrow, rbasiscol; */
/*   MatGetSize(spaceBasis, &row, &rbasiscol); */
/*   //--------------------------------------------------------------------------- */
/*   // Initialization */
/*   //--------------------------------------------------------------------------- */
/*   for (i=0; i<rbasiscol; i++){ */
/*     MatCreate(PETSC_COMM_SELF, &snapshot_phi[i]);//create snapshot matrix */
/*     MatSetSizes(snapshot_phi[i], col, reduced->numParamSet, col, */
/* 		reduced->numParamSet); */
/*     MatSetType(snapshot_phi[i], MATSEQDENSE);        //set snapshot matrix type */
/*     MatSetUp(snapshot_phi[i]); */
/*   } */
/*   VecCreate(PETSC_COMM_SELF, &phi); */
/*   VecSetSizes(phi, primal->self->systemSize, primal->self->systemSize); */
/*   VecSetType(phi, VECSEQ); */

/* VecCreate(PETSC_COMM_SELF, &psi_p); */
/* VecSetSizes(psi_p, n, n); */
/* VecSetType(psi_p, VECSEQ); */
/*       for (i=0; i<n; i++) */
/* 	colIndex[i] = i; */

/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   for (j = 0; j<rbasiscol; j++){ */
/*     ks_MatCol2Vec(spaceBasis, row, rowIndex, phi, j, INSERT_VALUES); */
/*     for (p=0; p<reduced->numParamSet; p++){ */
/*       MatGetSize(s_mu[p], &m, &n); */
/*       colIndex = (PetscInt *) malloc (n*sizeof(PetscInt)); */


/*       MatMultTranspose(s_mu[p], phi, psi_p); */
/*       ks_Vec2MatCol(snapshot_phi[j], n, colIndex, p, psi_p, INSERT_VALUES); */
/*       VecDestroy(&psi_p); */
/*       free(colIndex); */
/*     } */
/*     MatAssemblyBegin(snapshot_phi[j], MAT_FINAL_ASSEMBLY); */
/*     MatAssemblyEnd(snapshot_phi[j], MAT_FINAL_ASSEMBLY); */
/*   } */
/*   VecDestroy(&phi); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   /\* ks_readSolution(multiquation->equation, primal->self, 0); *\/ */
/*   /\* array2Vec(multiquation->equation, primal->self, state_0); *\/ */
/*   /\* for (i=0; i<row; i++) *\/ */
/*   /\*   rowIndex[i] = i; *\/ */
/*   /\* for (j = 0; j<reduced->nBasisFuncs; j++){ *\/ */

/*   /\*   ks_MatCol2Vec(reduced->rOBState, row, rowIndex, phi, j, INSERT_VALUES); *\/ */
/*   /\*   for (p=0; p<reduced->numParamSet; p++){ *\/ */
/*   /\*     ykflow->FomName(ykflow, multiquation->equation, reduced, p, fomBuffer); *\/ */
/*   /\*     chdir(fomBuffer); *\/ */
/*   /\*     getcwd(cwd, sizeof(cwd)); *\/ */
/*   /\*     //Read the state solutions for the initial time in that directory *\/ */
/*   /\*     //Add data to the Petsc Mat Type *\/ */
/*   /\*     for (i=0; i<snapshotCount; i++){      //begin creating snapshot matrix *\/ */
/*   /\* 	index_i = p*snapshotCount+i; *\/ */
/*   /\* 	if (i==snapshotCount-1) *\/ */
/*   /\* 	  snapshot_i = primal->self->time.count; //add the last state *\/ */
/*   /\* 	else if (i==0) *\/ */
/*   /\* 	  snapshot_i = 1;                //do not include the initial condition *\/ */
/*   /\* 	else *\/ */
/*   /\* 	  snapshot_i = snapshotStep*i; *\/ */
/*   /\* 	//extract the state solution and add to the snapshot matrix *\/ */
/*   /\* 	ks_readSolution(multiquation->equation, primal->self, snapshot_i); *\/ */
/*   /\* 	array2Vec(multiquation->equation, primal->self, state_i); *\/ */
/*   /\* 	VecAXPY(state_i, -1, state_0); *\/ */
/*   /\* 	VecDot(state_i, phi, &val); *\/ */
/*   /\* 	MatSetValues(*snapshot, 1, &i, 1, &p, &val, INSERT_VALUES); *\/ */
/*   /\* 	/\\* MatSetValue(snapshot, i, p, val, INSERT_VALUES); *\\/ *\/ */
/*   /\*     } *\/ */
/*   /\*     chdir("../"); *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */
/*   /\* MatView(snapshot, PETSC_VIEWER_STDOUT_SELF); *\/ */
/*   /\* getchar(); *\/ */
/*   free(rowIndex); */
/*   free(colIndex); */
/* } */

void yk_createTemporalSnapshotState(yk_PrimalSolver *ykflow,
				    Multiverse *multiquation, Cluster *primal,
				    Mat spaceBasis, Mat *snapshot_phi,
				    Mat *s_mu, int snapshotStep, int dim3,
				    Is_it *reduced){
  /* 				    Is_it *reduced){ */
  int i, j, p; //initialization for iteration
  Vec state_0;
  Vec state_i;
  Vec phi;
  PetscScalar val;
  PetscInt row=primal->self->systemSize;
  //Number of colum.ns in that snapshot matrix basically
  int col = primal->self->time.count;
  PetscInt m,n;
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

  /* PetscInt *colIndex; */
  for (i=0; i< row; i++)
    rowIndex[i] = i;
  PetscInt rbasisrow, rbasiscol;
  MatGetSize(spaceBasis, &row, &rbasiscol);

  /* int i,j,p; */
  /* char cwd[1024]; */
  /* char fomBuffer[50]; */
  /* PetscInt index_i; */
  /* PetscInt snapshot_i; */
  /* Vec psi_p; */
  /* PetscInt *rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt)); */
  PetscInt *colIndex = (PetscInt *) malloc (col*sizeof(PetscInt));
  for (i=0; i< row; i++)
    rowIndex[i] = i;
  for (i=0; i<col; i++)
    colIndex[i] = i;
  /* PetscInt rbasisrow, rbasiscol; */
  /* MatGetSize(spaceBasis, &row, &rbasiscol); */
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  for (i=0; i<rbasiscol; i++){
    MatCreate(PETSC_COMM_SELF, &snapshot_phi[i]);//create snapshot matrix
    MatSetSizes(snapshot_phi[i], col, dim3, col, dim3);
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
  PetscInt mm, nn;
  MatGetSize(s_mu[0], &mm, &nn);


  for (j = 0; j<rbasiscol; j++){
    ks_MatCol2Vec(spaceBasis, row, rowIndex, phi, j, INSERT_VALUES);
    for (p=0; p< dim3; p++){
      MatMultTranspose(s_mu[p], phi, psi_p);
      ks_Vec2MatCol(snapshot_phi[j], col, colIndex, p, psi_p, INSERT_VALUES);
    }
    MatAssemblyBegin(snapshot_phi[j], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(snapshot_phi[j], MAT_FINAL_ASSEMBLY);
  }
  free(rowIndex);
  free(colIndex);
  VecDestroy(&psi_p);
  VecDestroy(&phi);
  VecDestroy(&psi_p);
}

void yk_properOrthogonalDecomposeGroup(Mat *snapshot, int numSnapshots,
				       int systemSize, PetscInt *index,
				       PetscInt *numSingularValues,
				       PetscScalar engyValues, Mat *A){
  int i;
  /* for (i=0; i<numSnapshots; i++){
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

      SVDGetSingularTriplet(svd, count, &sigma, singular, NULL);
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

/* void yk_properOrthogonalDecompose(Mat *snapshot, int systemSize, */
/* 				  PetscInt *index, PetscInt *numSingularValues, */
/* 				  PetscScalar engyValues, Mat *A){ */
/*   //------------------------------------//------------------------------------- */
/*   // Variables here                     // Comments section */
/*   //------------------------------------//------------------------------------- */
/*   int i, j;                             //initialization for iteration */
/*   int flag = 0; */
/*   double denomE = 0; */
/*   double Etotal = 0; */
/*   double EtotalPercent = 0; */
/*   PetscInt nsv; */
/*   PetscInt nconv; */
/*   PetscInt count = 0; */
/*   PetscInt m, n; */
/*   PetscReal sigma;                      //eigenvalues */
/*   Vec su; */
/*   Vec singular;                         //singular vectors rotation */
/*   Mat snapshotCopy; */
/*   Mat adjustedSnapshot; */
/*   SVD svd;                              //create svd object */
/*   //--------------------------------------------------------------------------- */
/*   // Initialization */
/*   //--------------------------------------------------------------------------- */
/*   MatGetSize(*snapshot, &m, &n); */
/*   yk_VecCreateSeq(&su, systemSize); */
/*   //Singular vector U */
/*   yk_VecCreateSeq(&singular, n); */
/*   //create a copy of snapshot matrix */
/*   yk_MatCreateSeqDense(&snapshotCopy, m, n); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   MatCopy(*snapshot, snapshotCopy, SAME_NONZERO_PATTERN); */
/*   MatTransposeMatMult(*snapshot, snapshotCopy , MAT_INITIAL_MATRIX, */
/* 		      PETSC_DEFAULT, &adjustedSnapshot); */
/*   //--------------------------------------------------------------------------- */
/*   // Create the Singular Value solver and set various options and solve */
/*   //--------------------------------------------------------------------------- */
/*   SVDCreate(PETSC_COMM_SELF, &svd); */
/*   SVDSetOperator(svd, adjustedSnapshot); */
/*   SVDSetFromOptions(svd); */
/*   SVDSolve(svd); */
/*   SVDGetDimensions(svd,&nsv,NULL,NULL); */
/*   SVDGetConverged(svd, &nconv); */
/*   //--------------------------------------------------------------------------- */
/*   // If nBasis is undetermined, than calculate witht he energy desire thing */
/*   //--------------------------------------------------------------------------- */
/*   if (engyValues != 0 && *numSingularValues == 0){ */
/*     for (i=0; i<nconv; i++){ */
/*       SVDGetSingularTriplet(svd, i, &sigma, NULL, NULL); */
/*       denomE+=sigma; */
/*     } */
/*     while (EtotalPercent < engyValues){ */
/*       SVDGetSingularTriplet(svd, count, &sigma, singular, NULL); */
/*       Etotal +=sigma; */
/*       EtotalPercent = Etotal/denomE; */
/*       count++; */
/*     } */
/*   } */
/*   //--------------------------------------------------------------------------- */
/*   // Create Basis Matrix here */
/*   //--------------------------------------------------------------------------- */
/*   yk_MatCreateSeqDense(A, systemSize, count); */

/*   for (i=0; i<count; i++){ */
/*     SVDGetSingularTriplet(svd, i, &sigma, singular, NULL); //extract rotation */
/*     MatMult(*snapshot, singular, su); */
/*     VecScale(su, 1.0/sqrt(sigma)); */
/*     ks_Vec2MatCol(*A, systemSize, index, i, su, INSERT_VALUES); */
/*   } */
/*   MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); */

/*   //--------------------------------------------------------------------------- */
/*   // Clean up */
/*   //--------------------------------------------------------------------------- */
/*   SVDDestroy(&svd); */
/*   VecDestroy(&su); */
/*   VecDestroy(&singular); */
/*   MatDestroy(&snapshotCopy); */
/*   MatDestroy(&adjustedSnapshot); */
/* } */


void yk_createSpaceTimeBasis(Mat *spacebasis, Mat *timebasis, Mat *spaceTime){

  int i, j;                                //initialization for iteration
  Vec phi;
  int init=0;
  PetscInt s_row, t_row;
  PetscInt s_col, t_col;
  PetscInt *rowIndex;
  PetscInt st_col = 0;
  PetscInt *rowstIndex;

  double **kronpiece=NULL;
  PetscInt *colstIndex = NULL;
  /*  PetscInt *colstIndex = (PetscInt *) malloc (colTime*sizeof(PetscInt)); */
 /*  PetscInt *colIndex = (PetscInt *) malloc (col*colTime*sizeof(PetscInt)); */
 /*  PetscInt *index_i = (PetscInt *) malloc (row*sizeof(PetscInt)); */

 /*  PetscScalar *st_basisA = (PetscScalar *) malloc */
 /*    (row*col*colTime*sizeof(PetscScalar)); */
  //---------------------------------------------------------------------------
  //Initialization
  //---------------------------------------------------------------------------
  MatGetSize(*spacebasis, &s_row, &s_col);
  rowIndex = (PetscInt *) malloc (s_row*sizeof(PetscInt));
  for (i=0; i< s_row; i++){rowIndex[i] = i;}


  //Vector representing the ith basis in the space basis matrix
  VecCreate(PETSC_COMM_SELF, &phi);
  VecSetSizes(phi, s_row, s_row);
  VecSetType(phi, VECSEQ);

  for (i=0; i<s_col; i++){
    MatGetSize(timebasis[i], &t_row, &t_col);
    st_col+=t_col;
  }
  rowstIndex = (PetscInt *) malloc (s_row*t_row*sizeof(PetscInt));

  for (i=0; i<s_row*t_row; i++){rowstIndex[i] = i;}

  MatCreate(PETSC_COMM_SELF, spaceTime);       //Set up the Matrix for POD
  MatSetSizes(*spaceTime, s_row*t_row, st_col, s_row*t_row, st_col);
  MatSetType(*spaceTime, MATSEQDENSE);
  MatSetUp(*spaceTime);

  //---------------------------------------------------------------------------
  //Implementation
  //---------------------------------------------------------------------------
  for (i = 0; i<s_col; i++){
    MatGetSize(timebasis[i], &t_row, &t_col);
    colstIndex = yk_reallocPetscInt(colstIndex, t_col);
    for (j=0; j<t_col; j++)
      colstIndex[j] = init+j;
    //Get the first basis function vector
    ks_MatCol2Vec(*spacebasis, s_row, rowIndex, phi, i, INSERT_VALUES);
    //Take the kronbecker product of the basis vector and the time basis
    yk_kron(timebasis[i], phi, &kronpiece);
    MatSetValues(*spaceTime, s_row*t_row, rowstIndex, t_col, colstIndex,
  		 *kronpiece, INSERT_VALUES);
    free(kronpiece[0]);
    free(kronpiece);
    init+=t_col;
  }
  MatAssemblyBegin(*spaceTime, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*spaceTime, MAT_FINAL_ASSEMBLY);


  VecDestroy(&phi);
  free(rowIndex);
  free(rowstIndex);
  free(colstIndex);
}

void yk_formatSpaceTimeBasis(Cluster *primal, Mat *subspaceTime,
			     Mat *spaceTime, Mat *spaceTime_i, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k, m;                       // initialization for iteration
  PetscInt j_n;
  PetscInt *nnz;
  PetscInt rows, nBasis;
  PetscInt n_rows, n_spaceTimeBasis=0;
  PetscInt *rIndex;
  PetscInt *bIndex;
  PetscInt *st_bIndex;
  PetscInt nBasis_pre=0;
  PetscInt b_0 = 0;
  PetscInt *st_index;
  PetscScalar *varray;
  PetscInt *colIndex;
  PetscInt sizeN = primal->self->systemSize;
  PetscScalar *vXarray;
  PetscInt *index_i = (PetscInt *) malloc (sizeN*sizeof(PetscInt));
  PetscInt *spaceindex = (PetscInt *) malloc (sizeN*sizeof(PetscInt));
  PetscScalar *ex_array;
  PetscInt n_b=0;

  //initialize the size of the spaceTime (number of rows in that big ass matrix
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatGetSize(subspaceTime[0], &rows, &nBasis);
  n_spaceTimeBasis = nBasis;

  rIndex = (PetscInt *) malloc (rows*sizeof(PetscInt));
  for (i=0; i<rows; i++)
    rIndex[i] = i;
  //For the off diagonals

  //Calcualte the total number of rows in the entire big ass matrix
  n_rows = rows*reduced->nSubWindows;
  nnz = (PetscInt *) malloc (n_rows*sizeof(PetscInt));

  //Find the total number rows and colums for the future spaceTime matrix
  for (i=1; i<reduced->nSubWindows; i++){
    MatGetSize(subspaceTime[i], &rows, &nBasis);
    n_spaceTimeBasis+=nBasis;
  }
  colIndex = (PetscInt *) malloc (n_spaceTimeBasis*sizeof(PetscInt));
  for (i=0; i<n_spaceTimeBasis; i++)
    colIndex[i] = i;
  for (i=0; i<reduced->nSubWindows; i++){
    MatGetSize(subspaceTime[i], &rows, &nBasis);
    n_b += nBasis;
    //Need to Matind the number of zeros for each row in the uber space-time
    for (j=0; j<rows; j++)
      nnz[i*rows+j] = n_b;
    for (j=0; j<primal->self->time.count; j++){
      j_n = i*primal->self->time.count+j;
      MatCreate(PETSC_COMM_SELF, &spaceTime_i[j_n]);
      MatSetSizes(spaceTime_i[j_n], primal->self->systemSize, n_spaceTimeBasis,
		  primal->self->systemSize, n_spaceTimeBasis);
      MatSetType(spaceTime_i[j_n], MATSEQAIJ);
      MatSeqAIJSetPreallocation(spaceTime_i[j_n], n_b, NULL);
      MatZeroEntries(spaceTime_i[j_n]);
    }
  }
  //Now that we have all the information, we can create the petsc MAT object
  MatCreate(PETSC_COMM_SELF, spaceTime);
  MatSetSizes(*spaceTime, n_rows, n_spaceTimeBasis, n_rows, n_spaceTimeBasis);
  MatSetType(*spaceTime, MATSEQAIJ); //Let's make it sparse
  MatSeqAIJSetPreallocation(*spaceTime, NULL, nnz);
  MatZeroEntries(*spaceTime);
 //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------

  //Now build the damn matrix, have to iterate through all the fucking windows
  for (i=0; i<reduced->nSubWindows; i++){
    //-------------------------------------------------------------------------
    // Prepare and initialize the data and such
    //-------------------------------------------------------------------------
    //Having to call MatSize all the time is so fucking annoying
    MatGetSize(subspaceTime[i], &rows, &nBasis);
    st_index = (PetscInt *) malloc (rows*sizeof(PetscInt));
    ex_array = (PetscScalar *) malloc
      (primal->self->systemSize*nBasis*sizeof(PetscScalar));
    //Need to properally allocate since each subwindow won't nec have samebasis
    bIndex = (PetscInt *) malloc (nBasis*sizeof(PetscInt));
    st_bIndex = (PetscInt *) malloc (nBasis*sizeof(PetscInt));
    varray = (PetscScalar *) malloc (nBasis*rows*sizeof(PetscScalar));

    for (j=0; j<nBasis; j++){
      bIndex[j] = j;
      st_bIndex[j] = b_0+j;
    }
    for (j=0; j<rows; j++)
      st_index[j] = i*rows+j;
    //-------------------------------------------------------------------------
    // Create the giant fuckng matrix
    //-------------------------------------------------------------------------
    //Set the diagonals here
    MatGetValues(subspaceTime[i], rows, rIndex, nBasis, bIndex, varray);
    MatSetValues(*spaceTime, rows, st_index, nBasis, st_bIndex, varray,
    		 INSERT_VALUES);

    for (j=0; j<primal->self->time.count; j++){
      for (k=0; k<primal->self->systemSize; k++)
 	index_i[k] = j*primal->self->systemSize+k;
      MatGetValues(subspaceTime[i], primal->self->systemSize, index_i,
 		   nBasis, bIndex, ex_array);
      MatSetValues(spaceTime_i[i*primal->self->time.count+j],
 		   primal->self->systemSize, primal->self->index, nBasis,
		   st_bIndex, ex_array,INSERT_VALUES);
    }
    /* //Set the lower off diagonal stuff here */
    for (j=0; j<primal->self->systemSize; j++)
      spaceindex[j] = rows-sizeN+j;

    MatGetValues(subspaceTime[i], primal->self->systemSize, spaceindex,
    		 nBasis, bIndex, ex_array);

    for (j=i+1; j<reduced->nSubWindows; j++){
      for (m=0; m<primal->self->time.count; m++){
    	for (k=0; k<primal->self->systemSize; k++)
    	  index_i[k]=j*(primal->self->time.count*sizeN)+m*sizeN+k;

    	MatSetValues(*spaceTime, primal->self->systemSize, index_i,
    		     nBasis, st_bIndex, ex_array, INSERT_VALUES);
    	MatSetValues(spaceTime_i[j*primal->self->time.count+m],
    		     primal->self->systemSize, primal->self->index, nBasis,
    		     st_bIndex,
    		     ex_array, INSERT_VALUES);
      }
    }
    for (j=0; j<primal->self->time.count; j++){
      MatAssemblyBegin(spaceTime_i[i*primal->self->time.count+j],
		       MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(spaceTime_i[i*primal->self->time.count+j],
		     MAT_FINAL_ASSEMBLY);
    }
    b_0+=nBasis;
    nBasis_pre = nBasis;
    free(ex_array);
    free(st_index);
    free(st_bIndex);
    free(bIndex);
    free(varray);
  }
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatAssemblyBegin(*spaceTime, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*spaceTime,MAT_FINAL_ASSEMBLY);
  free(spaceindex);
  free(colIndex);
  free(index_i);
  free(rIndex);
  free(nnz);
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
