#include "ks_greedyAlgorithm.h"

PetscInt *yk_reallocPetscInt(PetscInt *oldArray, int newSize){
  PetscInt *newArray;
  if (oldArray == NULL)
    newArray = (PetscInt *) malloc (newSize*sizeof(PetscInt));
  else
    newArray = (PetscInt *) realloc (oldArray, newSize*sizeof(PetscInt));
  if (newArray == NULL){
    printf("ERROR: error in reallocating the array");
    return NULL;
  }
  return newArray;
}

PetscScalar *yk_reallocPetscScalar(PetscScalar *oldArray, int newSize){
  PetscScalar *newArray;
  if (oldArray == NULL)
    newArray = (PetscScalar *) malloc (newSize*sizeof(PetscScalar));
  else
    newArray = (PetscScalar *) realloc (oldArray, newSize*sizeof(PetscScalar));
  if (newArray == NULL){
    printf("ERROR: error in reallocating the array");
    return NULL;
  }
  return newArray;
}


void yk_spaceTimeGappyError(Cluster *primal, int *numTempSam_r, Mat Z,
			    int count, Vec gappyError, PetscInt r_i, Is_it *reduced,
			    int dimenFlag){
  int i, j; //initialization for iteration
  PetscInt jj;
  PetscInt r_ilong;
  PetscInt colSize;
  PetscInt rowSize;
  PetscInt row, col, nbar_r;
  PetscInt *rowIndex = NULL;
  PetscScalar valnorm;
  Vec epsilon, epsilon0, epsilon1;
  Vec epsilonVec;
  Mat epsilonMat, epsilonMatT;
  Mat C, D, E, Cplus;
  Mat annoyBasis_r;
  Mat Identity;
  Vec phi_i;
  PetscInt *nbar_i =  NULL;
  PetscInt *nbarR = NULL;
  PetscInt m, n;
  PetscInt sysSize = primal->self->systemSize;
  PetscInt sysTime = primal->self->time.count; //not include t0
  PetscInt *index_i = (PetscInt *) malloc (sysTime *sizeof(PetscInt));
  PetscScalar *st_basis_i = NULL;
  PetscScalar *phi_error=NULL;
  PetscScalar *phi_e=NULL;
  //---------------------------------------------------------------------------
  // Initlization
  //---------------------------------------------------------------------------
  MatGetSize(reduced->ST_rOBResidual, &row, &nbar_r);

  if (dimenFlag == 1){ //Space
    colSize = sysSize;
    rowSize = sysTime;
  }else{
    colSize = sysTime;
    rowSize = sysSize;
  }

  rowIndex = (PetscInt *) malloc (row*sizeof(PetscInt));
  for (i=0; i<row; i++)
    rowIndex[i] = i;

  for (i=0; i<sysTime; i++)
    index_i[i] = i;

  st_basis_i = (PetscScalar *) malloc (row*sizeof(PetscScalar));

  VecCreate(PETSC_COMM_SELF, &epsilonVec);
  VecSetSizes(epsilonVec, rowSize, rowSize);
  VecSetFromOptions(epsilonVec);

  VecCreate(PETSC_COMM_SELF, &phi_i);
  VecSetSizes(phi_i, row, row);
  VecSetFromOptions(phi_i);

    //Create the identity matrix for the whole space-time
  MatCreate(PETSC_COMM_SELF, &Identity);
  MatSetSizes(Identity, sysTime*sysSize, sysTime*sysSize, sysTime*sysSize,
	      sysSize*sysTime);
  MatSetType(Identity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(Identity, 1, NULL);
  MatZeroEntries(Identity);
  MatAssemblyBegin(Identity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Identity, MAT_FINAL_ASSEMBLY);
  MatShift(Identity, 1);

   //Begin Greedy Iteration
  MatCreate(PETSC_COMM_SELF, &epsilonMat);
  MatSetSizes(epsilonMat, sysSize, sysTime, sysSize, sysTime);
  MatSetType(epsilonMat, MATSEQDENSE);
  MatSetUp(epsilonMat);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  if (numTempSam_r[r_i] != 0){
    if (r_i==0){
      MatGetValues(reduced->ST_rOBResidual, row, rowIndex, 1, &r_i, st_basis_i);
    }else{
      MatGetSize(Z, &m, &n);
      //Find the partial phi
      phi_error = (PetscScalar *) realloc (phi_e, r_i*row*sizeof(PetscScalar));
      //Used to calculate the parts of the error vector
      MatCreate(PETSC_COMM_SELF, &annoyBasis_r);
      MatSetSizes(annoyBasis_r, row, r_i, row, r_i);
      MatSetType(annoyBasis_r, MATSEQDENSE);
      MatSetFromOptions(annoyBasis_r);
      MatSetUp(annoyBasis_r);

      //Create and fill up the annoyBasis_r [phi_r,1 .... phi_r, i-1]
      //Number of column basis to look, create index
      nbar_i = (PetscInt *) realloc (nbarR, r_i *sizeof(PetscInt));
      for (j=0; j<r_i; j++)
	nbar_i[j] = j;
      r_ilong = r_i;
      MatGetValues(reduced->ST_rOBResidual,row, rowIndex, r_i, nbar_i, phi_error);

      MatSetValues(annoyBasis_r, row, rowIndex, r_i, nbar_i, phi_error,
		   INSERT_VALUES);
      MatAssemblyBegin(annoyBasis_r, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(annoyBasis_r, MAT_FINAL_ASSEMBLY);


      //calculate the error in the gappy POS approximation of phi_i (epsilon)
      //Z*[phi_r,1 .... phi_r, i-1]^+*Z
      MatMatMult(Z, annoyBasis_r, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
      moorePenrosePseudoInv(C, m, r_i, &Cplus);
      /* MatView(C, PETSC_VIEWER_STDOUT_SELF); */
      /* getchar(); */

      MatMatMult(Cplus, Z, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D);

      VecCreate(PETSC_COMM_SELF, &epsilon);
      VecSetSizes(epsilon, row, row);
      VecSetFromOptions(epsilon);

      VecCreate(PETSC_COMM_SELF, &epsilon0);
      VecSetSizes(epsilon0, r_i, r_i);
      VecSetFromOptions(epsilon0);

      VecCreate(PETSC_COMM_SELF, &epsilon1);
      VecSetSizes(epsilon1, row, row);
      VecSetFromOptions(epsilon1);

      MatGetColumnVector(reduced->ST_rOBResidual, phi_i, r_i);

      MatMult(D, phi_i, epsilon0);
      MatMult(annoyBasis_r, epsilon0, epsilon1);
      VecScale(epsilon1, -1);
      MatMultAdd(Identity, phi_i, epsilon1, epsilon);

      VecGetValues(epsilon, row, rowIndex, st_basis_i);
      MatDestroy(&annoyBasis_r);
      VecDestroy(&epsilon);
      VecDestroy(&epsilon0);
      VecDestroy(&epsilon1);
      MatDestroy(&D);
      MatDestroy(&C);
      MatDestroy(&Z);
      MatDestroy(&Cplus);
    }

    MatSetValues(epsilonMat, sysSize, primal->self->index, sysTime, index_i,
		 st_basis_i, INSERT_VALUES);
    MatAssemblyBegin(epsilonMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(epsilonMat, MAT_FINAL_ASSEMBLY);
    if (dimenFlag == 1) //Space
      MatTranspose(epsilonMat, MAT_INITIAL_MATRIX, &epsilonMatT);
    else
      MatConvert(epsilonMat, MATSAME, MAT_INITIAL_MATRIX, &epsilonMatT);
    //Calculate the epsilon error to use to find the indicies
    for (j=0; j<colSize; j++){
      jj = j;
      MatGetColumnVector(epsilonMatT, epsilonVec, jj);
      VecNorm(epsilonVec, NORM_2, &valnorm);
      VecSetValues(gappyError, 1, &jj, &valnorm, INSERT_VALUES);
    }
    VecAssemblyBegin(gappyError);
    VecAssemblyEnd(gappyError);
    MatDestroy(&epsilonMatT);

  }
  // Need to add information about the temporal shit in reduced


  free(rowIndex);
  free(index_i);
  free(nbar_i);
  free(st_basis_i);
  VecDestroy(&phi_i);
  MatDestroy(&Identity);
  MatDestroy(&epsilonMat);
  free(phi_error);
  VecDestroy(&epsilonVec);

}


void yk_greedyAlgorithm_spatialSet(Cluster *primal, PetscInt **spatialSet,
                                   PetscInt **temporalSet, Mat *Z,
				   Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k;                           //initialization for iteration
  PetscInt count0=0, count =0;               //Number of temporal indices keeping
  int flag;
  PetscInt timeCount;
  int numBasis = primal->self->basis.nodes;
  PetscInt *numTempSam_r;
  PetscInt sysSize = primal->self->systemSize;
  PetscInt sysTime = primal->self->time.count; //not include t0
  PetscInt *keepTrack = (PetscInt *) malloc (primal->self->space.elem.count*sizeof(PetscInt));
  PetscInt row, col, nbar_r;
  PetscInt locMax;
  PetscInt nbar_s= reduced->nSampleNodes; //Force the beginning time only
  PetscScalar maxValue;
  Vec temporal_samSet;                  //temporal sample set
  Vec gappyError;
  Mat spaceIdentity,  timeIdentity;
  Mat Zin;
  PetscInt one = 1;
  int *tempNodeUSet = (int *) malloc     // array info on which element to keep
    (primal->self->space.elem.count*sizeof(int));
  int countwhat = 0;
  int countyuki = 0;;
  PetscInt *spatialFullSet=NULL;
  PetscInt left; //everything is in terms of elements now

  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  *spatialSet = NULL;
  for (i=0; i<primal->self->space.elem.count; i++)
    tempNodeUSet[i] = 0;
  for (i=0; i<primal->self->space.elem.count; i++)
    keepTrack[i] = 0;
  MatGetSize(reduced->ST_rOBResidual, &row, &nbar_r);

  numTempSam_r = (PetscInt *) malloc (nbar_r*sizeof(PetscInt));

  VecCreate(PETSC_COMM_SELF, &gappyError);
  VecSetSizes(gappyError, sysSize, sysSize);
  VecSetFromOptions(gappyError);

  //Create the idenity spatial matrix
  timeCount = reduced->nSampleTime-1;
  MatCreate(PETSC_COMM_SELF, &timeIdentity);
  MatSetSizes(timeIdentity, timeCount, sysTime, timeCount, sysTime);
  MatSetType(timeIdentity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(timeIdentity, 1, NULL);
  //Fill in the identity matrix based on the time indices we want to keep
  for (j=0; j<timeCount; j++)
    MatSetValue(timeIdentity, j, (*temporalSet)[j], 1, INSERT_VALUES);
  MatAssemblyBegin(timeIdentity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(timeIdentity, MAT_FINAL_ASSEMBLY);

  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Determine number of temporal samples to compute at each greedy iteration
  for (i=0; i<nbar_r; i++){
    if (i<nbar_s % nbar_r)
      numTempSam_r[i] = floor(nbar_s/nbar_r)+1;
    else
      numTempSam_r[i] = floor(nbar_s/nbar_r);
  }
  /* //Iterate through the number of residual basis */
  for (i=0; i<nbar_r; i++){
    if (numTempSam_r[i] > 0){
      count0 = count;  //Save the number of time nodes from previous iteration
      count += numTempSam_r[i];
      yk_spaceTimeGappyError(primal,numTempSam_r, Zin, count, gappyError,i,
			     reduced, 1);
      spatialFullSet = yk_reallocPetscInt(spatialFullSet, count*12);
      *spatialSet = yk_reallocPetscInt(*spatialSet,count);
      for (j=0; j<numTempSam_r[i]; j++){
	flag = 0;
	while (flag == 0){
	  /* 	printf("wo\n"); */
	  VecMax(gappyError, &locMax, &maxValue);
	  left = (int)(locMax/12.0);
	  if (keepTrack[left] ==0){
	    keepTrack[left]  = 1;
	    (*spatialSet)[countwhat]=left;
	    tempNodeUSet[left] = 1;
	    for (k=0; k<12; k++){
	      spatialFullSet[countyuki] = 12*left+k;
	      countyuki++;
	    }
	    countwhat++;
	    flag =1 ;
	  }
	  for (k=0; k<12; k++)
	    VecSetValue(gappyError, 12*left+k, 0, INSERT_VALUES);
	}
      }
      PetscSortInt(count*12, spatialFullSet);
      PetscSortInt(count, *spatialSet);

      /*   //Create the identity matrix for time adfter adding more nodes */
      MatCreate(PETSC_COMM_SELF, &spaceIdentity);
      MatSetSizes(spaceIdentity, count*12, sysSize, count*12, sysSize);
      MatSetType(spaceIdentity, MATSEQAIJ);
      MatSeqAIJSetPreallocation(spaceIdentity, 1, NULL);
      /*   //Fill in the identity matrix based on the time indices we want to keep */
      for (j=0; j<count*12; j++)
	MatSetValue(spaceIdentity, j, spatialFullSet[j], 1, INSERT_VALUES);
      MatAssemblyBegin(spaceIdentity, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(spaceIdentity, MAT_FINAL_ASSEMBLY);
      /*   //Create the Z matrix given the time and space indces */
      yk_kkron(timeIdentity, spaceIdentity, &Zin);
      MatDestroy(&spaceIdentity);
    }
  }
  MatDuplicate(Zin, MAT_COPY_VALUES, Z);
  //---------------------------------------------------------------------------
  // Create the reducedMesh with elems that are the sample elems
  //---------------------------------------------------------------------------
  reduced->reducedMesh.elem.count = reduced->nSampleNodes;
  reduced->reducedMesh.node.count = reduced->nSampleNodes*numBasis;
  initArrayInt(&reduced->reducedMesh.elem, reduced->nSampleNodes);
  initArrayInt(&reduced->reducedMesh.node, reduced->nSampleNodes*numBasis);
  count = 0;


  for (i=0; i<primal->self->space.elem.count; i++){
    if (tempNodeUSet[i] == 1){ //If this element exists in reduced mesh
      reduced->reducedMesh.elem.array[count] = i;
      for (j=0; j<numBasis; j++)
        reduced->reducedMesh.node.array[count*numBasis+j] = i*numBasis+j;
      count++;
    }
  }
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  free(numTempSam_r);
  free(keepTrack);
  free(spatialFullSet);
  VecDestroy(&gappyError);
  MatDestroy(&timeIdentity);
  MatDestroy(&Zin);
  free(tempNodeUSet);
}



void yk_greedyAlgorithm_temporalSet(Cluster *primal, PetscInt **temporalSet,
				    Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             //initialization for iteration
  PetscInt count0=0, count =0;               //Number of temporal indices keeping
  int flag;
  PetscInt *numTempSam_r;
  PetscInt sysSize = primal->self->systemSize;
  PetscInt sysTime = primal->self->time.count; //not include t0
  PetscInt *keepTrack = (PetscInt *) malloc ((sysTime)*sizeof(PetscInt));
  PetscInt row, col, nbar_r;
  PetscInt locMax;
  PetscInt nbar_t= reduced->nSampleTime-1; //Force the beginning time only
  PetscScalar maxValue;
  Vec temporal_samSet;                  //temporal sample set
  Vec gappyError;
  Mat Z;
  Mat spaceIdentity,  timeIdentity;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  *temporalSet = NULL;

  for (i=0; i<sysTime; i++)
    keepTrack[i] = 0;

  MatGetSize(reduced->ST_rOBResidual, &row, &nbar_r);
  numTempSam_r = (PetscInt*) malloc (nbar_r*sizeof(PetscInt));

  VecCreate(PETSC_COMM_SELF, &gappyError);
  VecSetSizes(gappyError, sysTime, sysTime);
  VecSetFromOptions(gappyError);

  //Create the idenity spatial matrix
  MatCreate(PETSC_COMM_SELF, &spaceIdentity);
  MatSetSizes(spaceIdentity, sysSize, sysSize, sysSize, sysSize);
  MatSetType(spaceIdentity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(spaceIdentity, 1, NULL);
  MatZeroEntries(spaceIdentity);
  MatAssemblyBegin(spaceIdentity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(spaceIdentity, MAT_FINAL_ASSEMBLY);
  MatShift(spaceIdentity, 1);
 //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Determine number of temporal samples to compute at each greedy iteration
  for (i=0; i<nbar_r; i++){
    if (i<nbar_t % nbar_r)
      numTempSam_r[i] = floor(nbar_t/nbar_r)+1;
    else
      numTempSam_r[i] = floor(nbar_t/nbar_r);
  }

  //Begin the greedy iteration
  //Iterate through the number of residual basis
  for (i=0; i<nbar_r; i++){
    if (numTempSam_r[i] > 0){
      yk_spaceTimeGappyError(primal,numTempSam_r, Z, count, gappyError,i,
			     reduced, 0);
      //Calculate the total number of indicies we are adding here
      count0 = count;  //Save the number of time nodes from previous iteration
      count += numTempSam_r[i];

      //change size of the time indices to take to account of new indices added
      *temporalSet = yk_reallocPetscInt(*temporalSet, count);

      //Find the time indicies to add
      for (j=0; j<numTempSam_r[i]; j++){
	flag = 0;
	while (flag == 0){
	  VecMax(gappyError, &locMax, &maxValue);
	  if (keepTrack[locMax] ==0){
	    keepTrack[locMax]  = 1;
	    (*temporalSet)[count0+j] = locMax;
	    flag = 1;
	  }
	  VecSetValue(gappyError, locMax, 0, INSERT_VALUES);
	}
      }
      //Sort the values so that you can create the Z matrix correctly
      PetscSortInt(count, *temporalSet);

      //Create the identity matrix for time adfter adding more nodes
      MatCreate(PETSC_COMM_SELF, &timeIdentity);
      MatSetSizes(timeIdentity, count, sysTime, count, sysTime);
      MatSetType(timeIdentity, MATSEQAIJ);
      MatSeqAIJSetPreallocation(timeIdentity, 1, NULL);

      //Fill in the identity matrix based on the time indices we want to keep
      for (j=0; j<count; j++)
	MatSetValue(timeIdentity, j, (*temporalSet)[j], 1, INSERT_VALUES);
      MatAssemblyBegin(timeIdentity, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(timeIdentity, MAT_FINAL_ASSEMBLY);
      //Create the Z matrix given the time and space indces
      yk_kkron(timeIdentity, spaceIdentity, &Z);
      MatDestroy(&timeIdentity);
    }
  }


  initArrayInt(&reduced->reducedTime, reduced->nSampleTime-1);
  for (i=0; i<reduced->nSampleTime-1; i++)
    reduced->reducedTime.array[i] = (*temporalSet)[i];
  for (i=0; i<sysTime; i++){
    if (keepTrack[i] ==1 && i>0)
      keepTrack[i-1] = 1;
    if (keepTrack[i] ==1 && i>1)
      keepTrack[i-2] = 1;
  }
  int count_Os = 0;
  for (i=0; i<sysTime; i++)
    count_Os +=keepTrack[i];
  int count_i = 0;
  initArrayInt(&reduced->reducedTime_Os, count_Os);
  for (i=0; i<sysTime; i++){
    if (keepTrack[i]==1){
      reduced->reducedTime_Os.array[count_i] = i;
      count_i++;
    }
  }

  free(numTempSam_r);
  free(keepTrack);
  VecDestroy(&gappyError);
  MatDestroy(&spaceIdentity);
  MatDestroy(&Z);
}


/* void yk_greedyAlgorithmST_space(Cluster *primal,int seedCount,Is_it *reduced){ */
/*   //Is_It *reduced contains information on the residual basis */
/*   //reduced->nSampleNodes */
/*   //Start with the temporal sample set all to 1. Initialize to all ones. */
/*   int nAddNodes = reduced->nSampleNodes-seedCount; */
/*   int spatialNodes = primal->self->solution[0].nodes.count; */
/*   Vec sptialSampleSet; */
/*   Vec errorVec; */
/*   PetscInt row, col; */
/*   /\* int *spatialSampleSet = (int *) malloc (spatialNodes*sizeof(int)); *\/ */
/*   int *numSpaceSamples = (int *) malloc (reduced->nBasisFuncsRJ*sizeof(int)); */
/*   //--------------------------------------------------------------------------- */
/*   // Initialization */
/*   //--------------------------------------------------------------------------- */
/*   MatGetSize(reduced->rOBResidual, &row, &col); */


/*   VecCreate(PETSC_COMM_SELF, &sptialSampleSet); */
/*   VecSetSizes(spatialSampleSet, spatialNodes, spatialNodes); */
/*   VecSetType(spatialSampleSet, VECSEQ); */

/*   VecCreate(PETSC_COMM_SELF, &errorVec); */
/*   VecSetSizes(errorVec, row, row); */
/*   VecSetType(errorVec, VECSEQ); */

/*   VecSetEntires(spatialSampleSet, 0); */
/*   //reduced->nBasisFuncsRJ is the number of greedy iterations */
/*   //Determine number of spatial samples to compute at each greedy iterration */
/*   for (i=0; i<nAddNodes % reduced->nBasisFuncsRJ; i++) */
/*     numSpaceSamples[i] = floor(nAddNodes/reduced->nBasisFuncsRJ)+1; */
/*   for (i=nAddNodes % reduced->nBasisFuncsRJ; i<reduced->nBasisFuncsRJ; i++) */
/*     numSpaceSamples[i] = floor(nAddNodes/reduced->nBasisFuncsRJ); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   // Peform Greedy Iteration here */
/*   for (i=0; i<reduced->nBasisFuncsRJ; i++){ */
/*     if (i==0) */
/*       ks_MatCol2Vec(reduced->rOBResidual, row, rowIndex, errorVec, i, */
/* 		    INSERT_VALUES); */
/*     else */
/*       Dotheother stuff; */
/*     for (j=0; j<numSpaceSamples[i]; j++){ */

/*     } */
/*   } */
/* } */



void ks_greedyAlgorithm(int seedCount, PetscInt *nodeSet, Is_it *reduced){
  int i, j;                            //initialization for iteration
  int nAddNodes = reduced->nSampleNodes-seedCount; //num of additional nodes
  int nb = 0;               //num working basis vectors used: coute
  int nGIterations=minInt(reduced->nBasisFuncsRJ, nAddNodes);
  int nRHS = ceil((double) reduced->nBasisFuncsRJ/nAddNodes);  //number of r
  int nCi;
  int nCiMin = floor((double)reduced->nBasisFuncsRJ/nGIterations);
  int nAi;
  int nAiMin = floor((double) nAddNodes*nRHS/reduced->nBasisFuncsRJ);
  int count = seedCount;               //counter to keep track nodeSetValues
  int *tempNodeSet;                    //Used to avoid duplicates
  int *rowIndex;                        //Contains index values for n (PETSC)
  int *nCiIndex;                       //Contains index values for Ci (PETSC
  PetscInt row, col;                   //Size of the system
  PetscInt maxErrorLoc;                //location of the error
  PetscScalar maxR;                    //Max error value from esidual
  PetscScalar maxJ;                    //Max error value fromJacobian
  PetscScalar nError;                  //maxR + maxJ
  PetscReal maxErrorValue;             //Contains the value
  Vec rOBVec;                          //A general vector ize
  Vec resJacTemp;                      //Size nCi
  Vec resHatvec;                       //size count
  Vec greekVec;                        //size nb
  Vec phiQ;                            //ny
  Vec nErrorVec;                       //where nErroris saved
  Mat iZ;                              //Indentity mtrix is saved here
  Mat Residual;                        //Residual Mtrix for greedy
  Mat Jacobian;                        //Jacobian atrix for greedy
  Mat rOBResidualnb;                   //Intermedate values of ROB
  Mat rOBJacobiannb;                   //Intermeiate values of ROB
  Mat resHat;                          //Intermdiate values of ROB
   //--------------------------------------------------------------------------
  // Initialization
  //--------------------------------------------------------------------------
  MatGetSize(reduced->rOBResidual, &row, &col);
  tempNodeSet = (int *) malloc (row*sizeof(int));
  rowIndex = (int *) malloc (row*sizeof(int));
  for (i=0; i<row; i++)
    rowIndex[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, row, &rOBVec);
  VecCreateSeq(PETSC_COMM_SELF, row, &nErrorVec);
  VecCreateSeq(PETSC_COMM_SELF, row, &phiQ);
  for (i=0; i<row; i++)
    tempNodeSet[i] = 0;
  tempNodeSet[0] = 1;
  tempNodeSet[row-1] = 1;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<nGIterations; i++){
    nCi = nCiMin;
    nAi = nAiMin;
    if (i+1 <= reduced->nBasisFuncsRJ % nGIterations)
      nCi++;
    if (nRHS == 1 && i+1<= nAddNodes % reduced->nBasisFuncsRJ)
      nAi++;
    MatCreateSeqDense(PETSC_COMM_SELF, row, nCi, NULL, &Residual);
    MatCreateSeqDense(PETSC_COMM_SELF, row, nCi, NULL, &Jacobian);
    VecCreateSeq(PETSC_COMM_SELF, nCi, &resJacTemp);
    nCiIndex = (int *) malloc (nCi*sizeof(int));
    for (j=0; j<nCi; j++)
      nCiIndex[j] = j;
    //-------------------------------------------------------------------------
    // Creating the Z identity Matrix
    //-------------------------------------------------------------------------
    if (i==0){
      for (j=0; j<nCi; j++){ //rOBResidualTemp has nCi columns
        ks_MatCol2Vec(reduced->rOBResidual, row, rowIndex, rOBVec, j,
		      INSERT_VALUES);
        ks_Vec2MatCol(Residual, row, rowIndex, j, rOBVec, INSERT_VALUES);
        ks_MatCol2Vec(reduced->rOBJacobian, row, rowIndex, rOBVec, j,
		      INSERT_VALUES);
        ks_Vec2MatCol(Jacobian, row, rowIndex, j, rOBVec, INSERT_VALUES);
      }
    }else{
      MatCreateSeqAIJ(PETSC_COMM_SELF, count, row, 1, NULL, &iZ);
      for (j=0; j<count; j++)
        MatSetValue(iZ, j, nodeSet[j], 1, INSERT_VALUES);
      MatAssemblyBegin(iZ, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(iZ, MAT_FINAL_ASSEMBLY);
      MatCreateSeqDense(PETSC_COMM_SELF, row, nb, NULL, &rOBResidualnb);
      MatCreateSeqDense(PETSC_COMM_SELF, row, nb, NULL, &rOBJacobiannb);
      MatCreateSeqDense(PETSC_COMM_SELF, count, nb, NULL, &resHat);
      VecCreateSeq(PETSC_COMM_SELF, count, &resHatvec);
      VecCreateSeq(PETSC_COMM_SELF, nb, &greekVec);
      for (j=0; j<nb; j++){
        ks_MatCol2Vec(reduced->rOBResidual, row, rowIndex, rOBVec, j,
		      INSERT_VALUES);
        ks_Vec2MatCol(rOBResidualnb, row, rowIndex, j, rOBVec, INSERT_VALUES);
        ks_MatCol2Vec(reduced->rOBJacobian, row, rowIndex, rOBVec, j,
		      INSERT_VALUES);
        ks_Vec2MatCol(rOBJacobiannb, row, rowIndex, j, rOBVec, INSERT_VALUES);
      }
      MatAssemblyBegin(rOBResidualnb, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(rOBResidualnb, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(rOBJacobiannb, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(rOBJacobiannb, MAT_FINAL_ASSEMBLY);
      for (j=0; j<nCi; j++){
        //---------------------------------------------------------------------
        // Residual Thing
        //--------------------------------------------------------------------
        ks_MatCol2Vec(reduced->rOBResidual, row, rowIndex, rOBVec, nb+j,
		      INSERT_VALUES);
        MatMatMult(iZ, rOBResidualnb, MAT_REUSE_MATRIX, PETSC_DEFAULT, &resHat);
        MatMult(iZ, rOBVec, resHatvec);
        linearLeastSquares(resHat, resHatvec, greekVec);
        MatMult(rOBResidualnb, greekVec, phiQ);
        VecAXPY(rOBVec, -1, phiQ);
        ks_Vec2MatCol(Residual, row, rowIndex, j, rOBVec, INSERT_VALUES);
	//---------------------------------------------------------------------
        // Jacobian Thing
        //--------------------------------------------------------------------
        ks_MatCol2Vec(reduced->rOBJacobian, row, rowIndex, rOBVec, nb+j,
		      INSERT_VALUES);
        MatMatMult(iZ, rOBJacobiannb, MAT_REUSE_MATRIX, PETSC_DEFAULT, &resHat);
        MatMult(iZ, rOBVec, resHatvec);
        linearLeastSquares(resHat, resHatvec, greekVec);
        MatMult(rOBJacobiannb, greekVec, phiQ);
        VecAXPY(rOBVec, -1, phiQ);
        ks_Vec2MatCol(Jacobian, row, rowIndex, j, rOBVec, INSERT_VALUES);
      }
      VecDestroy(&resHatvec);
      VecDestroy(&greekVec);
      MatDestroy(&iZ);
      MatDestroy(&resHat);
      MatDestroy(&rOBResidualnb);
      MatDestroy(&rOBJacobiannb);
    }
    MatAssemblyBegin(Residual, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Residual, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Jacobian, MAT_FINAL_ASSEMBLY);
    for (j=0; j<row; j++){
      ks_MatRow2Vec(Residual, nCi, nCiIndex, resJacTemp, j, INSERT_VALUES);
      VecDot(resJacTemp, resJacTemp, &maxR);
      ks_MatRow2Vec(Jacobian, nCi, nCiIndex, resJacTemp, j, INSERT_VALUES);
      VecDot(resJacTemp, resJacTemp, &maxJ);
      nError = maxR + maxJ;
      VecSetValue(nErrorVec, j, nError, INSERT_VALUES);
    }
    MatDestroy(&Residual);
    MatDestroy(&Jacobian);
    VecDestroy(&resJacTemp);
    free(nCiIndex);
    VecAssemblyBegin(nErrorVec);
    VecAssemblyEnd(nErrorVec);
    for (j=0; j<nAi; j++){
      do{
        VecMax(nErrorVec, &maxErrorLoc, &maxErrorValue);
        tempNodeSet[maxErrorLoc]++;
        VecSetValue(nErrorVec, maxErrorLoc, 0, INSERT_VALUES);
      }while(tempNodeSet[maxErrorLoc]>1);
      nodeSet[count] = maxErrorLoc;
      count++;
    }
    partialSortInt(nodeSet, 0, count);
    nb += nCi;
  }
  free(tempNodeSet);
  free(rowIndex);
  VecDestroy(&phiQ);
  VecDestroy(&rOBVec);
  VecDestroy(&nErrorVec);
}

void partialSortInt(PetscInt array[], int first, int last){
  int i;
  PetscInt *inew = (PetscInt *) malloc ((last-first)*sizeof(PetscInt));;
  for (i=first; i<last; i++)
    inew[i-first] = array[i];
  PetscSortInt(last-first, inew);
  for (i=first; i<last; i++)
    array[i] = inew[i];
  free(inew);;
}
