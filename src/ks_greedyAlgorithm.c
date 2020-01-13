
#include "ks_greedyAlgorithm.h"

/* void yk_greedyAlgorithmST_space(Cluster *primal, int seedCount, */
/* 				Is_it *reduced){ */
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
