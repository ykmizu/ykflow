
#include "tools.h"


void yk_kron(Mat A,  Vec B, double ***kron){
  PetscInt i, j, k; //initialization for iteration

  PetscScalar val[1];
  PetscInt row, col;

  PetscInt vecSize;
  MatGetSize(A, &row, &col);
  VecGetSize(B, &vecSize);
  PetscScalar *feelings = (PetscScalar *) malloc
    (vecSize*sizeof(PetscScalar));
  PetscInt *vecSizeIndex = (PetscInt *) malloc (vecSize*sizeof(PetscInt));
  for (i=0; i<vecSize; i++)
    vecSizeIndex[i] = i;
  int *ptr;
  /* kron = (double *) malloc(row*vecSize*col*sizeof(double)); */
  int r = row*vecSize;
  int c = col;

  *kron = (double **) malloc (r*sizeof(double*));
  (*kron)[0] = (double *) malloc (r*c*sizeof(double));
  for (i = 1; i<r; i++)
    (*kron)[i] = (*kron)[i-1] +c;
  /* int len = sizeof(int *) * r + sizeof(int) * c * r; */
  /* *kron = (int **)malloc(len); */

  /* // ptr is now pointing to the first element in of 2D array */
  /* ptr = (int *)(*kron + r); */

  /* // for loop to point rows pointer to appropriate location in 2D array */
  /* for(i = 0; i < r; i++) */
  /*   (*kron)[i] = (ptr + c * i); */

  VecGetValues(B, vecSize, vecSizeIndex, feelings);
  for (i=0; i<row; i++){
    for (j=0; j<col; j++){

      MatGetValues(A, 1, &i, 1, &j, val);
      for (k=0; k< vecSize; k++)
	(*kron)[i*vecSize+k][j]  = val[0]*feelings[k];
    }
  }
  free(feelings);
  free(vecSizeIndex);

}

/* void yk_kron(Mat A, Vec B, Mat *kron){ */
/*   int i, j; */
/*   PetscInt row, col, vecSize; */
/*   //Get the sizes of the Matrices and Vector */
/*   MatGetSize(A, &row, &col); */
/*   VecGetSize(B, &vecSize); */
/*   //Create MatBlock here of rthe A matrix */
/*   MatCreate(PETSC_COMM_SELF, kron);//create snapshot matrix */
/*   MatSetSizes(*kron, vecSize, 1, row*vecSize, col); */
/*   MatSetBlockSizes(*kron, vecSize,1); */

/*   MatSetType(*kron, MATSEQBAIJ);           //set snapshot matrix type */
/*   MatSetFromOptions(*kron); */
/*   MatSeqAIJSetPreallocation(*kron,col, NULL); */


/*   PetscScalar *ptr; */
/*   /\* kron = (double *) malloc(row*vecSize*col*sizeof(double)); *\/ */
/*   int len = sizeof(PetscScalar *) * row + sizeof(PetscScalar) * col * row; */
/*   PetscScalar **Avalues = (PetscScalar **)malloc(len); */

/*   // ptr is now pointing to the first element in of 2D array */
/*   ptr = (PetscScalar *)(Avalues + row); */

/*   // for loop to point rows pointer to appropriate location in 2D array */
/*   for(i = 0; i < row; i++) */
/*     Avalues[i] = (ptr + col * i); */

/*   printf("ok okokokoko\n"); */
/*   getchar(); */
/*   //iterate through all the values int he A Matrix */
/*   PetscInt *rowIndex; */
/*   PetscInt *colIndex; */
/*   for (i=0; i<row; i++) */
/*     rowIndex[i] = i; */
/*   for (i=0; i<col; i++) */
/*     colIndex[i] = i; */
/*   MatGetValues(A, row, &rowIndex, col, &colIndex, *Avalues); */
/*   printf("lets seee\n"); */
/*   getchar(); */

/*   /\* for (i=0; i<row; i++){ *\/ */
/*   /\*   for (j=0; j<col; j++){ *\/ */

/*   /\*   } *\/ */
/*   /\* } *\/ */

/* } */

double dot(double* a, double*b, int size){
  int i; //initialization for iteration
  double sum = 0;
  for (i=0; i<size; i++)
    sum += a[i]*b[i];
  return sum;
}

void projection(Universe eqn, Vec f, Vec vp){
  double pnt, numerator, denominator;      //projection values
  VecDot(f, vp, &numerator);     //dot product the numerator
  VecDot(f, f, &denominator);    //dot product the denominator
  pnt = numerator/denominator; //Now calculate the projection of input vector
  VecAXPY(vp, -pnt, f);
  /* for (i=0; i<eqn.numStates*U->space.size; i++) */
  /*   vp[i] -= - pnt*f[i]; */
}

void getSegArray(Universe equation, Galaxy *U, double *R, int seg, Vec vp){
  int i;                         //initialization for iteration
  int systemSize = equation.numStates*U->space.node.count;
  double vpValue;
  for (i=0; i<systemSize; i++){
    vpValue = R[seg*systemSize+i];
    VecSetValue(vp, i, vpValue, INSERT_VALUES);
  }
}



void array2Vec(Universe eqn, Galaxy *U, Vec vcon){
  int i, j, k;
  PetscInt location;
  for (i=0; i<U->space.elem.count; i++){
    for (j=0; j<U->basis.nodes; j++){
      for (k=0; k<eqn.numStates; k++){
	location = (i*eqn.numStates*U->basis.nodes)+(j*eqn.numStates)+k;
	VecSetValues(vcon, 1, &location,
		     &U->solution[k].array[i*U->basis.nodes+j],
		     INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(vcon);
  VecAssemblyEnd(vcon);
}

void vec2Array(Universe eqn, Galaxy *U, Vec vcon){
  int i;
  int systemSize = U->space.node.count*eqn.numStates;
  PetscInt *index = (PetscInt *) malloc (systemSize*sizeof(PetscInt));
  PetscScalar *values = (PetscScalar *)malloc(systemSize*sizeof(PetscScalar));
  for (i=0; i<systemSize; i++)
    index[i] = i;
  VecGetValues(vcon, systemSize, index, values); //extract value to S
  data2Array(eqn, U, values);
  free(index);
  free(values);
}



void data2Array(Universe eqn, Galaxy *U, double *u){
  int i, j, k;
  int location;
  for (i=0; i<U->space.elem.count; i++){
    for (j=0; j<U->basis.nodes; j++){
      for (k=0; k<eqn.numStates; k++){
	location = (i*eqn.numStates*U->basis.nodes)+(j*eqn.numStates)+k;
        //location = (i*eqn.numStates+j)*U->basis.nodes+k;
        U->solution[k].array[i*U->basis.nodes+j] = u[location];
      }
    }
  }
}

void array2Data(Universe eqn, Galaxy *U, double *u){
  int i, j, k;
  PetscInt location;
  for (i=0; i<U->space.elem.count; i++){
    for (j=0; j<U->basis.nodes; j++){
      for (k=0; k<eqn.numStates; k++){
    //    for (j=0; j<eqn.numStates; j++){
    // for (k=0; k<U->basis.nodes; k++){
	location = (i*eqn.numStates*U->basis.nodes)+(j*eqn.numStates)+k;
	u[location] = U->solution[k].array[i*U->basis.nodes+j];
      }
    }
  }

  /* for (i=0; i<U->space.elem.count; i++){ */
  /*   for (j=0; j<eqn.numStates; j++){ */
  /*     for (k=0; k<U->basis.p+1; k++){ */
  /*       location = (i*eqn.numStates+j)*(U->basis.p+1)+k; */
  /*       u[location] = U->solution[j].array[i*(U->basis.p+1)+k];  */
  /*     } */
  /*   } */
  /* } */
}


void ks_Mat2MatBlock(Mat A, PetscInt m, PetscInt idxm[], PetscInt n,
		     PetscInt idxn[], Mat ABlock, InsertMode iora){
  int i;
  PetscScalar **Aij = (PetscScalar **) malloc (m*sizeof(PetscScalar *));
  Aij[0] = (PetscScalar *) malloc (m*n*sizeof(PetscScalar));
  for (i=0; i<m; i++)
    Aij[i] = (*Aij + n*i);

  MatGetValues(A, m, idxm, n, idxn, *Aij);
  MatSetValues(ABlock, m, idxm, n, idxn, *Aij, iora);
  MatAssemblyBegin(ABlock, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(ABlock, MAT_FINAL_ASSEMBLY);
  free(Aij[0]);
  free(Aij);
}

void ks_MatRow2Vec(Mat M, PetscInt m, PetscInt idxm[], Vec V, PetscInt rowi,
                   InsertMode iora){
  PetscScalar V_i[m];
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  MatGetValues(M, 1, &rowi, m, idxm, V_i);
  VecSetValues(V, m, idxm, V_i, iora);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
}

void ks_MatCol2Vec(Mat M, PetscInt m, PetscInt idxm[], Vec V, PetscInt coli,
                   InsertMode iora){
  PetscScalar V_i[m];
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  MatGetValues(M, m, idxm, 1, &coli, V_i);
  VecSetValues(V, m, idxm, V_i, iora);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
}

void ks_Vec2MatCol(Mat M, PetscInt m, PetscInt idxm[], PetscInt col, Vec V,
                   InsertMode iora){
  PetscScalar V_i[m];
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  VecGetValues(V, m, idxm, V_i); //extract value to S
  MatSetValues(M, m, idxm,  1, &col, V_i, iora);
}

void ks_Vec2MatRow(Mat M, PetscInt m, PetscInt idxm[], PetscInt row, Vec V,
                   InsertMode iora){
  PetscScalar V_i[m];

  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  VecGetValues(V, m, idxm, V_i);
  MatSetValues(M, 1, &row, m, idxm, V_i, iora);
}

int *createMeshMap(Mesh *original){
  int i;
  int totalNumElements = original->elem.count;
  int originalElement_i;
  int lastValue = original->elem.array[totalNumElements-1];
  int *getMeshIndex = (int *) malloc ((lastValue+1)*sizeof(int));
  for (i=0; i<totalNumElements; i++){
    originalElement_i = original->elem.array[i];
    getMeshIndex[originalElement_i] = i;
  }
  return getMeshIndex;
}


void solveEquations(Mat A, Vec b, Vec x){
  KSP ksp;
  PC pc;
  PetscInt iterNum;
  VecCopy(b, x);
  KSPCreate(PETSC_COMM_SELF, &ksp); //Create the KSP object
  //--------------------------------------------------------------------------
  // Perform Newton Solve to find the Won)
  //--------------------------------------------------------------------------
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, KSPGMRES);
  KSPSetFromOptions(ksp);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCBJACOBI);
  KSPSolve(ksp, x, x);
  KSPGetIterationNumber(ksp, &iterNum);
}


int minInt(int A, int B){
  if (A > B)
    return B;
  return A;
}

/* void findCognateState(Multiverse multiquation, Cluster *object, int node, */
/*                       Is_it *reduced){ */
/*   Universe _eqnOfInterest; */
/*   if (reduced->reducedSolution == 1) */
/*     _eqnOfInterest = multiquation.equationReduced; */
/*   else */
/*     _eqnOfInterest = multiquation.equation; */
/*   //--------------------------------------------------------------------------- */
/*   // For the reduced case, primalApprox!= primal               */
/*   //--------------------------------------------------------------------------- */
/*   if (reduced->reducedSolution==1){ //if we are solving HROM-LSS equations    */
/*     ks_readSolution(multiquation.equation, object->stateApprox, node); */
/*   }else{ */
/*     ks_readSolution(_eqnOfInterest, object->stateApprox, node); */
/*     ks_copySolutions(multiquation.equation, object->stateApprox, object->state); */
/*   } */
/* } */
//NEED TO REWRITE THIS SHIT. SIGH. FML. THIS SUCKS.
void MatMatMultBAIJ(Mat A, Mat B, Mat C){
  int i;
  PetscInt rowi = 0;
  PetscInt coli;
  PetscInt rowiA = 0;
  PetscInt coliA;
  Vec Bi;
  Vec Ci;
  PetscInt *rowiIndex = NULL;
  PetscInt *rowiAIndex = NULL;
  MatGetSize(B, &rowi, &coli);
  MatGetSize(A, &rowiA, &coliA);
  rowiIndex =  (PetscInt *) malloc (rowi*sizeof(PetscInt));
  rowiAIndex =  (PetscInt *) malloc (rowiA*sizeof(PetscInt));
  for (i=0; i<rowi; i++)
    rowiIndex[i] = i;
  for (i=0; i<rowiA; i++)
    rowiAIndex[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, rowi, &Bi);
  VecCreateSeq(PETSC_COMM_SELF, rowiA, &Ci);
  for (i=0; i<coli; i++){
    ks_MatCol2Vec(B, rowi, rowiIndex, Bi, i, INSERT_VALUES);
    MatMult(A, Bi, Ci);
    ks_Vec2MatCol(C, rowiA, rowiAIndex, i, Ci, INSERT_VALUES);
  }
  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  free(rowiIndex);
  free(rowiAIndex);
  VecDestroy(&Bi);
  VecDestroy(&Ci);
}

void inverseMatLU(Mat A, Mat *invA){
  int i;
  IS row, col;
  PetscInt rowA, colA;
  Mat identity;
  MatFactorInfo info;
  MatGetSize(A, &rowA, &colA);

  MatCreateSeqDense(PETSC_COMM_SELF, rowA, colA, NULL, &identity);
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, colA, NULL, invA);
  for (i=0; i<colA; i++)
    MatSetValue(identity, i, i, 1, INSERT_VALUES);
  MatAssemblyBegin(identity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(identity, MAT_FINAL_ASSEMBLY);
  MatGetOrdering(A, MATORDERINGNATURAL, &row, &col);
  MatLUFactor(A, row, col, &info);
  MatMatSolve(A, identity , *invA);
  MatSetUnfactored(*invA);

  MatDestroy(&identity);
  ISDestroy(&row);
  ISDestroy(&col);
}

void MatBlockMult(Mat A, Vec x, Vec b){
  int i, j;
  int *Bindex;
  int numBlocks;
  PetscInt row, col;
  PetscInt vecSize;
  PetscInt *index;
  PetscInt *bindex;
  PetscInt *rowIndex;
  double *tempApprox;
  double *AtempApprox;
  Vec xBlock;
  Vec bApproxVec;
  //----------------------------------------------------------------------------
  // Initialization
  //----------------------------------------------------------------------------
  MatGetSize(A, &row, &col); //For now row and col should equal
  VecGetSize(x, &vecSize);
  numBlocks = vecSize/col;
  index = (int *) malloc (col*sizeof(int));
  rowIndex = (int *) malloc (col*sizeof(int));
  Bindex = (int *) malloc (row*sizeof(int));
  bindex = (int *) malloc (row*sizeof(int));
  tempApprox = (double *) malloc (col*sizeof(double));
  AtempApprox = (double *) malloc (row*sizeof(double));
  for (i=0; i<col; i++)
    rowIndex[i] = i;
  for (i=0; i<row; i++)
    bindex[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, col, &xBlock);
  VecCreateSeq(PETSC_COMM_SELF, row, &bApproxVec);
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  for (i=0; i<numBlocks; i++){ //number of elements
    for (j=0; j<col; j++)
      index[j] = i*col+j;
    for (j=0; j<row; j++)
      Bindex[j] = i*row+j;
    VecGetValues(x, col, index, tempApprox);
    VecSetValues(xBlock, col, rowIndex, tempApprox, INSERT_VALUES);
    VecAssemblyBegin(xBlock);
    VecAssemblyEnd(xBlock);
    MatMult(A, xBlock, bApproxVec);
    VecGetValues(bApproxVec, row, bindex, AtempApprox);
    VecSetValues(b, row, Bindex, AtempApprox, INSERT_VALUES);
  }
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  //----------------------------------------------------------------------------
  // Destroy everything that you see here
  //----------------------------------------------------------------------------
  VecDestroy(&xBlock);
  VecDestroy(&bApproxVec);
  free(index);
  free(bindex);
  free(rowIndex);
  free(Bindex);
  free(tempApprox);
  free(AtempApprox);
}


void MatBlockMatMultBAIJ(Mat A, Mat B, Mat C, Galaxy *primal, Is_it *reduced){
  int i, j, k;
  PetscScalar **MatArray;
  PetscInt rowA, colA, rowB, colB;
  PetscInt *rowAIndex;
  PetscInt *colBIndex;
  PetscInt *index;
  Mat BBlock, CBlock;
  int numBlocks;
  int i_epsilon;
  int *getMeshIndex = createMeshMap(&primal->space);
  Mesh _meshOfInterest;
  if (reduced->hrom==1)
    _meshOfInterest = reduced->reducedMesh;
  else
    _meshOfInterest = primal->space;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatGetSize(A, &rowA, &colA); //Gives the size of A  = primal->basis.p+1
  MatGetSize(B, &rowB, &colB);
  numBlocks = rowB/rowA;
  index = (PetscInt *) malloc (rowA*sizeof(PetscInt));
  rowAIndex = (PetscInt *) malloc (rowA*sizeof(PetscInt));
  colBIndex = (PetscInt *) malloc (rowB*sizeof(PetscInt));
  MatArray = (PetscScalar **) malloc (rowA*sizeof(PetscScalar *));
  MatArray[0] = (PetscScalar *) malloc (rowA*rowA*sizeof(PetscScalar));
  for (i=0; i<rowA; i++){
    rowAIndex[i] = i;
    MatArray[i] = (*MatArray + rowA*i);
  }
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, rowA, NULL, &BBlock);
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, rowA, NULL, &CBlock);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<rowA; j++)//For the rows the index changes
      index[j] = i*rowA+j;
    for (j=-1; j<2; j++){
      if ((i_epsilon == 0 && j == -1)||
	(i_epsilon == primal->space.elem.count-1 && j == 1))
      	continue;
      for (k=0; k<rowA; k++)
	colBIndex[k] = (i_epsilon+j)*rowA+k;
      MatGetValues(B, rowA, index, rowA, colBIndex, *MatArray);
      MatSetValues(BBlock, rowA, rowAIndex, rowA, rowAIndex, *MatArray,
		   INSERT_VALUES);
      MatAssemblyBegin(BBlock, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(BBlock, MAT_FINAL_ASSEMBLY);

      MatMatMult(A, BBlock, MAT_REUSE_MATRIX, PETSC_DEFAULT, &CBlock);
      MatGetValues(CBlock, rowA, rowAIndex, rowA, rowAIndex, *MatArray);
      MatSetValues(C, rowA, index, rowA, colBIndex, *MatArray,
		   INSERT_VALUES);
    }
  }
  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
  MatDestroy(&BBlock);
  MatDestroy(&CBlock);
  free(MatArray[0]);
  free(MatArray);
  free(rowAIndex);
  free(colBIndex);
  free(index);
  free(getMeshIndex);
}
void MatBlockMatMultDense(Mat A, Mat B, Mat C){
  //The point of this function is multiple matrix B by A where A is a smaller
  //block like mass matrix.
  int i, j;
  int numBlocks = 0;
  PetscInt *index;
  PetscInt rowA, colA, rowB, colB;
  PetscInt *rowAIndex;
  PetscInt *colBIndex;
  PetscScalar **MatArray;
  Mat BBlock, CBlock;
  clock_t start, end;
  double cpu_time_used;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  start = clock();
  MatGetSize(B, &rowB, &colB);
  MatGetSize(A, &rowA, &colA);
  numBlocks = rowB/rowA;
  index = (PetscInt *) malloc (rowA*sizeof(PetscInt));
  rowAIndex = (PetscInt *) malloc (rowA*sizeof(PetscInt));
  colBIndex = (PetscInt *) malloc (colB*sizeof(PetscInt));
  for (i=0; i<rowA; i++)
    rowAIndex[i] = i;
  for (i=0; i<colB; i++)
    colBIndex[i] = i;
  MatArray = (PetscScalar **) malloc (rowA*sizeof(PetscScalar *));
  MatArray[0] = (PetscScalar *) malloc (rowA*colB*sizeof(PetscScalar));
  for (i=0; i<rowA; i++)
    MatArray[i] = (*MatArray + colB*i);
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, colB, NULL,  &BBlock);
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, colB, NULL, &CBlock);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  printf("CALL %d\n", numBlocks);
  getchar();
  for (i=0; i<numBlocks; i++){
    for (j=0; j<rowA; j++)
      index[j] = i*rowA+j;
    MatGetValues(B, rowA, index, colB, colBIndex, *MatArray);
    MatSetValues(BBlock, rowA, rowAIndex, colB, colBIndex, *MatArray,
		 INSERT_VALUES);
    MatAssemblyBegin(BBlock, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BBlock, MAT_FINAL_ASSEMBLY);
    MatMatMult(A, BBlock, MAT_REUSE_MATRIX, PETSC_DEFAULT, &CBlock);
    MatGetValues(CBlock, rowA, rowAIndex, colB, colBIndex, *MatArray);
    MatSetValues(C, rowA, index, colB, colBIndex, *MatArray, INSERT_VALUES);
  }
  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  // Destroy and free everything
  //---------------------------------------------------------------------------
  MatDestroy(&BBlock);
  MatDestroy(&CBlock);
  free(index);
  free(rowAIndex);
  free(colBIndex);
  free(MatArray[0]);
  free(MatArray);

}





/* void MatMatMult(Mat a, Mat b, Mat y){ */
/*   PetscInt rowa, cola; */
/*   PetscInt rowb, colb; */
/*   int* rowaIndex; */
/*   int *rowbIndex; */
/*   Vec bVec; */
/*   Vec yVec; */
/*   //---------------------------------------------------------------------------- */
/*   // Initialization */
/*   //---------------------------------------------------------------------------- */
/*   VecCreateSeq(PETSC_COMM_SELF, colb, &bVec); */
/*   VecCreateSeq(PETSC_COMM_SELF, cola, &yVec); */
/*   MatCreateSeqDense(PETSC_COMM_SELF, rowa, colb, NULL, &y); */
/*   MatGetSize(b, &rowb, &colb); */
/*   MatGetSIze(a, &rowa, &cola); */
/*   for (i=0; i<rowa; i++) */
/*     rowaIndex[i] = i;                   */
/*   for (i=0; i<rowb; i++) */
/*     rowbIndex[i] = i; */
/*   //---------------------------------------------------------------------------- */
/*   // Implementation */
/*   //---------------------------------------------------------------------------- */
/*   for (i=0; i<colb; i++){ */
/*     ks_MatCol2Vec(b, rowb, rowbIndex, bVec, i, INSERT_VALUES); */
/*     MatMult(a, bVec, yVec); */
/*     ks_Vec2MatCol(y, rowa, rowaIndex, yVec, i, ISNERT_VALUES); */
/*   } */
/* } */
