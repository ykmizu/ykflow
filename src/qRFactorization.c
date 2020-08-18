
#include "qRFactorization.h"

//||Ax-b||

void linearLeastSquaresHROM(Mat A, Vec b, Vec x){
  Mat Q;
  Mat R;
  Vec rightHS;
  KSP ksp;
  PC pc;
  PetscInt maxits = 1000;
  PetscReal rtol = pow(10,-12);
  PetscReal abstol = pow(10,-50);
  PetscReal dtol= 1000;
  PetscInt rowA;
  PetscInt colA;
  PetscScalar RHS;
  Mat C;
  Vec y;
  Vec cc;
  PetscReal val;
  //----------------------------------------------------------------------------
  // Initialization
  //----------------------------------------------------------------------------
  MatGetSize(A, &rowA, &colA);
  MatCreateSeqDense(PETSC_COMM_SELF, rowA, colA, NULL, &Q);
  MatCreateSeqDense(PETSC_COMM_SELF, colA, colA, NULL, &R);
  KSPCreate(PETSC_COMM_SELF, &ksp);
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  MatView(A, PETSC_VIEWER_STDOUT_SELF);
  VecView(b, PETSC_VIEWER_STDOUT_SELF);


  qRFactorization(A, Q, R);

  VecCreateSeq(PETSC_COMM_SELF, colA, &rightHS);
  /* getchar(); */
  /* MatTransposeMatMult(A,A ,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C); */
  MatMultTranspose(Q, b, rightHS);
  MatView(Q, PETSC_VIEWER_STDOUT_SELF);
  KSPSetOperators(ksp, R, R);
  KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);
  /* KSPSetType(ksp, KSPGMRES); */
  KSPSetFromOptions(ksp);
  /* KSPGetPC(ksp, &pc); */
  /* PCSetType(pc, PCBJACOBI); */
  KSPSolve(ksp, rightHS, x);
  /* VecCopy(b, x); //x is the one that minimizations Ax-b */
  /* VecDestroy(&rightHS); */
  KSPDestroy(&ksp);
}

void linearLeastSquares(Mat A, Vec b, Vec x){
  Mat Q;
  Mat R;
  Vec rightHS;
  KSP ksp;
  PC pc;
  PetscInt maxits = 1000;
  PetscReal rtol = pow(10,-12);
  PetscReal abstol = pow(10,-50);
  PetscReal dtol= 1000;
  PetscInt rowA;
  PetscInt colA;
  PetscScalar RHS;
  Mat C;
  Vec y;
  Vec cc;
  PetscReal val;
  //----------------------------------------------------------------------------
  // Initialization
  //----------------------------------------------------------------------------
  MatGetSize(A, &rowA, &colA);
  VecCreateSeq(PETSC_COMM_SELF, colA, &rightHS);
  /* MatCreateSeqDense(PETSC_COMM_SELF, rowA, colA, NULL, &Q); */
  /* MatCreateSeqDense(PETSC_COMM_SELF, colA, colA, NULL, &R); */
  KSPCreate(PETSC_COMM_SELF, &ksp);
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  MatTransposeMatMult(A,A ,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
  MatMultTranspose(A, b, rightHS);
  /* VecNorm(rightHS, NORM_2, &val); */
  /* printf("lmult J,R = %0.16e\n", val); */
  /* printf("next best thing\n"); */
  /* qRFactorization(A, Q, R); //QR decomposition */
  /* printf("sketch\n"); */
  /* MatMultTranspose(Q, b, rightHS); //Calculation of the right hand side */
  /* VecNorm(rightHS, 2, &RHS); */
  /* printf("sneeze\n"); */
  //  Solve for x that minimizes the value of ||Ax-b||2
  /* KSPSetOperators(ksp, R, R); */
  /* KSPSetTolerances(ksp, rtol, abstol, dtol, maxits); */
  /* KSPSetType(ksp, KSPGMRES); */
  /* KSPSetFromOptions(ksp); */
  /* KSPGetPC(ksp, &pc); */
  /* PCSetType(pc, PCBJACOBI); */
  /* printf("okokokok\n"); */
  /* KSPSolve(ksp, rightHS, rightHS); */
  /* printf("too slow"); */
  /* VecCopy(rightHS, x); //x is the one that minimizations Ax-b */
  /* VecDestroy(&rightHS); */
  /* KSPDestroy(&ksp); */
  /* MatDestroy(&Q); */
  /* MatDestroy(&R); */


  KSPSetOperators(ksp, C, C);
  KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);
  /* KSPSetType(ksp, KSPGMRES); */
  KSPSetFromOptions(ksp);
  /* KSPGetPC(ksp, &pc); */
  /* PCSetType(pc, PCBJACOBI); */
  KSPSolve(ksp, rightHS, x);
  /* VecCopy(b, x); //x is the one that minimizations Ax-b */
  /* VecDestroy(&rightHS); */
  KSPDestroy(&ksp);
  MatDestroy(&C);
  VecDestroy(&rightHS);
}


void qRFactorization(Mat A, Mat Q, Mat R){
  BV bvA;
  Mat tempQ;
  //  BVCreate(PETSC_COMM_SELF, &bvA);
  BVCreateFromMat(A, &bvA);
  BVSetType(bvA, BVMAT);
  BVOrthogonalize(bvA, R);
  BVCreateMat(bvA, &tempQ);
  MatCopy(tempQ, Q, SAME_NONZERO_PATTERN);
  MatDestroy(&tempQ);
  BVDestroy(&bvA);

}


/* void qRFactorization(Mat A, Mat Q, Mat R){ */
/*   int i, j;                //initialization for iteration */
/*   int *numRowsIndex; */
/*   int numRows; */
/*   int numCols; */
/*   PetscScalar r_ij; */
/*   PetscScalar r_ii; */
/*   PetscScalar norm_w; */
/*   Vec q_i; */
/*   Vec qrA_i; */
/*   Vec qrA_j; */
/*   Mat qrA; */
/*   //---------------------------------------------------------------------------- */
/*   // Initialization */
/*   //---------------------------------------------------------------------------- */
/*   MatGetSize(A, &numRows, &numCols); */
/*   numRowsIndex = (int *) malloc (numRows*sizeof(int)); */
/*   VecCreateSeq(PETSC_COMM_SELF, numRows, &q_i); */
/*   VecCreateSeq(PETSC_COMM_SELF, numRows, &qrA_i); */
/*   VecCreateSeq(PETSC_COMM_SELF, numRows, &qrA_j); */
/*   MatCreateSeqDense(PETSC_COMM_SELF, numRows, numCols, NULL, &qrA); */
/*   for (i=0; i<numRows; i++) */
/*     numRowsIndex[i] = i; */
/*   //---------------------------------------------------------------------------- */
/*   // Implementation */
/*   //---------------------------------------------------------------------------- */
/*   MatCopy(A, qrA, SAME_NONZERO_PATTERN); */
/*   for (i=0; i<numCols; i++){ */
/*     ks_MatCol2Vec(qrA, numRows, numRowsIndex, qrA_i, i, INSERT_VALUES); */
/*     VecNorm(qrA_i, NORM_2, &r_ii); */
/*     MatSetValue(R, i, i, r_ii, INSERT_VALUES);    //Set r11 to norm */
/*     norm_w = 1/r_ii; */
/*     VecScale(qrA_i, norm_w); */
/*     ks_Vec2MatCol(Q, numRows, numRowsIndex, i, qrA_i, INSERT_VALUES); */
/*     MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY); */
/*     MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY); */
/*     for (j = i+1; j<numCols; j++){  //mSnapindex needs to edited. ug */
/*       ks_MatCol2Vec(Q, numRows, numRowsIndex, q_i, i, INSERT_VALUES); */
/*       ks_MatCol2Vec(qrA, numRows, numRowsIndex, qrA_j, j,INSERT_VALUES); */
/*       VecDot(q_i, qrA_j, &r_ij); */
/*       MatSetValue(R, i, j, r_ij, INSERT_VALUES);    //Set r11 to norm */
/*       VecAXPY(qrA_j, -r_ij, q_i); */
/*       ks_Vec2MatCol(qrA, numRows, numRowsIndex, j, qrA_j,INSERT_VALUES); */
/*       MatAssemblyBegin(qrA, MAT_FINAL_ASSEMBLY); */
/*       MatAssemblyEnd(qrA, MAT_FINAL_ASSEMBLY); */
/*     } */
/*   } */
/*   MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY); */
/*   free(numRowsIndex); */
/*   VecDestroy(&q_i); */
/*   VecDestroy(&qrA_j); */
/*   VecDestroy(&qrA_i); */
/*   MatDestroy(&qrA); */
/* } */
