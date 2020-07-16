#include "ks_podv2.h"

void createSnapshotState(Universe equation, Galaxy *primal, Mat *snapshot,
			 int dss){
  int i;
  int snapshotCount = floor(primal->time.count/dss)+1; //number of snapshots
  int snapshot_i;
  int systemSize = primal->space.node.count*equation.numStates;

  int* index = (int *) malloc (systemSize *sizeof(int));  //array of index vale
  PetscScalar* states =(PetscScalar *) malloc (systemSize*sizeof(PetscScalar));
  Vec primal_i;                        //state: values at node i
  //---------------------------------------------------------------------------
  // I n i t i a l i z a t i o n
  //---------------------------------------------------------------------------
  for (i=0; i<systemSize; i++)         //fill in the index values
    index[i] = i;
  MatCreateSeqDense(PETSC_COMM_SELF, systemSize, snapshotCount, NULL, snapshot);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primal_i);
  //---------------------------------------------------------------------------
  // I m p l e n t a t i o n
  //---------------------------------------------------------------------------
  for (i=0; i<snapshotCount; i++){     //loop through all the num of snapshots
    if (i==snapshotCount-1)
      snapshot_i = primal->time.count;
    else if (i==0)
      snapshot_i = 1;
    else
      snapshot_i = dss*i;  //first snapshot file to put into a S matrix/ski
    ks_readSolution(equation, primal, snapshot_i);     //extrac
    array2Data(equation, primal, states); //to vector form
    VecSetValues(primal_i, systemSize, index, states, INSERT_VALUES);  //Set U_
    ks_Vec2MatCol(*snapshot, systemSize, index, i, primal_i, INSERT_VALUES);
  }
  MatAssemblyBegin(*snapshot, MAT_FINAL_ASSEMBLY);  //No need to change snapsho
  MatAssemblyEnd(*snapshot, MAT_FINAL_ASSEMBLY);    //No need to change snapshot
  //---------------------------------------------------------------------------
  // T e r m i n a t i o n
  //---------------------------------------------------------------------------
  free(states);
  free(index);
  VecDestroy(&primal_i);
}

/* void approxState(Universe equation, Galaxy *primal, int dss, int snapshot_i, */
/*                  Mat reducedOrderBasis){ */
/*   int i;                               //initialization\1;95;0c of iterations */
/*   int lineCount = 0;                   //line counter for Ur reading */
/*   int systemSize = primal->space.node.count*equation.numStates; */
/*   int snapshotCount = floor(primal->time.count/dss); //Number of state files t */
/*   char outputFile[1000]; */
/*   char buffer[1024]; */
/*   char* current_line; */
/*   char* lines; */
/*   FILE* solutionFile; */
/*   Galaxy* primalReduced = (Galaxy *) malloc (sizeof(Galaxy)); */
/*   PetscInt index[systemSize]; */
/*   PetscScalar primalk[systemSize]; */
/*   Vec primalReduced_i;                 //Saves the reduced solution to a Vec */
/*   Vec primal_i; */
/*   Vec primal_0; */
/*   int snapshotCountIndex[snapshotCount]; */
/*   //--------------------------------------------------------------------------- */
/*   // Malloc and initialize everything before implementation         */
/*   //--------------------------------------------------------------------------- */
/*   for (i=0; i<snapshotCount; i++) */
/*     snapshotCountIndex[i] = i; */
/*   for (i=0; i<systemSize; i++) */
/*     index[i] = i; */
/*   primalReduced->solution = (Array *) malloc (sizeof(Array));    */
/*   initArray(primalReduced->solution, snapshotCount); */
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &primal_0); */
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &primal_i); */
/*   VecCreateSeq(PETSC_COMM_SELF, snapshotCount, &primalReduced_i); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   //--------------------------------------------------------------------------- */
/*   // Read in the reduced states for these equations and save to Vecs */
/*   //--------------------------------------------------------------------------- */
/*   sprintf(outputFile, "%s_%s_%d_%d_reduced_%d.dat", equation.nameEqn, */
/*           primal->utypeId, primal->basis.p, primal->time_method, */
/*           snapshot_i); */
/*   solutionFile = fopen(outputFile, "r"); */
/*   if (solutionFile == NULL){ */
/*     perror("Error"); */
/*     exit(EXIT_FAILURE); */
/*   }else{ */
/*     while ((current_line = fgets(buffer, sizeof(buffer), solutionFile))){ */
/*       if (current_line[0]!= '\n'){ */
/*         lines =strtok(current_line, " "); */
/*         primalReduced->solution[0].array[lineCount] = atof(lines); */
/*       } */
/*       lineCount ++; */
/*     } */
/*   } */
/*   fclose(solutionFile); */
/*   VecSetValues(primalReduced_i, snapshotCount, snapshotCountIndex, */
/*                primalReduced->solution[0].array, INSERT_VALUES); */
/*   VecAssemblyBegin(primalReduced_i); */
/*   VecAssemblyEnd(primalReduced_i); */
/*   //---------------------------------------------------------------------------- */
/*   // Calculate the final state approximation */
/*   //---------------------------------------------------------------------------- */
/*   ks_readSolution(equation, primal, 0);//extract the initia */
/*   array2Vec(equation, primal, primal_0);  //Put into a Vec    */
/*   MatMultAdd(reducedOrderBasis, primalReduced_i, primal_0, primal_i); */
/*   VecGetValues(primal_i, systemSize, index, primalk); */
/*   data2Array(equation, primal, primalk); */
/*   primal->time.node = snapshot_i; */
/*   //---------------------------------------------------------------------------- */
/*   // Free and Destroy everything */
/*   //---------------------------------------------------------------------------- */
/*   VecDestroy(&primal_i); */
/*   VecDestroy(&primalReduced_i); */
/*   delArray(primalReduced->solution);   */
/*   free(primalReduced); */
/* } */



void moorePenrosePseudoInvSparse(Mat A, PetscInt rowSize, PetscInt colSize, Mat *Aplus){
  int i, j, k;
  PetscInt *rowIndex = (PetscInt *) malloc (rowSize*sizeof(PetscInt));
  PetscInt *colIndex = (PetscInt *) malloc (colSize*sizeof(PetscInt));
  PetscScalar *singular_i = (PetscScalar *)malloc(rowSize*sizeof(PetscScalar));
  Vec singular;
  Vec decomposition;
  Vec decomptemp;
  SVD svd;
  PetscReal sigma;
  PetscInt wtf;
  PetscScalar *V_i = (PetscScalar *) malloc (colSize*sizeof(PetscScalar));

  //----------------------------------------------------------------------------
  //  I n i t i a l i z a t i o n
  //---------------------------------------------------------------------------
  for (i=0; i<rowSize; i++)
    rowIndex[i] = i;
  for (i=0; i<colSize; i++)
    colIndex[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, rowSize, &singular);
  VecCreateSeq(PETSC_COMM_SELF, colSize, &decomposition);
  VecCreateSeq(PETSC_COMM_SELF, colSize, &decomptemp);
  //---------------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Perform SVD
  //---------------------------------------------------------------------------
  SVDCreate(PETSC_COMM_SELF, &svd);
  SVDSetOperator(svd, A);
  SVDSetFromOptions(svd);
  SVDSetDimensions(svd, colSize,
                   PETSC_DEFAULT, PETSC_DEFAULT);
  SVDSolve(svd);

  SVDGetConverged(svd, &wtf);
  for (i=0; i<wtf; i++){ //Iterate through the number of converged iteratio
    SVDGetSingularTriplet(svd, i, &sigma, singular, decomposition);
    if (sigma >0)
      sigma = 1.0/sigma;
    else
      sigma = sigma;
    VecScale(decomposition, sigma);
    VecGetValues(singular, rowSize, rowIndex, singular_i);
    for (j=0; j<rowSize; j++){
      VecCopy(decomposition, decomptemp);
      VecScale(decomptemp, singular_i[j]);
      //---------------------------------------------------------------------------
      // Implementation
      //---------------------------------------------------------------------------
      VecGetValues(decomptemp, colSize, colIndex, V_i); //extract value to S
      for (k=0; k<colSize;k++){
        if (fabs(V_i[k])>pow(10,-8)){
	  MatSetValue(*Aplus,k,j,V_i[k],ADD_VALUES);
	  }
      }
      /* ks_Vec2MatCol(*Aplus, colSize, colIndex, j, decomptemp, ADD_VALUES); */
    }
  }

  MatAssemblyBegin(*Aplus, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Aplus, MAT_FINAL_ASSEMBLY);
  //----------------------------------------------------------------------------
  // T e r  m i n a t i o n
  //----------------------------------------------------------------------------
  VecDestroy(&singular);
  VecDestroy(&decomposition);
  VecDestroy(&decomptemp);
  SVDDestroy(&svd);
  free(rowIndex);
  free(colIndex);
  free(singular_i);
  free(V_i);
}




void moorePenrosePseudoInv(Mat A, PetscInt rowSize, PetscInt colSize, Mat *Aplus){
  int i, j;
  PetscInt *rowIndex = (PetscInt *) malloc (rowSize*sizeof(PetscInt));
  PetscInt *colIndex = (PetscInt *) malloc (colSize*sizeof(PetscInt));
  PetscScalar *singular_i = (PetscScalar *)malloc(rowSize*sizeof(PetscScalar));
  Vec singular;
  Vec decomposition;
  Vec decomptemp;
  SVD svd;
  PetscReal sigma;
  PetscInt wtf;
  //----------------------------------------------------------------------------
  //  I n i t i a l i z a t i o n
  //---------------------------------------------------------------------------
  for (i=0; i<rowSize; i++)
    rowIndex[i] = i;
  for (i=0; i<colSize; i++)
    colIndex[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, rowSize, &singular);
  VecCreateSeq(PETSC_COMM_SELF, colSize, &decomposition);
  VecCreateSeq(PETSC_COMM_SELF, colSize, &decomptemp);
  printf("BPOOOP\n");
  yk_MatCreateSeqDense(Aplus, colSize, rowSize);
  MatZeroEntries(*Aplus);
  printf("MOTHEIREURHER\n");
  //---------------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Perform SVD
  //---------------------------------------------------------------------------
  SVDCreate(PETSC_COMM_SELF, &svd);
  SVDSetOperator(svd, A);
  SVDSetFromOptions(svd);
  SVDSetDimensions(svd, colSize,
                   PETSC_DEFAULT, PETSC_DEFAULT);
  SVDSolve(svd);
  /* SVDConvergedReason reason; */
  /* SVDGetConvergedReason(svd, &reason); */
  /* printf("%d\n", reason); */
  /* printf("REASON\n"); */
  /* getchar(); */

  SVDGetConverged(svd, &wtf);
  for (i=0; i<wtf; i++){ //Iterate through the number of converged iteratio
    SVDGetSingularTriplet(svd, i, &sigma, singular, decomposition);
    /* if (sigma = 0){ */
    /*   fprintf(stderr, "Division by zero! Exiting...\n"); */
    /* } */
    /* printf("sigma %g\n", sigma); */
    /* getchar(); */
    if (sigma >0)
      sigma = 1.0/sigma;
    else
      sigma = sigma;
    VecScale(decomposition, sigma);
    VecGetValues(singular, rowSize, rowIndex, singular_i);
    for (j=0; j<rowSize; j++){
      VecCopy(decomposition, decomptemp);
      VecScale(decomptemp, singular_i[j]);
      ks_Vec2MatCol(*Aplus, colSize, colIndex, j, decomptemp, ADD_VALUES);
    }
  }
  MatAssemblyBegin(*Aplus, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Aplus, MAT_FINAL_ASSEMBLY);
  //----------------------------------------------------------------------------
  // T e r  m i n a t i o n
  //----------------------------------------------------------------------------
  VecDestroy(&singular);
  VecDestroy(&decomposition);
  VecDestroy(&decomptemp);
  SVDDestroy(&svd);
  free(rowIndex);
  free(colIndex);
  free(singular_i);

}

void properOrthogonalDecompose(Mat snapshot, int systemSize,
                               PetscInt *numSingularValues, Mat *A){
  int i;                               //initialization for iteration
  int index[systemSize];               //indices for row size
  PetscReal sigma;                     //eigenvalues
  PetscScalar A_i[systemSize];         //column values for singular matrix
  Vec singular;                        //singular vectors rotation
  SVD singularValueDecomposition;      //create svd object
  PetscInt row, col;

  MatGetSize(snapshot, &row, &col);
  //----------------------------------------------------------------------------
  // I n i t i a l i z a t i o n
  //----------------------------------------------------------------------------
  for (i=0; i<systemSize; i++)
    index[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &singular);
  //----------------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Create the Singular Value solver and set various options
  //----------------------------------------------------------------------------
  SVDCreate(PETSC_COMM_SELF, &singularValueDecomposition);
  SVDSetOperator(singularValueDecomposition, snapshot);
  SVDSetFromOptions(singularValueDecomposition);
  SVDSetDimensions(singularValueDecomposition, *numSingularValues,
                   PETSC_DEFAULT, PETSC_DEFAULT);
  //----------------------------------------------------------------------------
  // Solver the Singular Value System
  //----------------------------------------------------------------------------
  //SVDSetTolerances(singularValueDecomposition, tol, maxits);
  SVDSolve(singularValueDecomposition);
  //SVDGetConverged(singularValueDecomposition, numSingularValues);
  //----------------------------------------------------------------------------
  // Truncate the first its<=Ns left singular vectors to obtain the ROB matrix
  //----------------------------------------------------------------------------
  SVDGetDimensions(singularValueDecomposition, numSingularValues, NULL, NULL);
  MatCreateSeqDense(PETSC_COMM_SELF, systemSize, *numSingularValues, NULL, A);
  for (i=0; i<*numSingularValues; i++){ //Iterate through the numerged iteratio
    SVDGetSingularTriplet(singularValueDecomposition, i,
                          &sigma, singular, NULL);    //extract rotation
    VecGetValues(singular, systemSize, index, A_i); //Retrieve u_i
    MatSetValues(*A, systemSize, index, 1, &i, A_i,INSERT_VALUES); //Set to ROB
  }
  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); //No need to change the ROB Mat
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);   //No need to change the ROB Mat
  //----------------------------------------------------------------------------
  // T e r m i n a t i o n
  //----------------------------------------------------------------------------
  VecDestroy(&singular);
  SVDDestroy(&singularValueDecomposition);
}


/* void reduceMatRows(Universe equation, Galaxy *primal, Mat A, int column, */
/*                    Mat *Ahat, Mesh *altered){ */

/*   int i, j;                            //initialization for iteration */
/*   int systemSize = altered->node.count*equation.numStates; //size of altered sta */
/*   PetscScalar *A_i = (PetscScalar *) malloc (column*sizeof(PetscScalar)); */
/*   //holds temporary altered matrix A */
/*   int *columnIndex = (int *) malloc (column*sizeof(int)); */
/*   //index for number of columns in A */
/*   int fullNodeCount_i;                 //node value respect to the original mesh */
/*   int alteredNodeCount_i;              //node value respect to the altered mesh */
/*   //---------------------------------------------------------------------------- */
/*   // I n i t i a l i za t i o n */
/*   //---------------------------------------------------------------------------- */
/*   for (i=0; i<column; i++) */
/*     columnIndex[i] = i; */
/*   MatCreateSeqDense(PETSC_COMM_SELF, systemSize, column, NULL, Ahat); //ROB */
/*   //---------------------------------------------------------------------------- */
/*   // I m p l e m e n t a t i o n */
/*   //---------------------------------------------------------------------------- */
/*   for (i=0; i<altered->elem.count; i++){ */
/*     for (j=0; j<primal->basis.p+1; j++){ */
/*       fullNodeCount_i = altered->elem.array[i]*(primal->basis.p+1)+j; */

/*       alteredNodeCount_i = i*(primal->basis.p+1)+j; */
/*       MatGetValues(A, 1, &fullNodeCount_i, column, columnIndex, A_i); */

/*       MatSetValues(*Ahat, 1, &alteredNodeCount_i, column, columnIndex, */
/*                    A_i, INSERT_VALUES); */
/*     } */
/*   } */

/*   MatAssemblyBegin(*Ahat, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(*Ahat, MAT_FINAL_ASSEMBLY); */
/*   free(columnIndex); */
/*   free(A_i); */
/* } */


/* void ks_createROBjacobian(Universe eqn, Utype *U, */
/*                           Mat V, int rowV, int Ns, intArray* relem, Mat S){ */
/*   int i, j;                         //initialization of iteration */
/*   int Ns_i;                           //time node for the state solution  */
/*   int numfiles = floor(U->time.count/Ns);     //Number of state files to be read */
/*   int sysSize = U->space.size*eqn.numStates;  */
/*   int sysSizeR = relem->size*(U->basis.p+1)*eqn.numStates; //Number of rows in r */
/*   PetscInt indexR[sysSizeR];                  //Indicies for the rows of S */
/*   PetscScalar S_i[sysSizeR];                  //Contains the columns of S at i */
/*   Vec residual;                               //Save Residual values */
/*   Vec U_0; */
/*   Vec V_i; */
/*   Mat dRdU; */
/*   int s_index; */
/*   PetscScalar V_iarray[sysSizeR]; */
/*   Vec S_iv; */
/*   //---------------------------------------------------------------------------- */
/*   // Malloc and initialize everything before implementation */
/*   //---------------------------------------------------------------------------- */
/*   for (i=0; i<sysSizeR; i++) */
/*     indexR[i] = i; */
/*   VecCreateSeq(PETSC_COMM_SELF, sysSizeR, &residual); */
/*   VecCreateSeq(PETSC_COMM_SELF, sysSize, &U_0); */
/*   VecCreateSeq(PETSC_COMM_SELF, sysSize, &V_i); */
/*   MatCreateSeqDense(PETSC_COMM_SELF, sysSizeR, numfiles, NULL, &S); */
/*   MatCreateSeqBAIJ(PETSC_COMM_SELF, U->basis.p+1, sysSize, sysSize, 3, NULL, */
/*                    &dRdU);   */
/*   //---------------------------------------------------------------------------- */
/*   // Implemenation */
/*   //---------------------------------------------------------------------------- */
/*   //---------------------------------------------------------------------------- */
/*   // Read in the initial conditions */
/*   //---------------------------------------------------------------------------- */
/*   ks_readSolution(eqn, U, 0);             //Read in the initial conditions state */
/*   array2Vec(eqn, U, U_0);                 //Put into a Vec */
/*   for (i=0; i<numfiles; i++){         //Loop through all the number of snapshots */
/*     //-------------------------------------------------------------------------- */
/*     // Approximate state solution at Ns*i+1; */
/*     //------------------------------------------------------------------------- */
/*     U->time.node = Ns*i+1;         //First snapshot file to put into a S matrix */
/*     approxState(eqn, U, U_0, V, rowV); */
/*     //-------------------------------------------------------------------------- */
/*     // Now calculate the residual snapshot matrix given the approximate solution */
/*     //-------------------------------------------------------------------------- */
/*     ks_totaldRdU(eqn, U, dRdU, 1, relem); */
/*     for (j=0; j<Ns; j++){ */
/*       s_index = i*numfiles+j; */
/*       MatGetValues(V, sysSize, indexR, 1, &j, V_iarray); */
/*       VecSetValues(V_i, sysSize, indexR, V_iarray, INSERT_VALUES);  */
/*       MatMult(dRdU, V_i, S_iv); */
/*       MatSetValues(S, sysSizeR, indexR, 1, &s_index, S_i, INSERT_VALUES);  */
/*     } */
/*   } */
/*   MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY); */
/*   VecDestroy(&residual); */
/*   VecDestroy(&U_0); */

/* } */
