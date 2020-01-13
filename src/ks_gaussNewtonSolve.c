#include "ks_gaussNewtonSolve.h"
//Zr = B
//Zj = A
//Is U_rom the same its?  Vr is optional if not in use trype NULL
void ks_gaussNewtonSolve(Multiverse multiquation, Cluster *primalApprox,
			 Is_it *reduced, int tSolve, int *iterationCount){
  int i;                           //initialization for iteration
  int flag = 0;
  int systemSize;
  int meshSize;
  int *index;
  int gaussCount = 0;
  double *snapshotsRJ;
  double alfa;
  PetscScalar normp;
  Vec p;                              //iono what p is. Look this part up
  Vec residual=NULL, residualTemp=NULL, residualRed=NULL;  //Contains residual
  Vec primalReducedVec, primalApproxVec, snapJacobian;
  Mat dRdU, dRdUTemp=NULL, dRdUTempHat=NULL, dRdUbasis = NULL;     //Jacobian
  FILE *residualSnapshotFile;
  FILE *jacobianSnapshotFile;
  Mesh _meshOfInterest;
  Mat rOBState = NULL;
  Universe equation = multiquation.equation;
  Universe equationReduced = multiquation.equationReduced;
  //---------------------------------------------------------------------------
  // Mallocing things occur here/Initizalize things here
  //---------------------------------------------------------------------------
  systemSize = equation.numStates*primalApprox->self->space.node.count;
  index = (int *) malloc (systemSize*sizeof(int));
  snapshotsRJ = (double *) malloc (systemSize*sizeof(double)); 
  for (i=0; i<systemSize; i++)
    index[i] = i;
  if (reduced->hrom == 1){
    rOBState = reduced->rOBStateBar;                                        
    _meshOfInterest = reduced->reducedMesh;
    VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &residualRed);
    VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &residual);
    MatCreateSeqDense(PETSC_COMM_SELF, equation.numStates*
		      _meshOfInterest.node.count, reduced->nBasisFuncs, NULL,
		      &dRdUbasis);
    MatCreateSeqDense(PETSC_COMM_SELF, reduced->nSampleNodes,
		      reduced->nBasisFuncs, NULL, &dRdUTempHat);
    MatCreateSeqDense(PETSC_COMM_SELF, reduced->nSampleNodes,                  
                      reduced->nBasisFuncs, NULL, &dRdU); 
  }else if (reduced->hrom == 0){
    rOBState = reduced->rOBState;                                              
    _meshOfInterest = primalApprox->self->space; 
    VecCreateSeq(PETSC_COMM_SELF, systemSize, &residual);
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, reduced->nBasisFuncs, NULL,
		      &dRdU);  
  }
  meshSize = equation.numStates*_meshOfInterest.node.count;
  VecCreateSeq(PETSC_COMM_SELF, meshSize, &residualTemp); 
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &snapJacobian);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &p); //Always
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &primalReducedVec);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &primalApproxVec);
  MatCreateSeqBAIJ(PETSC_COMM_SELF, primalApprox->self->basis.p+1, meshSize,
		   systemSize, 3, NULL, &dRdUTemp);
  //---------------------------------------------------------------------------
  // Set the solution to the previous solution i-1
  //---------------------------------------------------------------------------
  array2Vec(equationReduced, primalApprox->reduced, primalReducedVec);
  do{
    MatMult(rOBState, primalReducedVec, primalApproxVec); //update 
    vec2Array(equation, primalApprox->self, primalApproxVec);
    //-------------------------------------------------------------------------
    // Calculate the Residual and Jacobian
    //-------------------------------------------------------------------------
    ks_totalResidual(multiquation, primalApprox, residualTemp, reduced,tSolve);
    ks_totaldRdU(multiquation, primalApprox, dRdUTemp, reduced, tSolve);
    //-------------------------------------------------------------------------
    // Gauss Newton with Approximated Tensors for hyper-reduced order modeling
    //-------------------------------------------------------------------------
    if (reduced->hrom == 1){
      MatMult(reduced->Zmult, residualTemp, residualRed);  //Calcualte 
      MatMult(reduced->B, residualRed, residual);
      MatMatMultBAIJ(dRdUTemp, rOBState, dRdUbasis);
      MatMatMult(reduced->Zmult, dRdUbasis, MAT_REUSE_MATRIX, PETSC_DEFAULT,
		 &dRdUTempHat);
      MatMatMult(reduced->A, dRdUTempHat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&dRdU);
    }else{
      VecCopy(residualTemp, residual);
      MatMatMultBAIJ(dRdUTemp, rOBState, dRdU);
    }
    if (reduced->hrom == 0){
      //-----------------------------------------------------------------------
      // Residual Snapshots
      //-----------------------------------------------------------------------
      VecGetValues(residual, systemSize, index, snapshotsRJ);
      residualSnapshotFile = fopen("residualSnapshot.dat", "a");
      fprintf(residualSnapshotFile, "%0.16f", snapshotsRJ[0]);
      for (i=1; i<systemSize; i++)  //Number Sample Rows
        fprintf(residualSnapshotFile, " %0.16f", snapshotsRJ[i]);
      fprintf(residualSnapshotFile, "\n");
      fclose(residualSnapshotFile);  
    }
    //-------------------------------------------------------------------------
    // Least Squares solve for the direction p
    //-------------------------------------------------------------------------
    VecScale(residual, -1);
    linearLeastSquares(dRdU, residual, p);
    (*iterationCount)++;
    //-----------------------------------------------------------------------
      // Jacobian Snapshots                                                   
      //-----------------------------------------------------------           
    if (reduced->hrom == 0){
      MatMult(dRdU, p, snapJacobian);
      VecGetValues(snapJacobian, systemSize, index, snapshotsRJ);
      jacobianSnapshotFile = fopen("jacobianSnapshot.dat", "a");
      fprintf(jacobianSnapshotFile, "%0.16f", snapshotsRJ[0]);
      for (i=1; i<systemSize; i++)
        fprintf(jacobianSnapshotFile, " %0.16f", snapshotsRJ[i]);
      fprintf(jacobianSnapshotFile, "\n");
      fclose(jacobianSnapshotFile);
    }
    //-------------------------------------------------------------------------
    // Perform line search or set \alpha to 1
    //-------------------------------------------------------------------------
    alfa = 1;
    VecNorm(p, 2, &normp);
    printf("Gauss-Newton's Iteration %d, Error %0.16e\n", gaussCount, normp);
    gaussCount++;
    if (normp <pow(10, -8))
      flag = 1;
    //-------------------------------------------------------------------------
    // Update Delu_k
    //-------------------------------------------------------------------------
    VecAXPY(primalReducedVec, alfa, p);
    vec2Array(equationReduced, primalApprox->reduced, primalReducedVec);     
  }while (flag == 0);  //Convergencecheck here
  //Convert the final solution the reduced solution at current t into Utype for
  array2Vec(equationReduced, primalApprox->reduced, primalReducedVec);
  MatMult(rOBState, primalReducedVec, primalApproxVec);
  vec2Array(equation, primalApprox->self, primalApproxVec);
  //---------------------------------------------------------------------------
  // Destroy all the malloc stuff and initialization stuff OMG
  //---------------------------------------------------------------------------
  VecDestroy(&p);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalApproxVec);
  VecDestroy(&residual);
  VecDestroy(&residualTemp);
  VecDestroy(&snapJacobian);
  MatDestroy(&dRdU);
  MatDestroy(&dRdUTemp);
  if (reduced->hrom == 1){
    VecDestroy(&residualRed);
    MatDestroy(&dRdUTemp);
    MatDestroy(&dRdUTempHat);
    MatDestroy(&dRdUbasis);
  }
  free(snapshotsRJ);
  free(index);
}
