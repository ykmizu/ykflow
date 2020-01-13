#include "ks_offlineMatrix.h"

void ks_offlineMatrix(Universe equation, Galaxy *object, Galaxy *objectApprox,
		       Is_it *reduced){
  int systemSize = equation.numStates*object->space.node.count;
  Mat Ba, Qc, C;
  Mat rOBResidualHatInv = NULL;
  Mat rOBJacobianHatInv = NULL;
  Mat rOBResidualHat = NULL;
  Mat rOBJacobianHat = NULL;
  //---------------------------------------------------------------------------
  // Second Round of Mallocing Due to New Information From Reduced Mesh
  //---------------------------------------------------------------------------
  MatCreateSeqDense(PETSC_COMM_SELF, reduced->nSampleNodes,                   
                    reduced->nSampleNodes, NULL, &reduced->A); //Rc           
  MatCreateSeqDense(PETSC_COMM_SELF, reduced->nSampleNodes,                   
                    reduced->nSampleNodes, NULL, &reduced->B);  
  MatCreateSeqDense(PETSC_COMM_SELF, systemSize, reduced->nSampleNodes, NULL,
		    &Qc);
  //---------------------------------------------------------------------------
  // Restricted Matrices: Build these using indices from nodeSet
  //---------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Build the restricted indentity matrices                                   
  //-------------------------------------------------------------------------- 
  //-------------------------------------------------------------------------- 
  // Creation of Restricted Matrices                                           
  //-------------------------------------------------------------------------- 
  MatMatMult(reduced->Z, reduced->rOBResidual, MAT_INITIAL_MATRIX,             
             PETSC_DEFAULT, &rOBResidualHat);                                
  MatMatMult(reduced->Z, reduced->rOBJacobian, MAT_INITIAL_MATRIX,             
             PETSC_DEFAULT, &rOBJacobianHat);
  //---------------------------------------------------------------------------
  // Offline computation of online matrices/Algorithim Rc, B online matrices
  //---------------------------------------------------------------------------
  moorePenrosePseudoInv(rOBResidualHat, reduced->nSampleNodes,
			reduced->nBasisFuncsRJ, &rOBResidualHatInv);
  moorePenrosePseudoInv(rOBJacobianHat, reduced->nSampleNodes,
			reduced->nBasisFuncsRJ, &rOBJacobianHatInv);

  MatMatMult(reduced->rOBJacobian, rOBJacobianHatInv, MAT_INITIAL_MATRIX,
	     PETSC_DEFAULT, &C);
  qRFactorization(C, Qc, reduced->A);
  MatTransposeMatMult(Qc, reduced->rOBResidual, MAT_INITIAL_MATRIX,
		      PETSC_DEFAULT, &Ba);
  MatMatMult(Ba, rOBResidualHatInv, MAT_REUSE_MATRIX, PETSC_DEFAULT,
	     &reduced->B);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&rOBResidualHat);
  MatDestroy(&rOBJacobianHat);
  MatDestroy(&rOBResidualHatInv);
  MatDestroy(&rOBJacobianHatInv);
  MatDestroy(&C);
  MatDestroy(&Qc);
  MatDestroy(&Ba);
}

void ks_offlineMatrixLSS(Universe equation, Galaxy *object, Is_it *reduced){
  PetscInt systemFullSize;
  PetscInt columnSize;
  
  MatGetSize(reduced->rOBJacobian, &systemFullSize, &columnSize);
  Mat C, D;
  Mat rOBResidualHat, rOBJacobianHat;
  Mat rOBResidualHatInv, rOBJacobianHatInv;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  MatCreateSeqDense(PETSC_COMM_SELF, reduced->nBasisFuncs,
                    reduced->nSampleNodes, NULL, &reduced->A); //Rc
  MatCreateSeqDense(PETSC_COMM_SELF, reduced->nBasisFuncs,
                    reduced->nSampleNodes, NULL, &reduced->B);
  MatCreateSeqDense(PETSC_COMM_SELF, systemFullSize, reduced->nSampleNodes,
		    NULL, &C);
  MatCreateSeqDense(PETSC_COMM_SELF, systemFullSize, reduced->nSampleNodes,
		    NULL, &D);
  //---------------------------------------------------------------------------
  // Restricted Matrices: Build these using indices from nodeSet
  //---------------------------------------------------------------------------
  MatMatMult(reduced->Z, reduced->rOBResidual, MAT_INITIAL_MATRIX,
             PETSC_DEFAULT, &rOBResidualHat);
  MatMatMult(reduced->Z, reduced->rOBJacobian, MAT_INITIAL_MATRIX,
             PETSC_DEFAULT, &rOBJacobianHat);
  moorePenrosePseudoInv(rOBResidualHat, reduced->nSampleNodes,
			reduced->nBasisFuncsRJ, &rOBResidualHatInv);
  moorePenrosePseudoInv(rOBJacobianHat, reduced->nSampleNodes,
			reduced->nBasisFuncsRJ, &rOBJacobianHatInv);
  //---------------------------------------------------------------------------
  // Residual offline Matrix for LSS
  //---------------------------------------------------------------------------
  MatMatMult(reduced->rOBResidual, rOBResidualHatInv, MAT_REUSE_MATRIX,
	     PETSC_DEFAULT, &C);
  MatTransposeMatMult(reduced->rOBState, C, MAT_REUSE_MATRIX, PETSC_DEFAULT,
		      &reduced->A);
  //  MatScale(reduced->A, -1);
  //---------------------------------------------------------------------------
  // Jacobian offline Matrix for LSS
  //---------------------------------------------------------------------------
  MatMatMult(reduced->rOBJacobian, rOBJacobianHatInv, MAT_REUSE_MATRIX,
	     PETSC_DEFAULT, &D);
  MatTransposeMatMult(reduced->rOBState, D, MAT_REUSE_MATRIX, PETSC_DEFAULT,
		      &reduced->B);
  //MatScale(reduced->B, -1);  //From the Residual having to go to the other sid
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&rOBResidualHat);
  MatDestroy(&rOBJacobianHat);
  MatDestroy(&rOBResidualHatInv);
  MatDestroy(&rOBJacobianHatInv);
  MatDestroy(&C);
  MatDestroy(&D);
  
}

void restrictedIdentity(Universe equation,Galaxy *object, Galaxy *objectApprox,
			PetscInt *nodeSet, Is_it *reduced){
  int i, j;
  intArray reducedN = reduced->reducedMesh.node;
  intArray primalN = object->space.node;
  intArray approxN = objectApprox->space.node;
  int nStates = equation.numStates;
  int reducedVecSize = nStates*reducedN.count;
  int primalVecSize = nStates*primalN.count;
  int approxVecSize = nStates*approxN.count;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatCreateSeqAIJ(PETSC_COMM_SELF,reduced->nSampleNodes, primalVecSize, 1,NULL,
		  &reduced->Z); 
  MatCreateSeqAIJ(PETSC_COMM_SELF, reducedVecSize, primalVecSize, 1, NULL,
		  &reduced->Zu);
  MatCreateSeqAIJ(PETSC_COMM_SELF, approxVecSize, primalVecSize, 1, NULL,
		  &reduced->Zp);
  //---------------------------------------------------------------------------
  //Z matrix with sample nodes
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->nSampleNodes; i++)
    MatSetValue(reduced->Z, i, nodeSet[i], 1, INSERT_VALUES);
  MatAssemblyBegin(reduced->Z, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(reduced->Z, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  //Z matrix for the approximate mesh
  //---------------------------------------------------------------------------
  for (i=0; i<approxN.count; i++)
    for (j=0; j<nStates; j++)
      MatSetValue(reduced->Zp, i*nStates+j, approxN.array[i]*nStates+j,
		  1,INSERT_VALUES);
  MatAssemblyBegin(reduced->Zp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(reduced->Zp, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  //Z matriz for the reduced mesh
  //---------------------------------------------------------------------------
  for (i=0; i<reducedN.count; i++)
    for (j=0; j<nStates; j++)
      MatSetValue(reduced->Zu, i*nStates+j, reducedN.array[i]*nStates+j,
		  1, INSERT_VALUES);
  MatAssemblyBegin(reduced->Zu, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(reduced->Zu, MAT_FINAL_ASSEMBLY);
  //---------------------------------------------------------------------------
  //Zmult here for finding hat when the size is of the reducedMesh size 
  //---------------------------------------------------------------------------
  MatMatTransposeMult(reduced->Z, reduced->Zu, MAT_INITIAL_MATRIX,
		      PETSC_DEFAULT, &reduced->Zmult);
}

void restrictedMatrices(Universe equation, Galaxy *object,
			Galaxy *objectApprox, Mat *rOBResidualHat,
			Mat *rOBJacobianHat, PetscInt *nodeSet,
			Is_it *reduced){
  //---------------------------------------------------------------------------
  // Build the restricted indentity matrices
  //---------------------------------------------------------------------------
  restrictedIdentity(equation, object, objectApprox, nodeSet, reduced);
  //---------------------------------------------------------------------------
  // Creation of Restricted Matrices
  //---------------------------------------------------------------------------
  MatMatMult(reduced->Zp, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
	     &reduced->rOBStateBar);
  MatMatMult(reduced->Zu, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
	     &reduced->rOBStateRed);
  MatMatMult(reduced->Z, reduced->rOBResidual, MAT_INITIAL_MATRIX,
	     PETSC_DEFAULT, rOBResidualHat);
  MatMatMult(reduced->Z, reduced->rOBJacobian, MAT_INITIAL_MATRIX,
	     PETSC_DEFAULT, rOBJacobianHat);
}
