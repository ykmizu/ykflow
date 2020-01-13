
#include "ks_function.h"

void ks_function(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		 Cluster *primal, Is_it *reduced, Vec f, int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section 
  //------------------------------------//-------------------------------------
  Vec nSamplefApprox;                   // fOr reduced case
  Vec function;
  int reducedSize;
  Universe _eqnOfInterest = multiquation->equation;
  Mesh _meshOfInterest;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    _meshOfInterest = reduced->reducedMesh;
    VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &nSamplefApprox); 
  }else if (reduced->reducedSolution == 0){
    _meshOfInterest = primal->state->space;
  }
  reducedSize = _meshOfInterest.node.count*_eqnOfInterest.numStates;
  VecCreateSeq(PETSC_COMM_SELF, reducedSize, &function);
  //---------------------------------------------------------------------------
  // Iplementation
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //Calculates M^-1*R = f
  //---------------------------------------------------------------------------
  //I think this part is slightly wrong. SHould change so work for FULL LSS
  ykflow->Function(ykflow, multiquation, primal, function, reduced, timeNode);
  //---------------------------------------------------------------------------
  // Extract the solution at the nodes refering to hat nodes
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    MatMult(reduced->Zmult, function, nSamplefApprox);
    MatMult(reduced->A, nSamplefApprox, f);
  }else{
    VecCopy(function, f);
  }
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  VecDestroy(&function);
  if (reduced->reducedSolution == 1)
    VecDestroy(&nSamplefApprox);
}


void ks_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *primal, Is_it *reduced, Mat dfdu, int timeNode){
  int elementSize;
  int systemSize;
  Mat dRdU;
  Mat dRdUTemp = NULL;
  Mat dRdUhat;
  Mat dRdUReduced = NULL;
  int columnSize;
  Universe _equation = multiquation->equation;
  elementSize = _equation.numStates*primal->state->basis.nodes;
  columnSize = _equation.numStates*primal->state->space.node.count;
  if (reduced->reducedSolution == 1)//meshSize
    systemSize = _equation.numStates*reduced->reducedMesh.node.count;
  else
    systemSize = _equation.numStates*primal->state->space.node.count;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
   MatCreateSeqBAIJ(PETSC_COMM_SELF, elementSize, systemSize, columnSize,
		   ykflow->numElemCol, NULL, &dRdU);
   if (reduced->reducedSolution == 1){
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, reduced->nBasisFuncs, NULL,
		      &dRdUTemp);
    MatCreateSeqDense(PETSC_COMM_SELF, systemSize, reduced->nBasisFuncs, NULL,
		      &dRdUReduced);
  }
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  ykflow->dFunctiondu(ykflow, multiquation, primal, dRdU, reduced, timeNode);  
  if (reduced->reducedSolution == 1){
    MatMatMultBAIJ(dRdU, reduced->rOBStateBar, dRdUTemp);
    MatMatMult(reduced->Zmult, dRdUTemp,MAT_INITIAL_MATRIX, PETSC_DEFAULT,
    	       &dRdUhat);  //dRdUhat
    MatMatMult(reduced->B, dRdUhat, MAT_REUSE_MATRIX, PETSC_DEFAULT, &dfdu);
  }else{
    MatCopy(dRdU, dfdu, DIFFERENT_NONZERO_PATTERN);
  }
  //---------------------------------------------------------------------------
  // Destroy and free everything
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    MatDestroy(&dRdUTemp);
    MatDestroy(&dRdUhat);
    MatDestroy(&dRdUReduced);
  }
  MatDestroy(&dRdU);
}

void yk_ykflow_function(yk_PrimalSolver *ykflow, Multiverse* multiquation,
			Cluster *primal, Vec vecObj, Is_it *reduced,
			int timeNode){                     
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------

  Mat MijBlock;                                                                
  //---------------------------------------------------------------------------
  // Initialization                                                            
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation->equation, primal->state, timeNode);
  yk_any1D_function(ykflow, multiquation, primal->state, vecObj, reduced);
}                                                                             

void yk_any1D_function(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Galaxy * primal, Vec vecObj, Is_it *reduced){
  ks_spatialStateResidual(multiquation->equation, primal, vecObj, 
                          reduced);                                
  VecScale(vecObj, -1);                     
  MatBlockMult(primal->Mij, vecObj, vecObj);
}


void yk_ykflow_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			   Cluster *primal, Mat matObj, Is_it *reduced,
			   int timeNode){
  Mat dRdUtemp;
  int elementSize;
  int columnSize;
  int systemSize;
  Universe _equation = multiquation->equation;
  //  MatDuplicate(matObj, MAT_DO_NOT_COPY_VALUES, &matCopy);
  elementSize = _equation.numStates*primal->state->basis.nodes;             
  columnSize = _equation.numStates*primal->state->space.node.count;            
  if (reduced->reducedSolution == 1)//meshSize                                 
    systemSize = _equation.numStates*reduced->reducedMesh.node.count;          
  else                                                                         
    systemSize = _equation.numStates*primal->state->space.node.count; 
  MatCreateSeqBAIJ(PETSC_COMM_SELF, elementSize, systemSize, columnSize,       
                   4, NULL, &dRdUtemp); 
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  ks_readSolution(multiquation->equation, primal->state, timeNode);  
  ks_spatialStatedRdU(multiquation->equation, primal->state, matObj, reduced);
  MatAssemblyBegin(matObj, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matObj, MAT_FINAL_ASSEMBLY);   
  MatScale(matObj, -1);
  MatBlockMatMultBAIJ(primal->state->Mij, matObj, dRdUtemp, primal->state, reduced);
  MatCopy(dRdUtemp, matObj, SAME_NONZERO_PATTERN);
  MatDestroy(&dRdUtemp);
}
