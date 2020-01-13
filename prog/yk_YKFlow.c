

#include "yk_sovlerCFD.h"

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void yk_findAdjacentElems1D(yk_PrimalSolver *solverObj, int elem, int *elemSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                
  //------------------------------------//-------------------------------------
  //1D
  elemSet[elem]= 1;
  elemSet[elem-1] = 1;
  elemSet[elem+1] = 1;
}

void yk_findBoundarySeeds1D(yk_PrimalSolver *primalSolverObj, int nSampNodes,
			    int *numSeeds, int *nodeSet){
  
  *numSeeds = 2;
  nodeSet[0] = 0;
  nodeSet[primalSolverObj->multiquation->equation.numStates] =
    primalSolverObj->systemSize-1;;
}

yk_ykflow_totalResidual(yk_PrimalSolver *primalSolverObj, Cluster *object,
			Vec residual, Mat Jacobian, Is_it *reduced,
			int timeSolveNum){
  ks_totalResidual(primalSolverObj, object, residual, reduced, timeSolveNum);
  ks_totaldRdU(primalSolverObj, object, Jacobian, reduced, timeSolveNum);
}

void yk_ykflow_function(yk_PrimalSolver *ykflow, Cluster *objectG, Vec vecObj,
			Is_it *reduced, int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                  
  //------------------------------------//-------------------------------------
  Mat MijBlock;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  ks_massBlock(multiquation.equation, primal->self, &MijBlock);
  MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
  inverseMatLU(MijBlock, &primal->self->Mij);
  ks_spatialStateResidual(multiquation.equation, primal->state, residuals,
			  reduced);
  VecScale(residuals, -1); 
  MatBlockMult(primal->state->Mij, residuals, vecObj);
}

void yk_ykflow_dfunctiondu(yk_PrimalSolver *ykflow, Clust *objectG,
			   Mat matObj, Is_it *reduced, int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  ks_spatialStatedRdU(multiquation.equation, primal->state, dRdU, reduced);
  MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY);
  MatScale(dRdU, -1);
  MatBlockMatMultDense(primal->state->Mij, dRdU, matObj);
}

void ks_jbaraverage(yk_PrimalSolver
