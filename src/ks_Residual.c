


#include "ks_Residual.h"
//-------------------------------------------------------------------------
// Calculates the Residual/ Add functions for different equations
// Each solution U will have a name and this will determine which   
// R(n+1)
// Yukiko Shimizu
// April 23, 2016
// Error Estimation for Chaotic Systems
//-------------------------------------------------------------------------
void ks_totalResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster* primal, Vec residual, Is_it *reduced,
		      int timeSolveNum){
  int i, j, k ,n;                      //Iteration for initialization pu
  int in;
  int temp = primal->self->space.dx;
  int systemSize;
  int reducedSize;
  double tempSum;                    //Temporary sum for residual calc
  double *c = getBDFCoef(timeSolveNum);//Get coefficient;
  Vec residualTemp, residualOne;       //Temporary Vecs for multiplication
  Mat Mij;
  int *getMeshIndex;
  Universe _eqnOfInterest;
  Cluster ** primalCalc =
    (Cluster **) malloc ((timeSolveNum+1)*sizeof(Cluster*));
  //---------------------------------------------------------------------------
  // I n i t i a l i z a t i o n
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)                                         
    _eqnOfInterest = multiquation->equationReduced;          
  else if (reduced->reducedSolution == 0 &&                                    
           strncmp(primal->self->utypeId, "psi", 3)==0)    
    _eqnOfInterest = multiquation->equationLSS;       
  else                                                                         
    _eqnOfInterest = multiquation->equation;     
  if (reduced->hrom == 1 && reduced->reducedSolution == 0) //gnat + rom       
    reducedSize = reduced->reducedMesh.node.count*_eqnOfInterest.numStates;
  else if (reduced->hrom == 1 && reduced->reducedSolution == 1) //gnat + lss
    reducedSize = primal->self->space.node.count*_eqnOfInterest.numStates; 
  else if (reduced->hrom == 0) //no gnat
    reducedSize = primal->self->space.node.count*_eqnOfInterest.numStates; 
  systemSize = primal->self->space.node.count*_eqnOfInterest.numStates;
  getMeshIndex = createMeshMap(&primal->self->space);
  for (i=0; i<timeSolveNum+1; i++){ //Create an array of struct pointers
    primalCalc[i] = (Cluster *) malloc (sizeof(Cluster));
    primalCalc[i]->self = (Galaxy *) malloc (sizeof(Galaxy));
  }
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &residualTemp);
  VecCreateSeq(PETSC_COMM_SELF, reducedSize, &residualOne);  
  MatCreateSeqAIJ(PETSC_COMM_SELF, reducedSize, systemSize,
		  primal->self->basis.p+1, PETSC_NULL, &Mij);
  //---------------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //---------------------------------------------------------------------------
  for (i=0; i<timeSolveNum+1; i++){
    ks_copyUtype(primal->self, primalCalc[i]->self);
    copySystem(_eqnOfInterest, primal->self, primalCalc[i]->self);
    primalCalc[i]->self->beta = primal->self->beta;
    primalCalc[i]->self->djdu = primal->self->djdu;
    strcpy(primalCalc[i]->clusterId, primal->clusterId);
  }
  primalCalc[0]->self->time.node = primal->self->time.node; //This is t
  if (strncmp(primal->self->utypeId, "psi1", 4) == 0)
    for (i=1; i<timeSolveNum+1; i++)
      primalCalc[i]->self->time.node = primal->self->time.node + i;
  else
    for (i=1; i<timeSolveNum+1; i++)
      primalCalc[i]->self->time.node = primal->self->time.node-i;
  //-----------------------------------------------------------------------
  // Make the initial guess for the Un+1 solution, Extract solution
  //-----------------------------------------------------------------------
  ks_copySolutions(_eqnOfInterest, primal->self, primalCalc[0]->self);
  for (i=1; i<timeSolveNum+1; i++){
    ks_readSolution(_eqnOfInterest, primalCalc[i]->self,
		    primalCalc[i]->self->time.node);  
  }
  //-----------------------------------------------------------------------
  // Mass Matrix
  //-----------------------------------------------------------------------
  ks_mass(_eqnOfInterest, primal->self, Mij, reduced); //Don't need the   
  //-----------------------------------------------------------------------
  // Find the Spatial Residual
  //-----------------------------------------------------------------------    
  ks_spatialResidual(ykflow, multiquation, primal, residual, reduced);
  MatScale(Mij, 1/primal->self->time.dt);
  //-----------------------------------------------------------------------
  // Find the Total Residual
  //-----------------------------------------------------------------------
  for (i=0; i<primal->self->space.elem.count; i++){ //number of elements
    for (j=0; j<_eqnOfInterest.numStates; j++){
      for (k=0; k<primal->self->basis.p+1; k++){
        in = (i*_eqnOfInterest.numStates+j)*(primal->self->basis.p+1)+k;
	tempSum = 0; 
	for (n=0; n<timeSolveNum+1; n++)
          tempSum += c[n]*primalCalc[n]->self->solution[j].array
	    [i*(primal->self->basis.p+1)+k];
	VecSetValue(residualTemp, in, tempSum, ADD_VALUES); 
      }
    }
  }
  VecAssemblyBegin(residualTemp);
  MatMult(Mij, residualTemp, residualOne);
  VecAXPY(residual, 1, residualOne);
  //-----------------------------------------------------------------------
  // Destroy Everything
  //-----------------------------------------------------------------------
  MatDestroy(&Mij);
  for (i=0; i<timeSolveNum+1; i++){
    destroySystem(_eqnOfInterest, primalCalc[i]->self);
    free(primalCalc[i]->self);
    free(primalCalc[i]);
  }
  free(primalCalc);
  free(c);
  VecDestroy(&residualTemp);
  VecDestroy(&residualOne);
  free(getMeshIndex); 
}


void ks_spatialResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			Cluster *primal, Vec residual, Is_it *reduced){
  Universe _equation = multiquation->equation;
  VecSet(residual, 0);
  if (strcmp(primal->self->utypeId, "state")==0)
    ks_spatialStateResidual(_equation, primal->self, residual, reduced);
  else if(strncmp(primal->self->utypeId, "psi2", 4)==0)
    ks_res_tangent(ykflow, multiquation, primal, residual, reduced); 
  else if(strncmp(primal->self->utypeId, "psi1",4)==0)
    ks_res_lagrange(ykflow, multiquation, primal, residual, reduced);
  VecAssemblyBegin(residual);                                               
  VecAssemblyEnd(residual);     

  //-----------------------------------------------------------------------
  // Copy Residuals into a Vec for Petsc use and then free
  //-----------------------------------------------------------------------
}

void ks_spatialStateResidual(Universe equation, Galaxy *object, Vec residual,
			    Is_it *reduced){
  VecSet(residual, 0);
  if (strncmp(equation.nameEqn, "meks", 4)==0){ 
    ks_res_advection(equation, object, residual, reduced);         
    ks_res_burgers(equation, object, residual, reduced);
    ks_res_diffusion(equation, object, residual, reduced);
    ks_res_fourth_IPDG(equation, object, residual, reduced);
  }
}


void yk_ykflow_totalResidual(yk_PrimalSolver *primalSolverObj,
			     Multiverse *multiquation, Cluster *object, 
			     Vec residual, Mat Jacobian, Is_it *reduced, 
			     int timeSolveNum){                       
  ks_totalResidual(primalSolverObj, multiquation, object, residual, reduced,
		   timeSolveNum);   
  if (Jacobian !=NULL)
    ks_totaldRdU(primalSolverObj, multiquation, object, Jacobian, reduced,
		 timeSolveNum);       
}  
