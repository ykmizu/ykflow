

#include "yk_calculateError.h"


void yk_estimateError(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *primal_h, Cluster *primal_H, Cluster *primal,
		      Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i, j, k, m;                             //initialization for iteration      
  int numfiles = primal->self->time.count/reduced->nTimeSegs;
  int timeSolve = 1;                    //designates which BDF solver to use   
  int left, right;                                                             
  int systemSize;                                                              
  int systemReducedSize;                                                       
  int hrom, IsItreduced;                                                       
  char adjointName[1000];                                                      
  PetscReal rtol = pow(10,-16);                                             
  PetscReal abstol = pow(10,-16);                               
  PetscReal dtol= 1000;    
  double timeFinal = primal->self->time.globalT_f; 
  double timeInitial = primal->self->time.globalT_0;
  double timelength = (timeFinal-timeInitial)/reduced->nTimeSegs;              
  double delJ = 0;                                                             
  PetscInt *index;                                                             
  double del = 0;
  PetscInt maxits = 1000; 
  double timedel = 0;
  Vec residual;                                                                
  Vec nSamplefApprox;
  Vec psi, truePsi, adjoint;                                                   
  Universe _eqnOfInterest;              //Current universe                     
  Galaxy *_solOfInterest;                                                      
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));                        
  Psi1->self = (Galaxy *) malloc (sizeof(Galaxy));                             
  Galaxy *Psi1Full = (Galaxy *) malloc (sizeof(Galaxy));
  Mat MijBlock;
  Vec f;
  Universe _eqnOfInterestFull = multiquation->equation;
  KSP ksp = NULL;
  PC pc;
  Vec reducedTime;
  Vec uPlus, uNeg;
  Mat rOBStatePlus;
  //---------------------------------------------------------------------------
  // Malloc and initialize here                                             
  //---------------------------------------------------------------------------
  KSPCreate(PETSC_COMM_SELF, &ksp);   
  hrom = reduced->hrom;                                                      
  IsItreduced = reduced->reducedSolution;                                    
  if (reduced->reducedSolution == 1){   //Choose which kind of LSS system      
    _eqnOfInterest = multiquation->equationReduced;  //ROM version
    _solOfInterest = primal->reduced;   //ROM version (ODE)                  
  }else{                                                                       
    _eqnOfInterest = multiquation->equationLSS; //primal can be be the
    _solOfInterest = primal->self;      //Original LSS version (ODE)         
  }                                                                            
  ks_copyUtype(primal->self, Psi1Full);
  createSystem(multiquation->equation, Psi1Full);
  systemReducedSize = _eqnOfInterest.numStates; //Size of the LSS system       
  VecCreateSeq(PETSC_COMM_SELF, reduced->reducedMesh.node.count*
  	       multiquation->equation.numStates, &residual); //Co
  //VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &residual);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &reducedTime);
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &psi);     //Reduced Adjoint
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &truePsi);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &uPlus);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &uNeg);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &adjoint);//Full
  VecCreateSeq(PETSC_COMM_SELF, reduced->nBasisFuncs, &f);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &nSamplefApprox);  
  MatCreateSeqDense(PETSC_COMM_SELF, reduced->nBasisFuncs,
		    primal_h->self->systemSize, NULL, &rOBStatePlus);
  //---------------------------------------------------------------------------
  // Implementation                                                          
  //---------------------------------------------------------------------------
  if (ykflow->numElemCol == 3){
    ks_massBlock(multiquation->equation, primal_h->self, &MijBlock);
    MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
    inverseMatLU(MijBlock, &primal_h->self->Mij);
  }

  moorePenrosePseudoInv(reduced->rOBState, primal_h->self->systemSize,
			reduced->nBasisFuncs, &rOBStatePlus);
     

  for (i=0; i<reduced->nTimeSegs; i++){ //each nu                              
    left = numfiles*i; //Refers to left node of the time segment               
    right = numfiles*(i+1);  //Refers to the right node of the time segment    
    ks_copyUtype(_solOfInterest, Psi1->self);                                  
    Psi1->self->time.t_0 = timelength*i; //Refers to the time t of the left n  
    Psi1->self->time.t_f = timelength*(i+1); //Refers to the time t of the     
    Psi1->self->basis.p = 0;                                                  
    Psi1->self->basis.nodes = 1;
    Psi1->self->space.elem.count = 1;  
    Psi1->self->space.node.count = 1;
    sprintf(Psi1->clusterId, "psi1_seg_%d",i);                                 
    sprintf(Psi1->self->utypeId, "psi1_seg_%d",i);                             
    createSystem(_eqnOfInterest, Psi1->self);                                  
    strcpy(adjointName, Psi1->self->id);                                
    for (j=left; j<right+1; j++){ //iterate through each time segment          
      //-----------------------------------------------------------------------
      // Left Residual Calculation                                           
      //-----------------------------------------------------------------------
      if (j > left){                                                          
	ks_readSolution(_eqnOfInterestFull, primal_H->self, j);
	ykflow->injectH2h(ykflow, multiquation, primal_H->self, primal_h->self,
			  reduced);
	array2Vec(_eqnOfInterestFull, primal_h->self, uPlus);
	ykflow->anyFunction(ykflow, multiquation, primal_h->self, residual,
			    reduced);
        if (reduced->hrom == 1){
	  MatMult(reduced->Zmult, residual, nSamplefApprox);    
	  MatMult(reduced->A, nSamplefApprox, f);  
	}

	ks_readSolution(_eqnOfInterestFull, primal_h->self, j-1);
	array2Vec(_eqnOfInterestFull, primal_h->self, uNeg);

	VecAXPY(uPlus, -1 , uNeg);
	VecScale(uPlus, 1.0/primal_h->self->time.dt);
	MatMult(rOBStatePlus, uPlus, reducedTime);
      }                                                                        
      //-----------------------------------------------------------------------
      // Adjoint Manipulation
      //-----------------------------------------------------------------------
      ks_readSolution(_eqnOfInterest, Psi1->self, j);
      array2Vec(_eqnOfInterest, Psi1->self, psi);
      //-----------------------------------------------------------------------
      // Error Estimation
      //-----------------------------------------------------------------------
      if (j>left){
	VecDot(psi, reducedTime, &timedel);
	VecDot(psi, f, &del);
        if (j == left || j == right)
          del = (del+timedel)*primal->self->time.dt/2.0;
        else
          del = (del+timedel)*primal->self->time.dt;
	//printf("%g\n", del);
	delJ -= del;
	if (j == left+1  && primal->self->time_method == 2) //Us
          timeSolve = 2;
        else if (j == left+1 && primal->self->time_method == 3)
          timeSolve = 2;
        else if (j == left+2  && primal->self->time_method == 3)
          timeSolve = 3; //After two iterations you can now go and use BDF3
        //---------------------------------------------------------------------
        // If use the reduced adjoint, then find the full adjoint form and prit
        //--------------------------------------------------------------------
      }
    }                                                                          
    destroySystem(_eqnOfInterest, Psi1->self);                                 
  }
  //---------------------------------------------------------------------------
  // Destroy everything                                                        
  //---------------------------------------------------------------------------
  destroySystem(multiquation->equation, Psi1Full);  
  VecDestroy(&residual);                                                       
  VecDestroy(&psi);                                                            
  VecDestroy(&truePsi);                                                        
  VecDestroy(&adjoint);                                                        
  free(Psi1->self);                                                            
  free(Psi1);                                                                  
  free(Psi1Full);
  reduced->delJ = delJ;                                                        
  VecDestroy(&f);
  VecDestroy(&uPlus);
  VecDestroy(&uNeg);
  VecDestroy(&reducedTime);
  VecDestroy(&nSamplefApprox);
  MatDestroy(&rOBStatePlus);
  
  if (ykflow->numElemCol == 3){ 
    MatDestroy(&primal_h->self->Mij);
    MatDestroy(&MijBlock);
  }
  if (reduced->hrom == 1)
    MatDestroy(&reduced->A);
}   
