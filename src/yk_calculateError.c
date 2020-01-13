

#include "yk_calculateError.h"
//will this work

void yk_estimateError(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *primal_h, Cluster *primal_H, Cluster *primal,
		      Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i, j;                             //initialization for iteration  
  int _size;
  int numfiles = primal->self->time.count/reduced->nTimeSegs;
  int timeSolve = 1;                    //designates which BDF solver to use   
  int left, right;                                                             
  int systemReducedSize;                                                       
  int hrom, IsItreduced;                                                       
  char adjointName[1000];      
  double timeFinal = primal->self->time.globalT_f; 
  double timeInitial = primal->self->time.globalT_0;
  double timelength = (timeFinal-timeInitial)/reduced->nTimeSegs;              
  double delJ_psi= 0;
  double delJ_dpsi = 0;
  double delJ_dpsi_H1 = 0;
  double del_psi= 0;
  double del_dpsi = 0;
  double del_dpsi_H1 = 0;
  Vec residual, residualReduced;
  Vec residualMinv, residualSamples;
  Vec psi, dpsi, dpsi_H1;                                                   
  Vec psifull, psihfinal, psihfinal_H1;
  Mat I_proj, I_inject, I_proj_H1, I_inject_H1;
  Mat Idel, Idel_H1, Imult;
  Mat MijBlock;
  Universe _eqnOfInterest;              //Current universe                     
  Universe _eqnOfInterestFull = multiquation->equation;
  Galaxy *_solOfInterest;
  Galaxy *Psi1Full = (Galaxy *) malloc (sizeof(Galaxy));
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));                        
  Psi1->self = (Galaxy *) malloc (sizeof(Galaxy));                             
  //---------------------------------------------------------------------------
  // Malloc and initialize here                                             
  //---------------------------------------------------------------------------
  hrom = reduced->hrom;                                                  
  IsItreduced = reduced->reducedSolution;                                    
  if (reduced->reducedSolution == 1){   //Choose which kind of LSS system      
    _eqnOfInterest = multiquation->equationReduced;  //ROM version
    _solOfInterest = primal->reduced;   //ROM version (ODE)                  
    _size = reduced->nBasisFuncs;
  }else{                                                                       
    _eqnOfInterest = multiquation->equationLSS; //primal can be be the
    _solOfInterest = primal->self;      //Original LSS version (ODE)         
    _size = primal_h->self->systemSize;
  }                                                                            
  ks_copyUtype(primal->self, Psi1Full);
  createSystem(multiquation->equation, Psi1Full);
  systemReducedSize = _eqnOfInterest.numStates; //Size of the LSS system       
  VecCreateSeq(PETSC_COMM_SELF, _size, &residualReduced); //Co
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &residual);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &residualMinv);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &residualSamples);
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &psi);     //Reduced Adjoint
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &dpsi);    //Reduced Adjoint
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &dpsi_H1); //Reduced Adjoint
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &psifull);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &psihfinal);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &psihfinal_H1);
  MatCreateSeqDense(PETSC_COMM_SELF, primal_h->self->basis.p+1,
  		    primal_h->self->basis.p+1, NULL, &Imult);
  //---------------------------------------------------------------------------
  // Implementation                                                          
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  //---------------------------------------------------------------------------
  // Original Psih-PsiH Idel calculation
  //---------------------------------------------------------------------------
  yk_Identity(&Idel, primal_h->self->basis.p+1);
  // Original Projection and injection
  yk_project_I(multiquation,primal_h->self,primal_H->self, &I_proj);
  yk_project_I(multiquation,primal_H->self,primal_h->self, &I_inject);
  //Calcualgtion of I-injection*projection
  MatMatMult(I_inject, I_proj, MAT_REUSE_MATRIX, PETSC_DEFAULT, &Imult);
  MatAXPY(Idel, -1, Imult, SAME_NONZERO_PATTERN);
  //---------------------------------------------------------------------------
  // H1 Projection Psih-PsiH Idel_H1 calculation
  //---------------------------------------------------------------------------
  yk_Identity(&Idel_H1, primal_h->self->basis.p+1);
   //H1 Projection and Injection
  yk_projectH1_I(multiquation, primal_h->self, primal_H->self, &I_proj_H1);
  yk_projectH1_I(multiquation, primal_H->self, primal_h->self, &I_inject_H1);
  //Calculation of I-injection*projection
  MatMatMult(I_inject_H1, I_proj_H1, MAT_REUSE_MATRIX, PETSC_DEFAULT, &Imult);
  MatAXPY(Idel_H1, -1, Imult, SAME_NONZERO_PATTERN);
  //---------------------------------------------------------------------------
  // Mass Matrix Calculation
  //---------------------------------------------------------------------------
  if (ykflow->numElemCol == 3){
    ks_massBlock(multiquation->equation, primal_h->self, &MijBlock);
    MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
    inverseMatLU(MijBlock, &primal_h->self->Mij);
  }
  //---------------------------------------------------------------------------
  // Error Rstimation
  //---------------------------------------------------------------------------
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
      if (j > 0){
  	ks_readSolution(_eqnOfInterestFull, primal_H->self, j);
  	ykflow->injectH2h(ykflow, multiquation, primal_H->self, primal_h->self,
  			  reduced);
  	primal_h->self->time.node = primal_H->self->time.node;
  	ykflow->Residual(ykflow, multiquation, primal_h, residual, NULL,
  			 reduced, timeSolve);
  	if (ykflow->numElemCol == 3){
  	  MatBlockMult(primal_h->self->Mij, residual, residualMinv);
  	}else{
  	  vec2Array(_eqnOfInterestFull, primal_h->self, residual);
  	  ykflow->multiplyByInvMass(ykflow, primal_h->self, residualMinv,
  				    reduced);
  	}
  	if (IsItreduced == 1){
  	  MatMultTranspose(reduced->rOBState, residualMinv, residualReduced);
	  //  	  MatMult(reduced->Z, residualMinv, residualSamples);
  	  //MatMult(reduced->A, residualSamples, residualReduced);
  	}else{
  	  VecCopy(residualMinv, residualReduced);
  	}
  	//---------------------------------------------------------------------
  	// Adjoint Manipulation
  	//---------------------------------------------------------------------
  	ks_readSolution(_eqnOfInterest, Psi1->self, j);
  	array2Vec(_eqnOfInterest, Psi1->self, psi);
	if (IsItreduced == 1)
	  MatMult(reduced->rOBState, psi, psifull);
	else
	  VecCopy(psi, psifull);
        MatBlockMult(Idel, psifull, psihfinal);
	MatBlockMult(Idel_H1, psifull, psihfinal_H1);
	if (IsItreduced == 1){
	  MatMultTranspose(reduced->rOBState, psihfinal, dpsi);
	  MatMultTranspose(reduced->rOBState, psihfinal_H1, dpsi_H1);
	}else{
	  VecCopy(psihfinal, dpsi);
	  VecCopy(psihfinal_H1, dpsi_H1);
	}
	//---------------------------------------------------------------------
  	// Error Estimation
  	//---------------------------------------------------------------------
  	VecDot(psi, residualReduced, &del_psi);            //Psi Error Estimate
	VecDot(dpsi, residualReduced, &del_dpsi);      //delPsi Error Estimate
	VecDot(dpsi_H1, residualReduced, &del_dpsi_H1);//delPsi h1 error estim
	
  	if (j == left || j == right){
          del_psi *= primal->self->time.dt/2.0;
	  del_dpsi *= primal->self->time.dt/2.0;
	  del_dpsi_H1 *= primal->self->time.dt/2.0;
	}else{
  	  del_psi *= primal->self->time.dt;
	  del_dpsi *= primal->self->time.dt;
	  del_dpsi_H1 *= primal->self->time.dt;
	}
	delJ_psi -= del_psi;
	delJ_dpsi -= del_dpsi;
	delJ_dpsi_H1 -= del_dpsi_H1;
	//Assign the next time method for BDF
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
  // Assigned all the error estiamtes to the reduced structure
  //---------------------------------------------------------------------------
  reduced->delJ_psi = delJ_psi;
  reduced->delJ_dpsi = delJ_dpsi;
  reduced->delJ_dpsi_H1 = delJ_dpsi_H1;
  //---------------------------------------------------------------------------
  // Destroy everything                                                        
  //---------------------------------------------------------------------------
  destroySystem(multiquation->equation, Psi1Full);  
  VecDestroy(&residual);
  VecDestroy(&residualSamples);
  VecDestroy(&residualMinv);
  VecDestroy(&residualReduced);
  VecDestroy(&psi);
  VecDestroy(&dpsi);
  VecDestroy(&psifull);
  VecDestroy(&psihfinal);
  VecDestroy(&dpsi_H1);
  VecDestroy(&psihfinal_H1);
  free(Psi1->self);                                                            
  free(Psi1);                                                                  
  free(Psi1Full);
  MatDestroy(&I_proj);
  MatDestroy(&I_inject);
  MatDestroy(&I_proj_H1);
  MatDestroy(&I_inject_H1);
  MatDestroy(&Idel_H1);
  MatDestroy(&Idel);
  MatDestroy(&Imult);
  if (ykflow->numElemCol == 3){ 
    MatDestroy(&primal_h->self->Mij);
    MatDestroy(&MijBlock);
  }
}   
