#include "ks_errorFind.h"

/* //Only worry about the full problem for this. Don't worry about the */      
/* //reduced order problem here for now */                                     
double yk_sensitivityLSScalc(Multiverse multiquation, Cluster *primal,         
                         int numTimeSegs){                                     
  int i, j;      //intialization for iteration                                 
  int left, right; //                                                          
  int numfiles = primal->self->time.count/numTimeSegs;                         
  int systemSize = multiquation.equationLSS.numStates;                         
  double dJbards = 0;                                                          
  double aveSpace = 0;                                                         
  double timeFinal = primal->self->time.globalT_f;                             
  double timeInitial=primal->self->time.globalT_0;                             
  double timelength = (timeFinal-timeInitial)/numTimeSegs;                     
  PetscScalar delS = 0;                                                        
  Vec vVec;                                                                    
  Vec f;                                                                       
  PetscScalar numer, denom;                                                    
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));                        
  Is_it *reduced = (Is_it *) malloc (sizeof(Is_it));                           
  Psi1->self = (Galaxy *) malloc (sizeof(Galaxy));                             
  Psi1->state = (Galaxy *) malloc (sizeof(Galaxy));                            
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &vVec);                            
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &f);                               
  ks_copyUtype(primal->self, Psi1->state);                                     
  copySystem(multiquation.equation, primal->self, Psi1->state);                
  for (i=0; i<numTimeSegs; i++){ //each nu                                     
    left = numfiles*i; //Refers to left node of the time segment               
    right = numfiles*(i+1);  //Refers to the right node of the time segment    
    ks_copyUtype(primal->self, Psi1->self);                                    
    Psi1->self->time.t_0 = timelength*i; //Refers to the time t of the left n  
    Psi1->self->time.t_f = timelength*(i+1); //Refers to the time t of the     
    Psi1->self->basis.p = 0;                                                   
    Psi1->self->space.dx = 0;                                                  
    sprintf(Psi1->clusterId, "psi2_seg_%d",i);                                 
    sprintf(Psi1->self->utypeId, "psi2_seg_%d",i);                             
    createSystem(multiquation.equationLSS, Psi1->self);                        
    reduced->reducedSolution = 0;                                              
    reduced->hrom = 0;                                                         
    //-------------------------------------------------------------------------
    // Retrieve the tangent solution for sensitivities                         
    //-------------------------------------------------------------------------
    for (j=left; j<right+1; j++){ //iterate through each time segment          
      //--------------------------------------------------------------------   
      // Left Residual Calculation                                             
      //---------------------------------------------------------------------  
      ks_readSolution(multiquation.equationLSS, Psi1->self, j);                
      array2Vec(multiquation.equationLSS, Psi1->self, vVec);                   
      VecDot(primal->self->djdu, vVec, &delS);                                 
      dJbards += (delS*primal->self->time.dt)/primal->self->time.globalT_f;
      if (j == right+1){                                                       
        ks_readSolution(multiquation.equation, Psi1->state, j);                
        ks_function(multiquation, Psi1, reduced, 0, f);                        
        VecDot(f, vVec, &numer);                                               
        VecDot(f, f, &denom);                                                  
        aveSpace=spatialOutputAverage(multiquation.equation, primal->self,     
                                      j, reduced);                             
        dJbards += (numer*primal->self->j_bar-aveSpace)/                       
          (primal->self->time.globalT_f*denom);                                
      }                                                                        
    }                                                                          
    destroySystem(multiquation.equationLSS, Psi1->self);                       
  }                                                                            
  destroySystem(multiquation.equation, Psi1->state);        
  VecDestroy(&vVec);                                                           
  VecDestroy(&f);                                                              
  free(Psi1->state);                                                           
  free(Psi1->self);                                                            
  free(Psi1);                                                                  
  free(reduced);                                                          
  return dJbards;                                                               
}

//For BDF1 and BDF2 kind of thing
void ks_errorD(Multiverse multiquation, Cluster *primal_h, Cluster *primal_H,
	       Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section       
  //------------------------------------//-------------------------------------
  int i, j, k;                             //initialization for iteration
  int numfiles = primal_h->self->time.count/reduced->nTimeSegs;
  int timeSolve = 1;                    //designates which BDF solver to use
  int left, right;
  int systemSize;
  int systemReducedSize;
  int hrom, IsItreduced;
  char adjointName[1000];
  double timeFinal = primal_h->self->time.globalT_f;
  double timeInitial=primal_h->self->time.globalT_0;
  double timelength = (timeFinal-timeInitial)/reduced->nTimeSegs;
  double delJ = 0;
  PetscInt *index; 
  double del = 0;
  Vec residual;
  Vec psi, truePsi, adjoint;
  Universe _eqnOfInterest;              //Current universe  
  Galaxy *_solOfInterest;                                                  
  Cluster *_primal_hOfInterest = (Cluster *) malloc (sizeof(Cluster));
  _primal_hOfInterest->self = (Galaxy *) malloc (sizeof(Galaxy));
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));
  Psi1->self = (Galaxy *) malloc (sizeof(Galaxy));
  Galaxy *Psi1Full;
  Mat MijBlock;
  
  //---------------------------------------------------------------------------
  // Malloc and initialize here
  //---------------------------------------------------------------------------
  hrom = reduced->hrom;
  IsItreduced = reduced->reducedSolution;
  if (reduced->reducedSolution == 1){   //Choose which kind of LSS system
    _eqnOfInterest = multiquation.equationReduced;  //ROM version
    _solOfInterest = primal_h->reduced;   //ROM version (ODE)
  }else{
    _eqnOfInterest = multiquation.equationLSS; //primal can be be the red
    _solOfInterest = primal_h->self;      //Original LSS version (ODE)
  }
  if (reduced->hrom == 1){               //Choose which full solution
    ks_copyUtype(primal_h->stateFull, _primal_hOfInterest->self);  
    copySystem(multiquation.equation, primal_h->stateFull,
	       _primal_hOfInterest->self);
    Psi1Full = (Galaxy *) malloc (sizeof(Galaxy));
    ks_copyUtype(primal_h->stateFull, Psi1Full);
    createSystem(multiquation.equation, Psi1Full); 	
  }else if (reduced->hrom == 0){
    ks_copyUtype(primal_h->self, _primal_hOfInterest->self);
    copySystem(multiquation.equation, primal_h->self,
	       _primal_hOfInterest->self);
  }
  strcpy(_primal_hOfInterest->clusterId, primal_h->clusterId);
  systemReducedSize = _eqnOfInterest.numStates; //Size of the LSS system
  systemSize = multiquation.equation.numStates* //Size of the full system
    _primal_hOfInterest->self->space.node.count;
  index = (int *) malloc (systemSize *sizeof(int)); //Index for full system
  for (i=0; i<systemSize; i++)
    index[i] = i;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &residual);   //Coarse solution R
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &psi); //Reduced Adjoint
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &truePsi);
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &adjoint);    //Full space adjoint
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Mass Matrix                                                          
  //-------------------------------------------------------------------------- 
  ks_massBlock(multiquation.equation, primal_h->self, &MijBlock);
  MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);
  inverseMatLU(MijBlock, &primal_h->self->Mij);
  for (i=0; i<reduced->nTimeSegs; i++){ //each nu
    left = numfiles*i; //Refers to left node of the time segment
    right = numfiles*(i+1);  //Refers to the right node of the time segment
    ks_copyUtype(_solOfInterest, Psi1->self);
    Psi1->self->time.t_0 = timelength*i; //Refers to the time t of the left n
    Psi1->self->time.t_f = timelength*(i+1); //Refers to the time t of the
    Psi1->self->basis.p = 0;
    Psi1->self->space.dx = 0;
    sprintf(Psi1->clusterId, "psi1_seg_%d",i);
    sprintf(Psi1->self->utypeId, "psi1_seg_%d",i);
    createSystem(_eqnOfInterest, Psi1->self);
    strcpy(adjointName, Psi1->self->id);
    reduced->reducedSolution = 0;
    reduced->hrom = 0;
    for (j=left; j<right+1; j++){ //iterate through each time segment
      if (j > left){
	//---------------------------------------------------------------------
	// Left Residual Calculation
	//---------------------------------------------------------------------
	ks_readSolution(multiquation.equation, primal_H->self, j);
	inject(multiquation.equation,primal_H->self,_primal_hOfInterest->self);
	ks_totalResidual(multiquation, _primal_hOfInterest, residual,
			 reduced, timeSolve);
      }
      //-----------------------------------------------------------------------
      // Adjoint Manipulation
      //-----------------------------------------------------------------------
      ks_readSolution(_eqnOfInterest, Psi1->self, j);
      array2Vec(_eqnOfInterest, Psi1->self, psi);
      if (hrom == 1)
	MatMult(reduced->rOBState, psi, truePsi);
      else
	VecCopy(psi, truePsi);
      MatBlockMult(primal_h->self->Mij, truePsi, adjoint);
      
      VecAssemblyBegin(adjoint);
      VecAssemblyEnd(adjoint);
      if (j>left){
	//---------------------------------------------------------------------
	// Error Estimation
	//---------------------------------------------------------------------
	VecDot(truePsi, residual, &del);
	//       	del = del*_primal_hOfInterest->self->time.dt;
	if ((_primal_hOfInterest->self->time_method == 2
	     || _primal_hOfInterest->self->time_method == 3) &&
	    (j == left || j == right))
	  del = del*_primal_hOfInterest->self->time.dt/2;
	else
	  del = del*_primal_hOfInterest->self->time.dt;
	delJ -= del;
	if (j == left+1  && _primal_hOfInterest->self->time_method == 2) //Us
	  timeSolve = 2;
	else if (j == left+1 && _primal_hOfInterest->self->time_method == 3)
	  timeSolve = 2;
	else if (j == left+2  && _primal_hOfInterest->self->time_method == 3)
	  timeSolve = 3; //After two iterations you can now go and use BDF3
	//---------------------------------------------------------------------
	// If use the reduced adjoint, then find the full adjoint form and prit
	//---------------------------------------------------------------------
      }
      if (hrom == 1){
	sprintf(Psi1Full->id, "%s_full", Psi1->self->id); 
	vec2Array(multiquation.equation, Psi1Full, adjoint);
	ks_printSolution(multiquation.equation, Psi1Full, j);
      }
    }
    destroySystem(_eqnOfInterest, Psi1->self);
  }
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  VecDestroy(&residual);
  VecDestroy(&psi);
  VecDestroy(&truePsi);
  VecDestroy(&adjoint);
  MatDestroy(&MijBlock);
  MatDestroy(&primal_h->self->Mij);
  destroySystem(multiquation.equation, _primal_hOfInterest->self); 
  free(Psi1->self);
  free(Psi1);
  free(_primal_hOfInterest->self);
  free(_primal_hOfInterest);
  free(index);			
  reduced->delJ = delJ;
}



void yk_tradErrorEst(Multiverse multiquation, Cluster *primal_h,
			       Cluster *primal_H, Is_it *reduced){
  int i;
  Vec residual;                                                               
  Vec psi;                                                                   
  double del;                                                                
  int j;
  int timeSolve = 1;
  Galaxy *Psi = (Galaxy *) malloc (sizeof(Galaxy));                   
  double tradDelJ = 0;  
  ks_copyUtype(primal_h->self, Psi);                                    
  sprintf(Psi->utypeId, "psi_unsteady");                                     
  createSystem(multiquation.equation, Psi);                                 
  reduced->hrom = 0;                                                          
  reduced->reducedSolution = 0;                                               
  VecCreateSeq(PETSC_COMM_SELF, 480, &residual);
  VecCreateSeq(PETSC_COMM_SELF, 480, &psi);                                  
  //--------------------------------------------------------------------------
  //
  //--------------------------------------------------------------------------
  if (!reduced->restart){
    
    printf("---------------------------------------------------------\n");     
    printf("|               Traditional Adjoint                     |\n");     
    printf("---------------------------------------------------------\n"); 
    ks_sensitivity(multiquation, primal_h);
    ks_adjoint_Unsteady(multiquation, primal_h);
  }

  for (i=1; i<primal_h->self->time.count+1; i++){                        
    ks_readSolution(multiquation.equation, primal_H->self, i);           
    inject(multiquation.equation, primal_H->self, primal_h->self);       
    ks_totalResidual(multiquation, primal_h, residual, reduced, timeSolve);  
    ks_readSolution(multiquation.equation, Psi, i);                      
    array2Vec(multiquation.equation, Psi, psi);                          
    VecDot(psi, residual, &del);                                              
    if (i==primal_h->self->time.count)                                  
      del = del*primal_h->self->time.dt/2;                              
    else                                                            
      del = del*primal_h->self->time.dt;                                 
    tradDelJ -= del;                                                        
    if (i == 1  && primal_h->self->time_method == 2) //                     
      timeSolve = 2;                                                         
    else if (i == 1 && primal_h->self->time_method == 3)                     
      timeSolve = 2;                                                     
    else if (i == 2  && primal_h->self->time_method == 3)                  
      timeSolve = 3; //After two iterations you can now go and use BDF3    
  }
  VecDestroy(&psi);                                                        
  VecDestroy(&residual);
  destroySystem(multiquation.equation, Psi);
  reduced->tradDelJ = tradDelJ;
}
