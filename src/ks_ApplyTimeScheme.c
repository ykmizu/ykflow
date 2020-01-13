
#include "ks_ApplyTimeScheme.h"
//-----------------------------------------------------------------------------
// Solves the primal solution (u_state)
//
// Yukiko Shimizu
// June 27, 2016
// Error Estimation for Chaotic Systems
//-----------------------------------------------------------------------------

void ks_ApplyTimeScheme(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			Cluster *primal, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section               
  //------------------------------------//-------------------------------------
  int i;                                //initialization for iteration
  int innertimesolver = 1;
  int printnum;
  int left = primal->self->time.t_0/primal->self->time.dt;
  int right = primal->self->time.t_f/primal->self->time.dt;
  Universe _eqnOfInterest;
  //---------------------------------------------------------------------------
  // Initializing things
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)                                           
    _eqnOfInterest = multiquation->equationReduced;
  else if (reduced->reducedSolution == 0 &&                                    
           strncmp(primal->self->utypeId, "psi", 3)==0)                       
    _eqnOfInterest = multiquation->equationLSS;              
  else                                                                         
    _eqnOfInterest = multiquation->equation;      
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=1; i<primal->self->time.count+1; i++){
    if (strncmp(primal->clusterId, "psi1", 4) == 0){ //If we're solving t
      primal->self->time.node = -i + right;                  
      printnum = primal->self->time.node;                    
    }else{                                  //If we're    
      primal->self->time.node = i + left ;                           
      printnum = primal->self->time.node;                 
    }                                                     
    ks_implicitTimeSolve(ykflow, multiquation, primal,reduced,innertimesolver);
    if (i == 1 && primal->self->time_method == 2) //Use BDF1 at the beginnin
      innertimesolver = 2;                                                    
    else if (i == 1 && primal->self->time_method == 3)
      innertimesolver = 2;                                                    
    else if (i == 2 && primal->self->time_method == 3)               
      innertimesolver = 3; //After two iterations you can 
    ks_printSolution(_eqnOfInterest, primal->self, printnum); //print o
  }      
}
