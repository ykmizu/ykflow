/*
  =============================================================================
  Name        : Kuramoto-Sivashinsky Implementation 
  Author      : Yukiko Sonya Shimizu
  Version     : February 18, 2016
  Copyright   : ykmizu inc
  Description : Implementaiton of Kuramoto-Sivashinsky PDE(4th order deriv)
  =============================================================================
*/


//Adjoint code: WILL IT WORK!?

#include <stdio.h>                                                        
#include <stdlib.h>                                                     
#include "math.h"                                                          
#include "struct_def.h"
#include <time.h>
#include "read_file.h"
#include "system_initialization.h"
#include "ks_initialConditions.h"
#include "ks_jbaraverage.h"
#include "injection.h"
#include "petsc.h"
#include "ks_errorEstimate.h"
#include "ks_errorFind.h"
#include "yk_createHROM.h"
#include "slepcsvd.h"
#include "mpi.h"
#include "yk_leastSquaresShadowing.h"
#include "ks_adjoint_Unsteady.h"
#include "yk_calculateError.h"


void yk_Init(int argc, char *argv[], yk_PrimalSolver *ykflow,
	     Multiverse *multiquation, Cluster *primal_H, Cluster *primal_h,
	     Is_it *reduced){
  int i;
  primal_h->self=(Galaxy *) malloc (sizeof(Galaxy));
  primal_H->self = (Galaxy *) malloc (sizeof(Galaxy));
  initIs_it(reduced);
  initUtype(primal_H->self);
  initUtype(primal_h->self);
  multiquation->equation = read_eqn(argc, argv);
  read_input(argc, argv, multiquation->equation, primal_h->self, primal_H->self
	     ,reduced);
  //---------------------------------------------------------------------------
  // Reduced multiquation equaion parameters
  //---------------------------------------------------------------------------
  strcpy(multiquation->equationReduced.nameEqn, "meksReduced");
  multiquation->equationReduced.numStates = reduced->nBasisFuncs;
  //---------------------------------------------------------------------------
  //Set U: Coarse solution at lower interpolation number
  //---------------------------------------------------------------------------
  strcpy(primal_H->clusterId, "state");
  strcpy(primal_H->self->utypeId, "state");
  primal_H->self->basis.nodes = primal_H->self->basis.p+1;
  // Need to add space.elem.count too.
  primal_H->self->space.elem.count =
    (primal_H->self->space.x_f-primal_H->self->space.x_0)/
    primal_H->self->space.dx;
  primal_H->self->space.node.count = primal_H->self->space.elem.count*
    primal_H->self->basis.nodes;
  primal_H->self->systemSize = primal_H->self->space.node.count*
    multiquation->equation.numStates;
  createSystem(multiquation->equation, primal_H->self);   //Tr
  ykflow->numElemCol = 3; //1D
  //---------------------------------------------------------------------------
  // Set Ufine: Fine Solution at higher interpolation number
  //---------------------------------------------------------------------------
  //ks_copyUtype(primal_H->self, primal_h->self);
  strcpy(primal_h->self->utypeId, primal_H->self->utypeId);
  primal_h->self->time.dt = primal_H->self->time.dt;                          
  primal_h->self->time.globalT_0 = primal_H->self->time.globalT_0;
  primal_h->self->time.globalT_f = primal_H->self->time.globalT_f;  
  primal_h->self->time.t_0 = primal_H->self->time.t_0;                   
  primal_h->self->time.t_f = primal_H->self->time.t_f;             
  primal_h->self->time.count = primal_H->self->time.count;
  primal_h->self->space.dx = primal_H->self->space.dx;    
  primal_h->self->space.x_0 = primal_H->self->space.x_0;  
  primal_h->self->space.x_f = primal_H->self->space.x_f;        
  primal_h->self->quad.n = primal_H->self->quad.n;  
  primal_h->self->basis.nodes = primal_h->self->basis.p+1;
  
  primal_h->self->space.elem.count = primal_H->self->space.elem.count;
  primal_h->self->space.node.count = primal_h->self->space.elem.count*
    primal_h->self->basis.nodes;
  strcpy(primal_h->clusterId, "state");
  strcpy(primal_h->self->utypeId, "state");
  createSystem(multiquation->equation, primal_h->self);
  primal_h->self->systemSize = primal_h->self->space.node.count*
    multiquation->equation.numStates;
  primal_h->self->index = (int *) malloc (primal_h->self->systemSize
					  *sizeof(int));
  for (i=0; i<primal_h->self->systemSize; i++)
    primal_h->self->index[i] = i;
}

void yk_forwardSolve(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		     Cluster *primal_H, Cluster *primal_h, Is_it *reduced){
  if (!reduced->restart){

    printf("----------------------------------------------------------\n");
    printf("|        Searching for Attractor/Trajectory Burn         |\n");
    printf("----------------------------------------------------------\n");
    ks_findAttractor(ykflow, multiquation, primal_H, reduced);
    ks_readSolution(multiquation->equation, primal_H->self, 0);
    inject(ykflow, multiquation, primal_H->self, primal_h->self, reduced);
    ks_printSolution(multiquation->equation, primal_h->self, 0);
    //-------------------------------------------------------------------------
    // descontinuous Galerkin Solver
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|                       p_H  State Solve                 |\n");
    printf("----------------------------------------------------------\n");
    ks_ApplyTimeScheme(ykflow, multiquation, primal_H, reduced);
    printf("----------------------------------------------------------\n");
    printf("|                         p_h State Solve                 |\n");
    printf("----------------------------------------------------------\n");
    ks_ApplyTimeScheme(ykflow, multiquation, primal_h, reduced);
  }
}
  

void yk_RunErrorEstChaos(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Cluster *primal_h, Cluster *primal_H,
			 Cluster *primalROM_h, Cluster *primalHROM_h,
			 char *argv[], Is_it *reduced){
  if (reduced->reducedLSS){
    //-------------------------------------------------------------------------
    // Reduced-Order Modeling
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|                 Reduced-Order Modeling                 |\n");
    printf("----------------------------------------------------------\n");
    yk_createReducedOrderModel(ykflow, multiquation, primal_h, primalROM_h,
			       reduced);
    if (!reduced->restart)
      yk_runReducedOrderModel(ykflow, multiquation, primal_h, primalROM_h,
    			      reduced);
    ks_jbaraverage(ykflow, multiquation, primalROM_h->self, reduced);
    
    //-------------------------------------------------------------------------
    // Hyper-Reduced_Order Modeling
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|             Hyper-Reduced-Order Modeling               |\n");
    printf("----------------------------------------------------------\n");
    yk_createHyperReducedOrderModel(ykflow, multiquation, primal_h,
    				    primalHROM_h, reduced);
    if (!reduced->restart)
      yk_runHyperReducedOrderModel(ykflow, multiquation, primal_h,
      				   primalHROM_h, reduced);
    ks_jbaraverage(ykflow, multiquation, primalHROM_h->stateFull, reduced);
    //-------------------------------------------------------------------------
    //reducced Least Square Shadowing Acted upon the Primal Solution
    //-------------------------------------------------------------------------
    printf("---------------------------------------------------------\n");
    printf("|              HROM-Least Squares Shadowing             |\n");
    printf("---------------------------------------------------------\n");
    reduced->reducedSolution = 1;
    yk_leastSquaresShadowing(ykflow, multiquation, primalHROM_h, argv,reduced);
    yk_estimateError(ykflow, multiquation, primal_h, primal_H, primalHROM_h,
    		     reduced);
  }else if(!reduced->reducedLSS){
    printf("---------------------------------------------------------\n");
    printf("|                Least Squares Shadowing                |\n");
    printf("---------------------------------------------------------\n");
    multiquation->equationLSS.numStates = primal_h->self->space.node.count;
    strcpy(multiquation->equationLSS.nameEqn, "meks");
    reduced->hrom = 0;
    reduced->reducedSolution = 0;
    yk_leastSquaresShadowing(ykflow, multiquation, primal_h, argv, reduced);
    yk_estimateError(ykflow, multiquation, primal_h, primal_H, primal_h,
  		     reduced);
  }
}

void yk_actualError(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *primal_H, Cluster *primal_h, Is_it *reduced){
  ks_jbaraverage(ykflow, multiquation, primal_H->self, reduced);
  ks_jbaraverage(ykflow, multiquation, primal_h->self, reduced);
  reduced->actualJbar = primal_H->self->j_bar - primal_h->self->j_bar;
 }

void yk_Finalize(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		 Cluster *primal_H, Cluster *primal_h, Cluster *primalROM_h,
		 Cluster *primalHROM_h, Is_it *reduced){
  printf("---------------------------------------------------------\n");
  printf("|                        Results                        |\n");
  printf("---------------------------------------------------------\n");
  printf("Actual Error = %0.16e\n", reduced->actualJbar);
  printf("HROM/LSS Estimated Error Psi = %0.16e\n", reduced->delJ_psi);
  printf("HROM/LSS Estimated Error dPsi = %0.16e\n", reduced->delJ_dpsi);
  printf("HROM/LSS Estimated Error dPsi H1= %0.16e\n", reduced->delJ_dpsi_H1);
  if (reduced->reducedLSS){
    MatDestroy(&reduced->rOBState);
    MatDestroy(&reduced->rOBResidual);
    MatDestroy(&reduced->rOBJacobian);
    MatDestroy(&reduced->rOBStateBar);
    MatDestroy(&reduced->rOBStateRed);
    MatDestroy(&reduced->Z);
    MatDestroy(&reduced->Zu);
    MatDestroy(&reduced->Zp);
    MatDestroy(&reduced->Zmult);
    MatDestroy(&reduced->A);
    MatDestroy(&reduced->B);
    destroySystem(multiquation->equationReduced, primalROM_h->reduced);
    destroySystem(multiquation->equationReduced, primalHROM_h->reduced);
    destroySystem(multiquation->equation, primalHROM_h->stateFull);
    destroySystem(multiquation->equation, primalROM_h->self);
    destroySystem(multiquation->equation, primalHROM_h->self); 
    free(primalROM_h->self->index);
    free(primalROM_h->self);
    free(primalROM_h->reduced);
    free(primalHROM_h->self);
    free(primalHROM_h->reduced);
    free(primalHROM_h->stateFull->index);
    free(primalHROM_h->stateFull);
  }
  delArrayInt(&reduced->reducedMesh.elem);
  delArrayInt(&reduced->reducedMesh.node);
  destroySystem(multiquation->equation, primal_H->self);
  destroySystem(multiquation->equation, primal_h->self);
  free(multiquation->equation.c); //Destroy the Universe
  free(primal_h->self->index);
  free(primal_h->self);
  free(primal_H->self);
  free(reduced);
  free(primalHROM_h);
  free(primalROM_h);
  free(primal_h);
  free(primal_H);
  free(multiquation);
}
int main(int argc, char *argv[]){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  double cpu_time_used;
  clock_t start, end;
  Multiverse *multiquation = (Multiverse *) malloc (sizeof(Multiverse));;
  Cluster *primal_H = (Cluster *) malloc (sizeof(Cluster));  //Have t
  Cluster *primal_h = (Cluster*) malloc (sizeof(Cluster));  //Have to malloc
  Cluster *primalROM_h = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primalHROM_h = (Cluster *) malloc (sizeof(Cluster));
  Is_it *reduced = (Is_it *) malloc (sizeof(Is_it));
   yk_PrimalSolver *ykflow = new_Ykflow();
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
   SlepcInitialize(&argc, &argv, NULL, NULL);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  yk_Init(argc, argv, ykflow, multiquation, primal_H, primal_h, reduced);
  start = clock();
  //---------------------------------------------------------------------------
  // Primal Solver and Forward Solution Setup
  //---------------------------------------------------------------------------
  yk_forwardSolve(ykflow, multiquation, primal_H, primal_h, reduced); //Prim
  //---------------------------------------------------------------------------
  // Output-Based Error Estimation for Chaos
  //---------------------------------------------------------------------------
  yk_actualError(ykflow, multiquation, primal_H, primal_h, reduced);
  yk_RunErrorEstChaos(ykflow, multiquation, primal_h, primal_H,
  		      primalROM_h, primalHROM_h, argv, reduced);
  //---------------------------------------------------------------------------
  // Finalize Everything
  //---------------------------------------------------------------------------
  yk_Finalize(ykflow, multiquation, primal_H, primal_h, primalROM_h,
  	      primalHROM_h, reduced);
  ykflow->delete(ykflow);
  SlepcFinalize();
  //---------------------------------------------------------------------------
  // Print Simulation Time
  //---------------------------------------------------------------------------
  end = clock();
  cpu_time_used = ((double) (end-start)) /CLOCKS_PER_SEC;
  printf("Time it took to run simulation = %0.16f\n", cpu_time_used);
  //---------------------------------------------------------------------------
  // Free And Destroy Everything
  //---------------------------------------------------------------------------
  return 0;
  
}
