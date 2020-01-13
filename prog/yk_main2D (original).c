

#include "xf_AllStruct.h"                                              
#include "xf_IO.h"                                                     
#include "xf_MPI.h"                                                  
#include "xf_String.h"                                                 
#include "xf_All.h"                                          
#include "xf_Mesh.h"                                               
#include "xf_MeshTools.h"                                              
#include "xf_Param.h"                                                  
#include "xf_Data.h"                                                        
#include "xf_EqnSet.h"                                                
#include "xf_EqnSetHook.h"                                                
#include "xf_Memory.h"                                          
#include "xf_Residual.h" 
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Output.h"
#include "xf_Solver.h"
#include "xf_SolverUnsteady.h"
#include "xf_Arg.h"
#include "xf_AdaptChaos_Common.h"
#include "xf_State.h"
#include "xf_ErrEst.h"
#include "xf_Adapt.h"
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "struct_def.h"                                                      
#include <time.h>                                                            
#include "read_file.h"                                                       
#include "system_initialization.h"                                           
#include "ks_initialConditions.h"                                            
#include "ks_jbaraverage.h"                                                  
#include "ks_sensitivity.h"                                                  
#include "injection.h"                                                       
#include "petsc.h"                                                           
#include "ks_errorEstimate.h"                                                
#include "ks_errorFind.h"                                                  
#include "yk_createHROM.h"                     
#include "slepcsvd.h"                       
#include "mpi.h"                                                          
#include "ks_adjoint_Unsteady.h"
#include "yk_solverCFD.h"
#include "xf_SolverMultiStepStage.h"

void yk_xfData2Dat(yk_PrimalSolver *ykflow, xf_VectorGroup *U,
		   char jobFile[xf_MAXSTRLEN]){
  int i,j;
  int length;
  char stateInput[xf_MAXSTRLEN];
  char newoutFile[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char numFilesString[xf_MAXSTRLEN];
  char dataString[xf_MAXSTRLEN];
  int numFiles = ykflow->primal->self->time.count;
  sprintf(stateInput, "%s_U", jobFile);
  sprintf(numFilesString, "%d", numFiles);
  char * argv[] = {"/Users/yukikoshimizu/xflow/bin/xf_Data2Text", "-inroot",
		   stateInput, "-batch", "0 1",  numFilesString, NULL};
  char * env[] = {NULL};
  int status;
  FILE *dataFile;
  FILE *newFile;
  int count = 0 ;
  pid_t pid = fork();
  int numStates;
  int stateCount = 0;
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver; 
  numStates = xflow->UGFine->Vector[0]->StateRank;
  if (pid == -1){
  }else if (pid > 0){ //If this is the parent process
    waitpid(pid, &status, 0);
  }else{
    execve(argv[0], argv, env);
    _exit(EXIT_FAILURE); //Needed incase execve fails
  }
  for (i=0; i<numFiles+1; i++){
    stateCount = 0;
    for (j=0; j<U->nVector; j++){
      sprintf(outFile, "%s%d_%s.txt", stateInput, i,
	      xfe_VectorRoleName[U->Type[j]]);
      sprintf(newoutFile, "%s_%d.dat", ykflow->primal->self->id, i); 
      dataFile = fopen(outFile, "r");
      newFile = fopen(newoutFile, "w");
      count = 0;
      while ((fgets(dataString, xf_MAXSTRLEN, dataFile))){
	if (dataString[0]!='\n' &&  dataString[0]!='%' && dataString[1]!=' '){
	  length = strlen(dataString);
	  dataString[length-1] = '\0';
	  fprintf(newFile, "%s", dataString);
	  //  fprintf(newFile, "\n");
	  count ++;
	  if (stateCount == numStates-1){
	    fprintf(newFile, "\n");
	    stateCount = -1;
	  }else{
	    	    fprintf(newFile, " ");}
	  stateCount ++;
	}
      }
      fclose(dataFile);
      fclose(newFile);
    }
  }
  //Converts the necessary solutions at particular times from xflow's Data file to ykflow's Dat file.
}


void yk_Init(yk_PrimalSolver *XFLOW,  Is_it *reduced){
  int i;
  int ierr;
  int nTimeStep = 0;
  char newFace[xf_MAXSTRLEN];
  char jobFile[xf_MAXSTRLEN];
  char eqnReduced[xf_MAXSTRLEN];
  char inputFile[xf_MAXSTRLEN];
  int argcIn = 2;
  char *argvIn[] = {"0", "9", NULL};
  yk_Xflow *xflow  = (yk_Xflow *) XFLOW->solver;
  //FOR REDUCED ORDER MODELING AND EVERYTHING RELATED WITH YKFLOW (WHAT A NAME)
  //---------------------------------------------------------------------------
  // Define the multiverse/equation parameters EQN
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));
  strcpy(xflow->ykflow->multiquation->equation.nameEqn, jobFile);
  xflow->ykflow->multiquation->equation.numStates =
    xflow->UGFine->Vector[0]->StateRank;
  //---------------------------------------------------------------------------
  // Initalize Parameters for hyper reduced order modeling INP
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "input", inputFile));
  sprintf(newFace, "%s.inp", inputFile);
  printf("%s\n", newFace);
  argvIn[1] = newFace;
  //Need to convert the xflow parameters into ykflow parameters initU
  read_input(argcIn, argvIn, xflow->ykflow->multiquation->equation, NULL, NULL,
  	     reduced);
  printf("pain %d\n", reduced->nBasisFuncs);
  getchar();
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  sprintf(eqnReduced, "%s%s", jobFile, "Reduced");
  strcpy(xflow->ykflow->multiquation->equationReduced.nameEqn, eqnReduced);
  XFLOW->multiquation->equationReduced.numStates = reduced->nBasisFuncs;
  printf("young %s\n", xflow->ykflow->multiquation->equationReduced.nameEqn);
  getchar();

  //---------------------------------------------------------------------------
  // Define the Cluster properties here
  //---------------------------------------------------------------------------
  initUtype(XFLOW->primal->self);
  XFLOW->primal->self->time_method=2;
  //---------------------------------------------------------------------------
  // Interpolation and quad points and basis shit
  //---------------------------------------------------------------------------
  xflow->ykflow->primal->self->basis.p = xflow->UGFine->Vector[0]->Order[0];


  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",       
                             &nTimeStep));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",
			       &XFLOW->primal->self->time.t_f));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time",     
                             &XFLOW->primal->self->time.t_0));  
  XFLOW->primal->self->time.dt =
    (XFLOW->primal->self->time.t_f-XFLOW->primal->self->time.t_0)/nTimeStep;
  printf("%d\n", nTimeStep);
  printf("%g\n", XFLOW->primal->self->time.t_f);
  printf("%g\n", XFLOW->primal->self->time.t_0);
  printf("shdjahdajsdh %g\n", XFLOW->primal->self->time.dt);
  getchar();
  //---------------------------------------------------------------------------
  //Set the total number of elements
  //---------------------------------------------------------------------------
  XFLOW->primal->self->space.elem.count = 0;
  for (i=0; i<xflow->All->Mesh->nElemGroup; i++)                               
    XFLOW->primal->self->space.elem.count +=
      xflow->All->Mesh->ElemGroup[i].nElem;                  
  //Set the number of basis nodes
  xflow->ykflow->primal->self->basis.nodes =
    (xflow->UGFine->Vector[0]->GenArray)->r/xflow->UGFine->Vector[0]->
    StateRank;
  //Set the total number of nodes                                              
  xflow->ykflow->primal->self->space.node.count =
    XFLOW->primal->self->basis.nodes*XFLOW->primal->self->space.elem.count;  
  //Set the systemSize
  xflow->ykflow->systemSize = (xflow->UGFine->Vector[0]->GenArray)->r*
    xflow->ykflow->primal->self->space.elem.count;
  XFLOW->index = (int *) malloc (XFLOW->systemSize *sizeof(int));
  for (i=0; i<XFLOW->systemSize; i++)
    XFLOW->index[i] = i;

  //clsuterId
  strcpy(XFLOW->primal->clusterId, "state"); 
  strcpy(XFLOW->primal->self->utypeId, "state");
  printf("WRONG MAN %d\n", xflow->ykflow->primal->self->space.node.count);
  getchar();
  
  createSystem(XFLOW->multiquation->equation, XFLOW->primal->self);

  // yk_xfData2Dat(XFLOW, xflow->UGFine, jobFile);
  
  printf("NUM ELEMS %d\n",  xflow->ykflow->primal->self->space.elem.count);
  printf("DOES THIS MATCH %d\n", XFLOW->primal->self->space.elem.count);
  printf("SYSTEMSIZE %d\n",  xflow->ykflow->systemSize);
  printf("ORDER P%d\n",  xflow->ykflow->primal->self->basis.p);
  printf("NUMBER OF TOTAL NODES%d\n",  xflow->ykflow->primal->self->space.node.count);
}

int yk_forwardSimulation(yk_Xflow *xflow){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int ierr;
  int OrderIncrement = 1;
  char jobFile[xf_MAXSTRLEN];
  xf_AdaptInfo AdaptInfo;
  real Time0;
  enum xfe_Bool DoRestart;
  //---------------------------------------------------------------------------
  // Find all the job and equation related information
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));
  xf_Call(xf_ReadAllFromJobFile(jobFile, xfe_True, &xflow->All));
  xf_Call(xf_LoadEqnSetLibrary(xflow->All->EqnSet->EqnSetLibrary)); //Equation
  xf_Call(xf_EqnSetRegister(xflow->All->EqnSet));
  //Primal State (DOUBLE CHECK THIS PART)
  xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart",
			     &DoRestart));
  xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG));
  xf_InitAdaptInfo(&AdaptInfo);
  xf_FillAdaptInfo(xflow->All, xflow->KeyValueArg, &AdaptInfo);
  //---------------------------------------------------------------------------
  // Implement Burn Time to find initial conditions for the coarse solution
  //---------------------------------------------------------------------------
  printf("BURN BURN BURN BURN\n");
  //  xf_Call(xf_RunBurn(xflow->All, &AdaptInfo));
  xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart",
			     &DoRestart));
  xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG));
  //---------------------------------------------------------------------------
  // Inject the new initial conditions from coarse onto fine mesh
  //---------------------------------------------------------------------------
  /* //Create the associated memeory */
  ierr =xf_Error(xf_FindSimilarVectorGroup
  		 (xflow->All, xflow->UG, "UGCoarse", xfe_True, xfe_False,
  		  NULL, &xflow->UGFine, NULL));
  /* //Set that allocated memory space to the vectorgroup of interest */
  ierr = xf_Error(xf_SetVectorGroup(xflow->UG, xfe_Set, xflow->UGFine));
  //Inject from coarse to fine space
  xf_Call(xf_ProjectVectorGroupInPlace_OrderIncrement
  	  (xflow->All->Mesh, NULL, xflow->UGFine, OrderIncrement));
    
  xf_Call(xf_CreateUniformTimeHistData
  	  (xflow->All, AdaptInfo.LogOutput, &xflow->TimeHistData));
  //---------------------------------------------------------------------------
  // Run the forward simulation in the fine space (Basis for the rest of sims)
  //---------------------------------------------------------------------------
  printf("RUNNING RUNNING EXECUTION\n");
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time", &Time0)); 
  /* xf_Call(xf_ApplyTimeScheme */
  /* 	  (xflow->All, AdaptInfo.SavePrefix,xfe_False, &xflow->UGFine, */
  /* 	   xflow->TimeHistData)); */
  printf("i hope it works\n");
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0)); 
  //Convert the output files into ykflow reading format
  //---------------------------------------------------------------------------
  // Destroy and Free anything tht should be dealt with
  //---------------------------------------------------------------------------
  xf_DestroyAdaptInfo(&AdaptInfo, xfe_False);
  return ierr;
}

void yk_RunErrorEstChaos(yk_PrimalSolver *ykflow,Cluster *primalROM_h,
			 Cluster *primalHROM_h, char *argv[], Is_it *reduced){
  printf("head\n");
  printf("tail\n");
  printf("%s\n", ykflow->primal->clusterId);
  printf("%s\n", ((yk_Xflow *)ykflow->solver)->ykflow->primal->clusterId);
  //  primalROM_h->self = (Galaxy *) malloc (sizeof(Galaxy));
  //primalHROM_h->self = (Galaxy *) malloc (sizeof(Galaxy));
  reduced->typeLSS = 1;
 
  if (reduced->reducedLSS){
    printf("----------------------------------------------------------\n");
    printf("|                 Reduced-Order Modeling                 |\n");
    printf("----------------------------------------------------------\n");
    yk_createReducedOrderModel(ykflow, primalROM_h, reduced);
    //if (!reduced->restart)
    yk_runReducedOrderModel(ykflow, primalROM_h, reduced);

  }
}





int main(int argc, char *argv[]){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                
  //------------------------------------//-------------------------------------
  int ierr;
  int myRank, nProc;
  char *ArgIn[] = {"job", "NULL", ".job file name to read (run parameters)",
		   "input", "NULL", ".inp file name to read (run ROM params)", 
                   "BurnTime", "0.0", "burn time before each simulation",
                   "AdaptMode", "Window", "Window/Average/Slope/Residual",
                   "nAdapt", "5", "number of chaos adaptive iterations",
                   "nWindow", "10", "starting number of time windows",     
                   "\0"};     
  Cluster *primalROM_h = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primalHROM_h = (Cluster *) malloc (sizeof(Cluster));
  Is_it *reduced =(Is_it*) malloc (sizeof(Is_it));
  yk_PrimalSolver *XFLOW = new_Xflow();
  yk_Xflow *xflow = (yk_Xflow *) XFLOW->solver;
  //---------------------------------------------------------------------------
  // Intialize everything related to the xflow stuff
  //---------------------------------------------------------------------------
  SlepcInitialize(&argc, &argv, NULL, NULL);                               
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  /* Initialize parallel-run (no effect in serial) */
  xf_Call(xf_MPI_Init(&argc, &argv));
  /* Determine myRank*/
  xf_Call(xf_MPI_GetRank(&myRank, &nProc));
  // initialize key-value
  xf_Call(xf_CreateKeyValue(&xflow->KeyValueArg));
  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, xflow->KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
  //---------------------------------------------------------------------------
  // Search for Attractor/Trajectory Burn Time
  //---------------------------------------------------------------------------
  yk_forwardSimulation(xflow);
  printf("words\n");
  getchar();
  //---------------------------------------------------------------------------
  // HROM and LSS STUFF HERE
  //---------------------------------------------------------------------------
  yk_Init(XFLOW,  reduced);

  printf("pyramid\n");
  yk_RunErrorEstChaos(XFLOW, primalROM_h, primalHROM_h, argv, reduced);
  //---------------------------------------------------------------------------
  // Destroy the Petsc and xflow stuff
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyTimeHistData(xflow->TimeHistData));
  xf_Call(xf_DestroyVectorGroup(xflow->UGFine, xfe_True, xfe_True));
  xf_CloseEqnSetLibrary(&xflow->All);
  xf_Call(xf_DestroyAll(xflow->All));
  xf_Call(xf_DestroyKeyValue(xflow->KeyValueArg));
  MatDestroy(&reduced->rOBState);

  SlepcFinalize(); 
  //---------------------------------------------------------------------------
  // Free and Unmalloc everything
  //---------------------------------------------------------------------------
  destroySystem(XFLOW->multiquation->equation, XFLOW->primal->self);
  destroySystem(XFLOW->multiquation->equation, primalROM_h->self); 
  destroySystem(XFLOW->multiquation->equationReduced,  primalROM_h->reduced);  
  free(XFLOW->index);
  free(primalROM_h->self);
  free(primalROM_h->reduced);
  //free(primalHROM_h->self);
  free(primalROM_h);
  free(primalHROM_h);
  free(reduced);
  free(xflow->ykflow->multiquation);
  free(xflow->ykflow->primal->self);
  free(xflow->ykflow->primal);
  free(xflow->ykflow);
  free(xflow);
}



