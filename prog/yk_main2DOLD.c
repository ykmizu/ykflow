
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

//Structure used for storing for xflow 
typedef struct yk_Xflow{                                                       
  yk_PrimalSolver *ykflow;                                                     
  int man;                                                                     
  xf_KeyValue *KeyValueArg;                                                    
  xf_All *All;                                                                 
  xf_UnsteadyData *UnsteadyData;                                               
  xf_VectorGroup *UG;                                                          
  xf_VectorGroup *UGFine;                                                      
  xf_VectorGroup *UGCoarse;                                                    
  xf_TimeHistData *TimeHistData;                                               
  xf_Data *D;                                                                  
}yk_Xflow; 

void yk_findBoundarySeeds(yk_PrimalSolver *primalSolverObj, int nSampNodes,
			  int *numSeeds, int *nodeSet);

void yk_xflow_totalResidual(yk_PrimalSolver *primalSolverObj, Cluster *object, 
                            Vec residual, Mat Jacobian, Is_it *reduced,       
                            int timeSolveNum);

void yk_findAdjacentElems(yk_PrimalSolver *solverObj, int elem, int *nodeSet);

yk_PrimalSolver* new_Xflow(){                                                 
  yk_Xflow * xflow = (yk_Xflow *) malloc (sizeof(yk_Xflow));                   
  yk_PrimalSolver* ykflow = new_Ykflow();
  ykflow->solver = xflow;                                                      
  xflow->ykflow = ykflow;                                                      
  ykflow->Residual = yk_xflow_totalResidual;                                   
  ykflow->boundarySeeds = yk_findBoundarySeeds;
  ykflow->adjacentElems = yk_findAdjacentElems;
  return ykflow;                                                               
}

//-----------------------------------------------------------------------------
//Convert the ykflow array format to xflow Vector format
//-----------------------------------------------------------------------------
void yk_array2Vector(xf_Vector *A, Galaxy *object){                            
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j, k, n;                                                              
  xf_GenArray *genA;                                                           
  int rOrd;  //Number of total basis nodes and states                          
  PetscInt gAInd;                                                              
  int count = 0;                                                               
  int numStates;                                                               
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  for (i=0; i<A->nArray; i++){                                                 
    genA = A->GenArray + i;                                                    
    for (j=0; j<genA->n; j++){ //Iterate through number of elements            
      rOrd = ((genA->vr == NULL) ? genA->r : genA->vr[j]);                     
      numStates = rOrd/object->basis.nodes;                                    
      if (A->GenArray[i].rValue != NULL)                                  
        for (k=0; k<object->basis.nodes; k++)      
          for (n=0; n<numStates; n++)        
            A->GenArray[i].rValue[j][(k*numStates)+n]=                         
              object->solution[n].array[count+j*object->basis.nodes+k];
      //else if (A->GenArray[i].iValue != NULL)   
    }                                                                          
    count += genA->n*object->basis.nodes;                                      
  }                                                                            
} 

//-----------------------------------------------------------------------------
//Convert the ykflow array format to Vector Group
//-----------------------------------------------------------------------------
void yk_array2VectorGroup(xf_VectorGroup *A, Galaxy *object){                
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i;                                                                       
  for (i=0; i<A->nVector; i++)                                                 
    yk_array2Vector(A->Vector[i], object);                                     
}                                                                              

//-----------------------------------------------------------------------------
// Setting a Mat Petsc Object xflow matrix
//-----------------------------------------------------------------------------
void yk_matrix2D2MatGroup(xf_All *All, xf_JacobianMatrix *dRdU, Mat matObj,
			  Is_it *reduced){
  int i, j, k, m, n, p;;
  int groupN, elemN, faceN;
  int negrp = dRdU->negrpself; //Number of groups
  int nface;
  int nn, nnN;
  real *A;
  int sr = All->EqnSet->StateRank;
  int *egrpOffset;
  int globalIndex;
  int r, rN, row, col;
  //---------------------------------------------------------------------------
  // Calculate the global index offsets since there are multiple groups
  //---------------------------------------------------------------------------
  xf_Error(xf_Alloc( (void **) &egrpOffset, negrp, sizeof(int)));   
  // loop over element groups, calculate global index offsets                 
  globalIndex = 0; 
  for (i=0; i<negrp; i++){                            
    egrpOffset[i] = globalIndex;                       
    for (j=0; j<dRdU->nElem[i]; j++){                  
      nn = xf_Jacobian_n(dRdU, i, j);                      
      globalIndex += nn*sr;                                
    }              
  }    
  //---------------------------------------------------------------------------
  // Find the spcieifc values
  //---------------------------------------------------------------------------
  for (i=0; i<negrp; i++){ //Loop through number of groups
    for (j=0; j<dRdU->nElem[i]; j++){ // iterate through the number of eleme
      nface = dRdU->nFace[i][j]; //gets the total number of faces
      for (k=-1; k<nface; k++){ //iterate through the number of
	if (k == -1){ //Itself inofmration
	  groupN = i;
	  elemN = j;
	  faceN = k;
	}else{
	  groupN = dRdU->egrpN[i][j][k];
	  elemN = dRdU->elemN[i][j][k];
	  faceN = dRdU->faceN[i][j][k];
	}
	if (groupN <0 ) continue;
	A = dRdU->Value[i][j][1+k];
	nn = xf_Jacobian_n(dRdU, i, j); //row element
	r = sr*nn;
	nnN = xf_Jacobian_n(dRdU, groupN, elemN); //column element
	rN = sr*nnN;
	row = egrpOffset[i];
	col = egrpOffset[groupN];
	for (m = 0; m<j; m++)
	  row +=sr*xf_Jacobian_n(dRdU, i, m);
	for (m = 0; m<elemN; m++)
	  col +=sr*xf_Jacobian_n(dRdU, groupN, m);
	for (m=0; m<nn; m++)
	  for (n=0; n<nnN; n++)
	    for (p=0; p<sr*sr; p++)
	      if (A[(m*nnN+n)*sr*sr+p]!= 0.0)
		MatSetValue(matObj, row+m*sr+(p/sr), col+n*sr+(p%sr),
			    A[(m*nnN+n)*sr*sr+p], INSERT_VALUES);	
      }
    }
  }
  MatAssemblyBegin(matObj, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matObj, MAT_FINAL_ASSEMBLY);
  xf_Release(egrpOffset);  
}
                                                                              
//-----------------------------------------------------------------------------
// Concert vector xflow format to Vec Petsc format                             
//-----------------------------------------------------------------------------
void yk_vector2Vec(xf_Vector *A, Vec vecObj){                                  
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k;                          // initialization for iteration        
  int count = 0;                                                               
  int rOrd;                                                                    
  xf_GenArray *genA;                                                           
  PetscScalar value;                                                           
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------

  for (i=0; i<A->nArray; i++){                                                 
    genA = A->GenArray + i;                                                    
    for (j=0; j<genA->n; j++){ //Iterate through number of elements            
      rOrd = ((genA->vr == NULL) ? genA->r : genA->vr[j]);                     
      if (A->GenArray[i].rValue != NULL){                                      
        for (k=0; k<rOrd; k++){                                                
          value = A->GenArray[i].rValue[j][k];                                 
          VecSetValue(vecObj, count+j*rOrd+k, value, INSERT_VALUES);           
        }                                                                      
      }else if (A->GenArray[i].iValue != NULL){                                
      }                                                                        
    }                                                                          
    count += genA->n*rOrd;                                                     
  }                                                                            
}                                                                              

//-----------------------------------------------------------------------------
// Convert vector xflow group format to Vec Petsc format
//-----------------------------------------------------------------------------
void yk_vector2VecGroup(xf_VectorGroup *A, Vec vecObj){
  int i;
  for (i=0; i<A->nVector; i++)
    yk_vector2Vec(A->Vector[i], vecObj);
}

//-----------------------------------------------------------------------------
// COnvert xflow text file to ykflow state file 
//-----------------------------------------------------------------------------
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
  //Converts the necessary solutions at particular times from xflow's Data fil
  //e to ykflow's Dat file.
}

void yk_findAdjacentElems(yk_PrimalSolver *solverObj, int elem, int *elemSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i, j;                             // initialization for iteration
  int ierr;
  int floorElem = 0;
  int nFace;
  int elemN;
  int newElem;
  enum xfe_Bool Found; 
  xf_JacobianMatrix *dRdU;              // structure for jacobian matrix
  yk_Xflow *xflow = (yk_Xflow *) solverObj->solver; //Set pointer to xflow 
  //---------------------------------------------------------------------------
  // Find the structre of the Jacobian matrix. Will not be calculating values
  //---------------------------------------------------------------------------
  xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, xflow->UGFine,
				NULL, xfe_True, NULL, NULL, &dRdU, &Found));
  //elem: current element of interest, want to find its adjacent elements
  //Need to translate the elem into its element number relative to its group #
  for (i=0; i<dRdU->negrp; i++){
    floorElem += dRdU->nElem[i]; //Calculates the number of elements
    if (elem <= floorElem-1){
      if (i>0)
	newElem = floorElem-elem-1;
      else
	newElem = elem;
      printf("group %d\n", i);
      printf("everyone makes mistakes %d\n", dRdU->nFace[i][newElem]);
      nFace = dRdU->nFace[i][newElem];
      //elem belongs to the ith element group
      printf("%d\n", elem);
      elemSet[elem] = 1;
      for (j=0; j<nFace; j++){
	elemN = dRdU->elemN[i][newElem][j];
	printf("%d\n", elemN);
	if (elemN>-1) //Don't count the boundary elements
	  elemSet[elemN] = 1;
      }
      floorElem = 0;
      break;
    }
  }
}
//Do by elements,use dRdU value use the -1 values
//xf_MeshStruct
//BoundaryFace Group nBFaceGroup
//Inside BFaceGroup list of boundary faces xf_BFGroup Elem xf_BFace<- elem on the boundary
//Redo Redo Redo Redo Redo Redo Redo




void yk_findBoundarySeeds(yk_PrimalSolver *primalSolverObj, int nSampNodes,    
                          int *numSeeds, int *nodeSet){                        
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j, k;                          //initialization for iteration
  int count = 0;                        //will give me total # of elems
  int elem;                             //saves the current element to add
  int numBFG;                           //total num of boundary group 
  int numBasis = primalSolverObj->primal->self->basis.nodes;
  int elemSize = numBasis*primalSolverObj->multiquation->equation.numStates;
  yk_Xflow *xflow = (yk_Xflow *) primalSolverObj->solver;                      
  xf_BFaceGroup *aBFG = xflow->All->Mesh->BFaceGroup;
  xf_ElemGroup *elemGroup = xflow->All->Mesh->ElemGroup;                       
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  numBFG = xflow->All->Mesh->nBFaceGroup; //total num of boundary group
  for (i=0; i<numBFG; i++){                                                    
    for (j=0; j<aBFG[i].nBFace; j++){   //iteration through num of bound elems
      elem = 0;                         //set to 0, want global elem number
      for (k=0; k< aBFG[i].BFace[j].ElemGroup; k++) //elemgroup              
        elem += elemGroup[k].nElem;     //find global number
      elem += aBFG[i].BFace[j].Elem;    //find actual elem value
      for (k=0; k<elemSize; k++)        //Iterate through total elem vec size 
	nodeSet[elemSize*count+k] = elemSize*elem+k; //node vector value
      count++;                          //sums up tot number of elems thus far
    }                                                                          
  }
  if (nSampNodes <= count*elemSize){    //Verfiy there's enough sample nodes
    printf("Need to specify at least %d more nSampleNodes\n", count*elemSize-
	   nSampNodes);        
    exit(EXIT_FAILURE);                 //quit if user needs to specify more
  }else{
    *numSeeds = count*elemSize;         //set the number of seeds to total node
  }
}  

//-----------------------------------------------------------------------------
// Residual Calculation for xflow
//-----------------------------------------------------------------------------
void yk_xflow_totalResidual(yk_PrimalSolver *primalSolverObj, Cluster *object, 
                            Vec residual, Mat Jacobian, Is_it *reduced,        
                            int timeSolveNum){                                 
  int ierr;
  enum xfe_Bool Found;
  int totNodes = primalSolverObj->primal->self->space.node.count*
    primalSolverObj->multiquation->equation.numStates;
  yk_Xflow *xflow = (yk_Xflow *) primalSolverObj->solver;                      
  xf_JacobianMatrix *dRdU;                                                     
  xf_VectorGroup *RG = NULL;                                                   
  xf_VectorGroup *UGi[2];                                                      
  xf_SolverData *SolverData=NULL;                                              
  Galaxy *postTime = NULL;
  Mat ZpT = NULL;
  Vec residualTemp;
  Mat dRdUTemp;
  Mat dRdUApprox = NULL;
  Cluster *temp = (Cluster *) malloc (sizeof(Cluster));                        
  temp->self = (Galaxy *) malloc (sizeof(Galaxy));
  Universe _equation = primalSolverObj->multiquation->equation;
  //need to go in here and modify it for the hyper reduced case
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Set up the State temp value to be used in the Residual                    
  //-------------------------------------------------------------------------- 
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UGFine, "state1",
				     xfe_True, xfe_True, NULL, &UGi[0], NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UGFine, "state2",
				     xfe_True, xfe_True, NULL, &UGi[1], NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UGFine, "RG",           
                                    xfe_True, xfe_True, NULL, &RG, NULL));     
  xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, xflow->UGFine,
                                NULL, xfe_True, NULL, NULL, &dRdU, &Found));   
  //create/allocate SOlverData 
  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));                       
  SolverData->DoProcessResidual = xfe_False; // do not process residual        
  SolverData->UseArtificialDt = xfe_False;                                     
  SolverData->c = 0.;
  //---------------------------------------------------------------------------
  //Create a temporary Cluster called temp to save timenode - 1 values         
  //---------------------------------------------------------------------------
  //Left side of the element
  if (reduced->hrom == 0)
    postTime = object->self;
  else if (reduced->hrom == 1){
    object->stateFull->time.node = object->self->time.node;
    postTime = object->stateFull;
    VecCreate(PETSC_COMM_SELF, &residualTemp);
    VecSetSizes(residualTemp, totNodes, totNodes);        
    VecSetType(residualTemp, VECSEQ);

    /* MatCreate(PETSC_COMM_SELF, &dRdUApprox); */
    /* MatSetSizes(dRdUApprox, totNodes, */
    /* 		_equation.numStates*object->self->space.node.count, */
    /* 		totNodes,_equation.numStates*object->self->space.node.count); */
    /* printf("KIMI TACHI %d\n",  _equation.numStates*reduced->reducedMesh.node.count); */
    /* MatSetType(dRdUApprox, MATSEQDENSE); */
    /* MatSetUp(dRdUApprox); */

    MatCreate(PETSC_COMM_SELF, &dRdUTemp);                                  
    MatSetSizes(dRdUTemp, totNodes, totNodes, totNodes, totNodes);
    MatSetType(dRdUTemp, MATSEQDENSE);
    MatSetUp(dRdUTemp);
    //MatSeqBAIJSetPreallocation(dRdUTemp, object->stateFull->basis.nodes*4,4,NULL); 
    
  }
  //remeber that stateFull is same size as self for hrom == 0
  ks_copyUtype(postTime, temp->self);
  copySystem(_equation, postTime, temp->self);
  temp->self->time.node = postTime->time.node - 1;
  strcpy(temp->clusterId, object->clusterId);
  ks_readSolution(_equation, temp->self, temp->self->time.node);
  yk_array2VectorGroup(UGi[0], temp->self);
  yk_array2VectorGroup(UGi[1], postTime);
  //WAIT I CAN REWRITE THIS BEETER
  //---------------------------------------------------------------------------
  // Residual and Jacobian                                                     
  //---------------------------------------------------------------------------
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time",
			     postTime->time.node*postTime->time.dt));
  //Calculate the Spatial Residual
  xf_Call(xf_CalculateResidual(xflow->All, UGi[1], RG, dRdU, SolverData));    
  //Temporary calculate for du u(n+1)-u(n) for now.... since eval acting weird
  xf_VectorGroupMultSet(UGi[0], -1.0, xfe_Add, UGi[1]);  
  //Calculate the total Unsteady Residual by adding the Mass Matrix Mdu/dt+R
  xf_AddMassMatrixGroup(xflow->All, 1/postTime->time.dt, NULL, UGi[1], RG,
			dRdU, SolverData);
  if (reduced->hrom == 0){
    yk_vector2VecGroup(RG, residual);
    yk_matrix2D2MatGroup(xflow->All, dRdU, Jacobian);
  }else if (reduced->hrom == 1){
    //residual
    yk_vector2VecGroup(RG, residualTemp);
    MatMult(reduced->Zu, residualTemp, residual);
    //Jacobian
    printf("WOWOWOWOW\n");
    yk_matrix2D2MatGroup(xflow->All, dRdU, dRdUTemp);
    printf("YAHAHHA\n");
    MatTranspose(reduced->Zp, MAT_INITIAL_MATRIX, &ZpT);
    MatMatMult(dRdUTemp, ZpT, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
    	       &dRdUApprox);
    //MatMatMultBAIJ(dRdUTemp, ZpT, dRdUApprox); 
    //MatView(dRdUApprox, PETSC_VIEWER_STDOUT_SELF);
    /* MatMatTransposeMult(dRdUTemp, reduced->Zp,  MAT_INITIAL_MATRIX, */
    /* 			PETSC_DEFAULT, &dRdUApprox); */
    printf("MAHAHHAHA\n");
    Mat mal;
    MatMatMult(reduced->Zu, dRdUApprox, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
	       &mal);
    MatView(mal, PETSC_VIEWER_STDOUT_SELF);
    printf("wanw anwanwnawnanwnaw\n");
    getchar();
  }
  //---------------------------------------------------------------------------
  // Jacobian                                                                  
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroySolverData(SolverData)); 
  destroySystem(_equation, temp->self);
  if (reduced->hrom == 1){
    VecDestroy(&residualTemp);
    MatDestroy(&dRdUTemp);
    MatDestroy(&dRdUApprox);
    MatDestroy(&ZpT);
  }
  free(temp->self);
  free(temp);
}


void yk_RunErrorEstChaos(yk_PrimalSolver *ykflow, Cluster *primalROM_h,
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
    //    yk_runReducedOrderModel(ykflow, primalROM_h, reduced);
    //-------------------------------------------------------------------------
    // Hyper-Reduced Order Modeling
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|             Hyper-Reduced-Order Modeling               |\n");
    printf("----------------------------------------------------------\n");
    yk_createHyperReducedOrderModel(ykflow, primalHROM_h, reduced);
    yk_runHyperReducedOrderModel(ykflow, primalHROM_h, reduced);
  }
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
  argvIn[1] = newFace;
  //Need to convert the xflow parameters into ykflow parameters initU
  read_input(argcIn, argvIn, xflow->ykflow->multiquation->equation, NULL, NULL,
  	     reduced);
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  sprintf(eqnReduced, "%s%s", jobFile, "Reduced");
  strcpy(xflow->ykflow->multiquation->equationReduced.nameEqn, eqnReduced);
  XFLOW->multiquation->equationReduced.numStates = reduced->nBasisFuncs;
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
  //tehcnically... thisi s wrong........
  xflow->ykflow->systemSize = (xflow->UGFine->Vector[0]->GenArray)->r*
    xflow->ykflow->primal->self->space.elem.count;
  XFLOW->index = (int *) malloc (XFLOW->systemSize *sizeof(int));
  for (i=0; i<XFLOW->systemSize; i++)
    XFLOW->index[i] = i;  
  //clsuterId
  strcpy(XFLOW->primal->clusterId, "state"); 
  strcpy(XFLOW->primal->self->utypeId, "state");
  createSystem(XFLOW->multiquation->equation, XFLOW->primal->self);
  
  //   yk_xfData2Dat(XFLOW, xflow->UGFine, jobFile);
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
  //xf_Call(xf_RunBurn(xflow->All, &AdaptInfo));
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
  printf("%s\n", AdaptInfo.SavePrefix);
  /* xf_Call(xf_ApplyTimeScheme */
  /* 	  (xflow->All, AdaptInfo.SavePrefix,xfe_False, &xflow->UGFine, */
  /* 	   xflow->TimeHistData)); */
  /* printf("i hope it works\n"); */
  /* getchar(); */
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0)); 
  //Convert the output files into ykflow reading format
  //---------------------------------------------------------------------------
  // Destroy and Free anything tht should be dealt with
  //---------------------------------------------------------------------------
  xf_DestroyAdaptInfo(&AdaptInfo, xfe_False);
  return ierr;
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
  //---------------------------------------------------------------------------
  // HROM and LSS STUFF HERE
  //---------------------------------------------------------------------------
  yk_Init(XFLOW,  reduced);
  yk_RunErrorEstChaos(XFLOW, primalROM_h, primalHROM_h, argv, reduced);
  //---------------------------------------------------------------------------
  // Destroy the Petsc and xflow stuff
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyTimeHistData(xflow->TimeHistData));
  xf_Call(xf_DestroyVectorGroup(xflow->UGFine, xfe_True, xfe_True));
  xf_CloseEqnSetLibrary(&xflow->All);
  xf_Call(xf_DestroyAll(xflow->All));
  xf_Call(xf_DestroyKeyValue(xflow->KeyValueArg));

  //---------------------------------------------------------------------------
  // Free and Unmalloc everything
  //---------------------------------------------------------------------------
  MatDestroy(&reduced->rOBState);                                             
  MatDestroy(&reduced->rOBStateBar);
  MatDestroy(&reduced->rOBResidual);
  MatDestroy(&reduced->rOBJacobian);
  MatDestroy(&reduced->Z);                                                    
  MatDestroy(&reduced->Zu);
  MatDestroy(&reduced->Zp);
  MatDestroy(&reduced->Zmult);
  destroySystem(XFLOW->multiquation->equation, XFLOW->primal->self);
  destroySystem(XFLOW->multiquation->equation, primalROM_h->self); 
  destroySystem(XFLOW->multiquation->equationReduced,  primalROM_h->reduced);  

  destroySystem(XFLOW->multiquation->equation, primalHROM_h->self);       
  destroySystem(XFLOW->multiquation->equationReduced,  primalHROM_h->reduced);
  destroySystem(XFLOW->multiquation->equation, primalHROM_h->stateFull);
  
  delArrayInt(&reduced->reducedMesh.elem);               
  delArrayInt(&reduced->reducedMesh.node); 
  SlepcFinalize(); 


  free(XFLOW->index);
  free(primalROM_h->self);
  free(primalROM_h->reduced);
  free(primalHROM_h->self);
  free(primalHROM_h->reduced);
  free(primalHROM_h->stateFull);
  free(primalROM_h);
  free(primalHROM_h);
  free(reduced);
  free(xflow->ykflow->multiquation);
  free(xflow->ykflow->primal->self);
  free(xflow->ykflow->primal);
  free(xflow->ykflow);
  free(xflow);
}
