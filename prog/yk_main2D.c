//-----------------------------------------------------------------------------
// YKFLOW: Output Based Error Estimation for Chaotic Systems                 
//                                                                           
// Yukiko Shimizu                                                            
// ykmizu@umich.edu                                                          
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// FILE: yk_main2D.c 
// Contains functions to use xflow for the primal forward solution
//-----------------------------------------------------------------------------
#include "yk_XFlow.h"

//-----------------------------------------------------------------------------
// Calls the xflow solver. ykflow solver inherits the xflow specific functions
//-----------------------------------------------------------------------------
yk_PrimalSolver* new_Xflow(){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  yk_Xflow * xflow = (yk_Xflow *) malloc (sizeof(yk_Xflow));                   
  yk_PrimalSolver* ykflow = new_Ykflow();
  ykflow->solver = xflow;                                                      
  xflow->ykflow = ykflow;
  ykflow->Function = yk_xflow_function;
  ykflow->dFunctiondu = yk_xflow_dfunctiondu;
  ykflow->Residual = yk_xflow_totalResidual;                                   
  ykflow->boundarySeeds = yk_findBoundarySeeds;
  ykflow->adjacentElems = yk_findAdjacentElems;
  ykflow->aveSpatialOutput = yk_uploadStateCalcSpatialAverageOutput;
  ykflow->spatialdJdU = yk_uploadStateCalcSpatialdJdU;
  ykflow->multiplyByInvMass = yk_MultInvMassMatrixArray;
  ykflow->injectH2h = yk_inject;
  ykflow->anyFunction = yk_anyArrayFunction;
  ykflow->delete = delete_Xflow;
  return ykflow;                                                               
}

void delete_Xflow(yk_PrimalSolver *ykflow){
  free(ykflow->solver);
  free(ykflow);
  
}
//-----------------------------------------------------------------------------
//Convert the ykflow array format to xflow Vector format
//-----------------------------------------------------------------------------
void yk_array2Vector(xf_Vector *A, Galaxy *object){                            
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j, k, n;                       // initialization for iteration
  int rOrd;                             // Num of  total basis nodes and states
  int count = 0;                        // Counter to save in correct indicies
  int numStates;                        // Number of totoal states in the sys 
  xf_GenArray *genA;                    // shortcode to group array
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  for (i=0; i<A->nArray; i++){                                                 
    genA = A->GenArray + i;                                                    
    for (j=0; j<genA->n; j++){ //Iterate through number of elements            
      rOrd = ((genA->vr == NULL) ? genA->r : genA->vr[j]);                     
      numStates = rOrd/object->basis.nodes;                                    
      if (A->GenArray[i].rValue != NULL) // Only going to look at real values
        for (k=0; k<object->basis.nodes; k++)      
          for (n=0; n<numStates; n++) //Save array values to xflow vectors
            A->GenArray[i].rValue[j][(k*numStates)+n]=                         
              object->solution[n].array[count+j*object->basis.nodes+k];
      //else if (A->GenArray[i].iValue != NULL)   //For int values
    }                                                                          
    count += genA->n*object->basis.nodes;                                      
  }                                                                            
} 

//-----------------------------------------------------------------------------
//Convert the ykflow array format to Vector Group
//-----------------------------------------------------------------------------
void yk_array2VectorGroup(xf_VectorGroup *A, Galaxy *pObj){                
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i;                                // Initialization for iteration
  for (i=0; i<A->nVector; i++)                                                 
    yk_array2Vector(A->Vector[i], pObj);                                     
}                                                                              

//-----------------------------------------------------------------------------
// Setting a Mat Petsc Object xflow matrix
//-----------------------------------------------------------------------------
void yk_matrix2D2MatGroup(xf_All *All, xf_JacobianMatrix *dRdU, Galaxy *pObj,
			  Mat matObj, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section  
  //------------------------------------//-------------------------------------
  int i, j, k, m, n, p;                 // Initialization for iteration
  int groupN, elemN, faceN;             // Actual values relative to xflow
  int negrp = dRdU->negrpself;          // Total number of groups
  int nface;                            // number of faces for each element
  int nn, nnN;
  int sr = All->EqnSet->StateRank;
  int *egrpOffset;
  int globalIndex;
  int r, rN, row, col;
  int floorElem = 0;
  int globalElem, localElem;
  int *getMeshIndex = createMeshMap(&pObj->space);
  int nBasis = pObj->basis.nodes;
  real *A;
  Mesh _meshOfInterest;                 // Mesh of interest depending of ROM
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->hrom==1)
    _meshOfInterest = reduced->reducedMesh;
  else
    _meshOfInterest = pObj->space;
  // initialize the group offset values
  xf_Error(xf_Alloc( (void **) &egrpOffset, negrp, sizeof(int)));
  //---------------------------------------------------------------------------
  // Calculate the global index offsets since there are multiple groups
  //---------------------------------------------------------------------------
  globalIndex = 0;
  for (i=0; i<negrp; i++){              // Calculate global index offsets
    egrpOffset[i] = globalIndex;
    for (j=0; j<dRdU->nElem[i]; j++){
      nn = xf_Jacobian_n(dRdU, i, j);
      globalIndex += nn*sr;
    }
  }
  //---------------------------------------------------------------------------
  // Routine to convert ykflow 2D array to Mat values
  //---------------------------------------------------------------------------
  for (j=0; j<_meshOfInterest.elem.count; j++){
    floorElem = 0;
    globalElem = _meshOfInterest.elem.array[j]; //find the global elem
    for (i=0; i<negrp; i++){            // calculate local elem in itsgroup
      floorElem += dRdU->nElem[i];      // Helps with element local value
      if (globalElem <= floorElem-1){
	localElem = globalElem;
	if (i>0)
	  for (m=0; m<i; m++)
	    localElem -= dRdU->nElem[m];
	break;
      }
    }
    nface = dRdU->nFace[i][localElem];  //gets the total number of faces
    for (k=-1; k<nface; k++){           //iterate through the number of faces
      if (k == -1){                     //Itself information
	groupN = i;
	elemN = localElem;
	faceN = k;
      }else{
	groupN = dRdU->egrpN[i][localElem][k]; //Surrounding element group
	elemN = dRdU->elemN[i][localElem][k];  //Surrounding element elem
	faceN = dRdU->faceN[i][localElem][k];  //Surrounding element face
      }
      if (groupN < 0 )
	continue;               //If boundary element skip rest
      A = dRdU->Value[i][localElem][1+k]; //Extract the block matrix values
      nn = xf_Jacobian_n(dRdU, i, localElem); //# of unknowns
      r = sr*nn;
      nnN = xf_Jacobian_n(dRdU, groupN, elemN); //# of unknowns;
      rN = sr*nnN;
      row = 0;
      col = 0;
      floorElem = 0;
      for (m=0; m<groupN; m++)
	floorElem += dRdU->nElem[m];
      elemN +=floorElem;
      for (m = 0; m<j; m++)
	row +=sr*nBasis;
      for (m = 0; m<getMeshIndex[elemN]; m++)
	col +=sr*nBasis;
      for (m=0; m<nn; m++)
	for (n=0; n<nnN; n++)
	  for (p=0; p<sr*sr; p++)
	    if (A[(m*nnN+n)*sr*sr+p]!= 0.0)
	      MatSetValue(matObj, row+m*sr+(p/sr), col+n*sr+(p%sr),
			  A[(m*nnN+n)*sr*sr+p], INSERT_VALUES);
    }
  }
  MatAssemblyBegin(matObj, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matObj, MAT_FINAL_ASSEMBLY);
  free(getMeshIndex);
  xf_Release(egrpOffset);
}

//-----------------------------------------------------------------------------
// Concert vector xflow format to Vec Petsc format                             
//-----------------------------------------------------------------------------
void yk_vector2Vec(xf_All *All, xf_Vector *A, Galaxy *pObj, Vec vecObj,
		   Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k, m;                          // initialization for iteration
  int rOrd;                                                                    
  int globalElem;
  int localElem;
  int floorElem;
  xf_GenArray *genA;                                                           
  PetscScalar value;
  Mesh _meshOfInterest;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->hrom==1)
    _meshOfInterest = reduced->reducedMesh;
  else
    _meshOfInterest = pObj->space;  
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  for (j=0; j<_meshOfInterest.elem.count; j++){
    floorElem = 0;                                                     
    globalElem = _meshOfInterest.elem.array[j]; //find the global elem 
    for (i=0; i<A->nArray; i++){            // calculate local elem in itsgroup
      genA = A->GenArray + i;
      floorElem += genA->n;      // Helps with element local value
      if (globalElem <= floorElem-1){ 
        localElem = globalElem;   
        if (i>0)             
          for (m=0; m<i; m++)                                     
            localElem -= A->GenArray[m].n;                    
        break;  
      }                                 
    }         
    rOrd = ((genA->vr == NULL) ? genA->r : genA->vr[localElem]);
    if (A->GenArray[i].rValue != NULL){                                      
      for (k=0; k<rOrd; k++){                                                
	value = A->GenArray[i].rValue[localElem][k];

	VecSetValue(vecObj, j*rOrd+k, value, INSERT_VALUES);
      }                                                                      
    } //else if (A->GenArray[i].iValue != NULL) int value 
  }                                                                            
}                                                                              


//-----------------------------------------------------------------------------
// Convert vector xflow group format to Vec Petsc format
//-----------------------------------------------------------------------------
void yk_vector2VecGroup(xf_All *All, xf_VectorGroup *A, Galaxy *pObj,
			Vec vecObj, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i;                                // Initialization for iteration
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<A->nVector; i++)
    yk_vector2Vec(All, A->Vector[i], pObj, vecObj, reduced);
}


void yk_add2VectorScalar(xf_Vector *A, double randomV){
  int i, j, n;
  for (n=0; n<A->nArray; n++)
    for (i=0; i< A->GenArray[n].n; i++)
      for (j=0; j<A->GenArray[n].r; j++)
	A->GenArray[n].rValue[i][j] +=randomV;
  
}

void yk_add2VectorGroupScalar(xf_VectorGroup *A, double randomV){
  int i;
  for (i=0; i<A->nVector; i++)
    yk_add2VectorScalar(A->Vector[i], randomV);
}

//-----------------------------------------------------------------------------
// Convert xflow text file to ykflow state file 
//-----------------------------------------------------------------------------
void yk_xfData2Dat(Galaxy *gObj, xf_VectorGroup *UG, char *jobFile){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                   
  //------------------------------------//-------------------------------------
  int i,j;                              // Initialization for iteration
  int length;
  int numStates;
  int stateCount = 0;
  int status;
  int count = 0 ;
  int numFiles = gObj->time.count;
  char stateInput[xf_MAXSTRLEN];
  char newoutFile[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char numFilesString[xf_MAXSTRLEN];
  char dataString[xf_MAXSTRLEN];
  sprintf(stateInput, "%s_U", jobFile);
  sprintf(numFilesString, "%d", numFiles);
  char * argv[] =  {"/Users/yukikoshimizu/xflow/bin/xf_Data2Text", " -inroot ",
  		    stateInput, " -batch", " 0 1 ",  numFilesString, NULL};
  char outputS[200];
  char * env[] = {NULL};
  FILE *dataFile;
  FILE *newFile;
  //pid_t pid = fork();
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  sprintf(outputS, "%s", argv[0]);
  for (i=1; i<6; i++)
    strcat(outputS, argv[i]);
  numStates = UG->Vector[0]->StateRank;
  system(outputS);
  /* if (pid == -1){ //Error */
  /*   perror("Error"); */
  /* }else if (pid > 0){ //If this is the parent process */
  /*   waitpid(pid, &status, 0); */
  /* }else{ */
  /*   execve(argv[0], argv, env); */
  /*   _exit(EXIT_FAILURE); //Needed incase execve fails */
  /* } */
  for (i=0; i<numFiles+1; i++){
    stateCount = 0;
    for (j=0; j<UG->nVector; j++){
      sprintf(outFile, "%s%d_%s.txt", stateInput, i,
	      xfe_VectorRoleName[UG->Type[j]]);
      //Should add yflow->primal_H need to add the coarse grid information
      sprintf(newoutFile, "%s_%d.dat", gObj->id, i); 
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

void yk_findAdjacentElems(yk_PrimalSolver *ykflow, Cluster *primal, int elem,
			  int *elemSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i, j, m;                             // initialization for iteration
  int ierr;
  int floorElem = 0;
  int nFace;
  int elemN;
  int newElem;
  int groupN;
  enum xfe_Bool Found; 
  xf_JacobianMatrix *dRdU;              // structure for jacobian matrix
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver; //Set pointer to xflow 
  //---------------------------------------------------------------------------
  // Find the structre of the Jacobian matrix. Will not be calculating values
  //---------------------------------------------------------------------------
  xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, xflow->UG_h,
				NULL, xfe_True, NULL, NULL, &dRdU, &Found));
  //elem: current element of interest, want to find its adjacent elements
  //Need to translate the elem into its element number relative to its group #
  for (i=0; i<dRdU->negrp; i++){
    floorElem += dRdU->nElem[i]; //Calculates the number of elements
    if (elem <= floorElem-1){
      newElem = elem;
      if (i>0)
	for (m=0; m<i; m++)
	  newElem -= dRdU->nElem[m];
      nFace = dRdU->nFace[i][newElem];
      //elem belongs to the ith element group
      elemSet[elem] = 1;
      for (j=0; j<nFace; j++){
	elemN = dRdU->elemN[i][newElem][j];
	groupN = dRdU->egrpN[i][newElem][j];
	for (m=0; m<groupN; m++)                                            
	  elemN += dRdU->nElem[m];  
	if (elemN>-1) //Don't count the boundary elements
	  elemSet[elemN] = 1;
      }
      floorElem = 0;
      break;
    }
  }
}

void yk_findBoundarySeeds(yk_PrimalSolver *ykflow, Multiverse* multiquation,
			  Cluster *primal, int nSampNodes, int *numSeeds,
			  int *nodeSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j, k;                          //initialization for iteration
  int count = 0;                        //will give me total # of elems
  int elem;                             //saves the current element to add
  int numBFG;                           //total num of boundary group 
  int numBasis = primal->self->basis.nodes;
  int elemSize = numBasis*multiquation->equation.numStates;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;                      
  xf_BFaceGroup *aBFG = xflow->All->Mesh->BFaceGroup;
  xf_ElemGroup *elemGroup = xflow->All->Mesh->ElemGroup;                       
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  numBFG = xflow->All->Mesh->nBFaceGroup; //total num of boundary group
  for (i=0; i<numBFG; i++){
    elem = 0;                         //set to 0, want global elem number
    j= aBFG[i].nBFace-1; //Just need one element on the boundary for each group
    for (k=0; k< aBFG[i].BFace[j].ElemGroup; k++) //elemgroup         
      elem += elemGroup[k].nElem;     //find global number
    elem += aBFG[i].BFace[j].Elem;    //find actual elem value
    for (k=0; k<multiquation->equation.numStates; k++){
      //rate through total elem vec size 
      //nodeSet[elemSize*count+k] = elemSize*elem+k; //node vector value
      nodeSet[multiquation->equation.numStates*count+k] =
	elemSize*elem+k;
      //nodeSet[elemSize*count] = elemSize*elem;
      //printf("%d\n", elemSize*count+k);
      // }
    }
    count++;                          //sums up tot number of elems thus far
  }
  if (nSampNodes <= count*multiquation->equation.numStates){
    //Verfiy there's enough sample nodes
    printf("Need to specify at least %d more nSampleNodes\n", count*
	   multiquation->equation.numStates- nSampNodes);        
    exit(EXIT_FAILURE);                 //quit if user needs to specify more
  }else
    *numSeeds = count*multiquation->equation.numStates;       
}
//-----------------------------------------------------------------------------
// Residual Calculation for xflow
//-----------------------------------------------------------------------------
void yk_xflow_totalResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, Vec residual, Mat Jacobian,
			    Is_it *reduced, int timeSolveNum){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int ierr;
  enum xfe_Bool Found;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;                      
  xf_JacobianMatrix *dRdU;                                                     
  xf_VectorGroup *RG = NULL;                                                   
  xf_VectorGroup *UGi[2];                                                      
  xf_SolverData *SolverData=NULL;                                              
  Galaxy *postTime = NULL;
  Vec primalReducedVec;
  Vec primalFullVec;
  Cluster *temp = (Cluster *) malloc (sizeof(Cluster));                        
  temp->self = (Galaxy *) malloc (sizeof(Galaxy));
  Universe _equation = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced; 
  //  reduced->hrom = 0;
  //need to go in here and modify it for the hyper reduced case
  //---------------------------------------------------------------------------
  // Implementation                                                            
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Set up the State temp value to be used in the Residual                    
  //-------------------------------------------------------------------------- 
  //Should probably make this UG_h part more... broad
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state1",
				     xfe_True, xfe_True, NULL, &UGi[0], NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state2",
				     xfe_True, xfe_True, NULL, &UGi[1], NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "RG",           
                                    xfe_True, xfe_True, NULL, &RG, NULL));       if (Jacobian != NULL)
    xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, xflow->UG_h,
				  NULL, xfe_True, NULL, NULL, &dRdU, &Found));
  else
    dRdU = NULL;
  //create/allocate SOlverData 
  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));                       
  SolverData->DoProcessResidual = xfe_False; // do not process residual        
  SolverData->UseArtificialDt = xfe_False;                                     
  SolverData->c = 0.;
  //---------------------------------------------------------------------------
  //Create a temporary Cluster called temp to save timenode - 1 values         
  //---------------------------------------------------------------------------
  //Left side of the element
  if (reduced->hrom == 0){
    postTime = primal->self;
  }else if (reduced->hrom == 1){
    VecCreate(PETSC_COMM_SELF, &primalReducedVec);                
    VecSetSizes(primalReducedVec, reduced->nBasisFuncs, reduced->nBasisFuncs);
    VecSetType(primalReducedVec, VECSEQ); 

    VecCreate(PETSC_COMM_SELF, &primalFullVec);
    VecSetSizes(primalFullVec, primal->stateFull->systemSize,
		primal->stateFull->systemSize);
    VecSetType(primalFullVec, VECSEQ);
    
    primal->stateFull->time.node = primal->self->time.node;
    array2Vec(_equationReduced, primal->reduced, primalReducedVec); 
    MatMult(reduced->rOBState, primalReducedVec, primalFullVec);
    vec2Array(_equation, primal->stateFull, primalFullVec);
    postTime = primal->stateFull;    
  }
  //remeber that stateFull is same size as self for hrom == 0
  ks_copyUtype(postTime, temp->self);
  copySystem(_equation, postTime, temp->self);
  temp->self->time.node = postTime->time.node - 1;
  strcpy(temp->clusterId, primal->clusterId);
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
  //Convert to Petsc Vec format (Residual)
  yk_vector2VecGroup(xflow->All, RG, primal->self, residual, reduced);
  //yk_vector2VecGroup(primalSolverObj, RG, residual, reduced);
  //Convert to Petsc Mat format (Jacobian)
  if (Jacobian != NULL)
    yk_matrix2D2MatGroup(xflow->All, dRdU, primal->self, Jacobian, reduced);
  //---------------------------------------------------------------------------
  // Jacobian                                                                  
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroySolverData(SolverData)); 
  destroySystem(_equation, temp->self);
  if (reduced->hrom == 1){
    VecDestroy(&primalReducedVec);
    VecDestroy(&primalFullVec);
  }
  free(temp->self);
  free(temp);
}

void yk_MultMatrixInvMassMatrixGroup(xf_All *All, xf_JacobianMatrix *R_U,
				     Galaxy *primal, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j, k, m;                          // initialization for iteration
  int ierr;
  real *A;
  real *C = NULL;
  enum xfe_BasisType Basis;   
  int globalElem, localElem;
  int n;
  int floorElem;
  //xf_Matrix *iM;
  real *iMM;
  real fac;
  int negrp = R_U->negrpself; //Number of groups 
  int sr = All->EqnSet->StateRank;
  Mesh _meshOfInterest;
  int elemN, faceN, groupN;
  int nface;
  int Order;
  if (reduced->hrom==1)                                                        
    _meshOfInterest = reduced->reducedMesh;                                    
  else                                                                      
    _meshOfInterest = primal->space;  
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  n =  (R_U->UG->Vector[0]->GenArray)->r/sr;
  xf_Call(xf_Alloc( (void **) &C, n*n*sr*sr, sizeof(real)));  
  for (j=0; j<_meshOfInterest.elem.count; j++){
    floorElem = 0;
    globalElem = _meshOfInterest.elem.array[j];
    for (i=0; i<negrp; i++){
      floorElem += R_U->nElem[i]; //Gives elemes
      if (globalElem <= floorElem-1){
        localElem = globalElem;
        if (i>0)
          for (m=0; m<i; m++)
            localElem -= R_U->nElem[m];
        break;
      }
    }
    xf_InterpOrderBasis(R_U->UG->Vector[0], i, localElem, &Order, &Basis);
     xf_Call(xf_ElemInvMassMatrix(All, i, localElem, Basis, Order, NULL,
				 NULL, &iMM, &fac));
    nface = R_U->nFace[i][localElem]; //gets the total number of faces
    for (k=-1; k<nface; k++){ //iterate through the number of
      if (k == -1){ //Itself inofmration
        groupN = i;
        elemN = localElem;
        faceN = k;
      }else{
        groupN = R_U->egrpN[i][localElem][k];
        elemN = R_U->elemN[i][localElem][k];
        faceN = R_U->faceN[i][localElem][k];
      }
      if (groupN>-1){
	A = R_U->Value[i][localElem][1+k];
	xf_MxBlockM(iMM, n, n, A, sr, n, xfe_Set, C);
	for (m=0; m<n*n*sr*sr; m++){
	  R_U->Value[i][localElem][1+k][m] = C[m];
	}
      }
    }
  }
  xf_Release( (void *) C);   
}

void yk_MultInvMassMatrixArray(yk_PrimalSolver *ykflow, Galaxy * galObj,
			       Vec adjoint, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int ierr;
  xf_VectorGroup * MA;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "intSol",
				    xfe_True, xfe_True, NULL, &MA, NULL));
  yk_array2VectorGroup(MA, galObj);
  xf_MultInvMassMatrixGroup(xflow->All, 1, xfe_VectorRoleElemState, MA); 
  yk_vector2VecGroup(xflow->All, MA, galObj, adjoint, reduced);
}

void yk_xflow_function(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Cluster *primal, Vec vecObj, Is_it *reduced,
		       int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  Galaxy *_solOfInterest;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)
    _solOfInterest = primal->stateFull;
  else
    _solOfInterest = primal->state;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation->equation, _solOfInterest, timeNode);
  yk_anyArrayFunction(ykflow, multiquation, _solOfInterest, vecObj, reduced);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
}

void yk_anyArrayFunction(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Galaxy *primal, Vec vecObj, Is_it *reduced){
  int ierr;
  xf_SolverData *SolverData=NULL;
  xf_VectorGroup *UG;
  xf_VectorGroup *RG = NULL; 
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver; 
  Vec nSamplefApprox;
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state1",     
				    xfe_True, xfe_True, NULL, &UG, NULL)); 
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "RG",        
                                    xfe_True, xfe_True, NULL, &RG, NULL));  
  //create/allocate SOlverData                                            
  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));                  
  SolverData->DoProcessResidual = xfe_False; // do not process residual  
  SolverData->UseArtificialDt = xfe_False;                   
  SolverData->c = 0.;
  
  yk_array2VectorGroup(UG, primal);                                 
  xf_Call(xf_CalculateResidual(xflow->All, UG, RG, NULL, SolverData));   
  xf_MultInvMassMatrixGroup(xflow->All, 1, xfe_VectorRoleElemState, RG);   
  yk_vector2VecGroup(xflow->All, RG, primal, vecObj, reduced);                
  VecScale(vecObj, -1);   
  xf_Call(xf_DestroySolverData(SolverData));
}

void yk_xflow_dfunctiondu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Cluster *primal, Mat matObj, Is_it *reduced,
			  int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section               
  //------------------------------------//-------------------------------------
  int ierr;
  enum xfe_Bool Found;  
  xf_SolverData *SolverData = NULL;
  xf_VectorGroup *UG;
  xf_VectorGroup *R;
  xf_JacobianMatrix *R_U;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  Galaxy *_solOfInterest; 
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)       
    _solOfInterest = primal->stateFull;    
  else                                                 
    _solOfInterest = primal->state;  

  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state1",
				    xfe_True, xfe_True, NULL, &UG, NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "residual",
				    xfe_True, xfe_True, NULL, &R, NULL));
  xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet,xflow->UG_h,
				NULL, xfe_True, NULL, NULL, &R_U, &Found));
   //create/allocate SOlverData                                    
  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));            
  SolverData->DoProcessResidual = xfe_False; // do not process residual 
  SolverData->UseArtificialDt = xfe_False;                                
  SolverData->c = 0.;                              
  //---------------------------------------------------------------------------
  // Implementation                                                         
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation->equation, _solOfInterest, timeNode);
  yk_array2VectorGroup(UG, _solOfInterest);
  xf_Call(xf_CalculateResidual(xflow->All, UG, R, R_U, SolverData));
  yk_MultMatrixInvMassMatrixGroup(xflow->All, R_U, _solOfInterest, reduced);
  yk_matrix2D2MatGroup(xflow->All, R_U, primal->state, matObj, reduced);
  MatScale(matObj, -1);   
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroySolverData(SolverData));  
}


void yk_calcTimeAverageOutput(yk_PrimalSolver *ykflow, Galaxy *primal,
			      xf_TimeHistData *TimeHistData){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section  
  //------------------------------------//-------------------------------------
  int i;                                // initialization for iteration
  int ierr;
  int navg = 0;                         // counter for non dimensionalizing
  int nTime;                            // number of time nodes
  int iOutput;                          // int value of the output of interest
  char AdaptOutput[xf_MAXSTRLEN];
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver; 
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  primal->j_bar = 0;               // output average of interest
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "AdaptOutput",
			 AdaptOutput));
  for (i=0; i<TimeHistData->nOutput; i++)
    if (strcmp(TimeHistData->OutputNames[i], AdaptOutput) == 0)
      break;
  if (i>=TimeHistData->nOutput)
    perror("Error: ");
  iOutput = i;
  nTime = TimeHistData->nTime;
  navg = 0;
  primal->j_bar = 0;
  for (i= 0; i<nTime; i++){
    if (i== 0 || i == nTime-1)
      primal->j_bar += TimeHistData->OutputValues[iOutput][i]*primal->time.dt/2;
    else
      primal->j_bar += TimeHistData->OutputValues[iOutput][i]*primal->time.dt;
    //navg ++;
  }
  primal->j_bar *= 1.0/(primal->time.globalT_f);
  //primal->j_bar /= ((real) navg);
}


void yk_uploadStateCalcSpatialAverageOutput(yk_PrimalSolver *ykflow,
					    Multiverse *multiquation,
					    Galaxy *primal, double *aveSpace,
					    int timeNode){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i;
  int ierr;
  xf_TimeHistData *TimeHistData;
  char LogOutput[xf_MAXSTRLEN];
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;               
  xf_VectorGroup *UGi;
  char AdaptOutput[xf_MAXSTRLEN];
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state1",    
                                    xfe_True, xfe_True, NULL, &UGi, NULL)); 
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_CreateUniformTimeHistData(xflow->All, LogOutput, &TimeHistData));
  ks_readSolution(multiquation->equation, primal, timeNode);   
  yk_array2VectorGroup(UGi, primal);
  xf_Call(xf_StoreTimeHistData(xflow->All, TimeHistData, timeNode, UGi));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "AdaptOutput",           
                         AdaptOutput));                                        
  for (i=0; i<TimeHistData->nOutput; i++)                                      
    if (strcmp(TimeHistData->OutputNames[i], AdaptOutput) == 0)                
      break;                                                                   
  if (i>=TimeHistData->nOutput)                                                
    perror("Error: ");                                                         
  *aveSpace = TimeHistData->OutputValues[i][timeNode];
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyTimeHistData(TimeHistData));  
}


void yk_uploadStateCalcSpatialdJdU(yk_PrimalSolver *ykflow,
				   Multiverse *multiquation, Galaxy *primal,
				   int timeNode, Vec vecObj, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int ierr;
  xf_VectorGroup *J_U, *UGi;
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  xf_TimeHistData *TimeHistData;
  char AdaptOutput[xf_MAXSTRLEN];
  char LogOutput[xf_MAXSTRLEN]; 
  PetscReal val;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "state1",
                                    xfe_True, xfe_True, NULL, &UGi, NULL));
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG_h, "J_U",
				    xfe_False, xfe_True, NULL, &J_U, NULL));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_CreateUniformTimeHistData(xflow->All, LogOutput, &TimeHistData));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "AdaptOutput",
                         AdaptOutput)); 
  ks_readSolution(multiquation->equation, primal, timeNode);
  yk_array2VectorGroup(UGi, primal);
  for (i=0; i<TimeHistData->nOutput; i++)
    if (strcmp(TimeHistData->OutputNames[i], AdaptOutput) == 0){
      xf_Call(xf_CalculateOutput(xflow->All, TimeHistData->OutputNames[i],
				 UGi, &TimeHistData->OutputValues[i][timeNode],
				 J_U, xfe_Set));
      break;
  }
  
  if (i>= TimeHistData->nOutput)
    perror("Error: ");
  yk_vector2VecGroup(xflow->All, J_U, primal,  vecObj, reduced);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyTimeHistData(TimeHistData));  
}
void yk_uploadStateCalcTimeAverageOutput(yk_PrimalSolver *ykflow,
					 Multiverse * multiquation,
					 xf_VectorGroup *UG_h,
					 Galaxy *primal){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;                                // initialization for iteration
  int ierr;
  xf_TimeHistData *TimeHistData;         // initialization for iteration
  char LogOutput[xf_MAXSTRLEN];
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  xf_VectorGroup *UGi;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, UG_h, "state1",
				    xfe_True, xfe_True, NULL, &UGi, NULL));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_CreateUniformTimeHistData(xflow->All, LogOutput, &TimeHistData));
  for (i=0; i<primal->time.count; i++){
    ks_readSolution(multiquation->equation, primal, i);
    yk_array2VectorGroup(UGi, primal);
    xf_Call(xf_StoreTimeHistData(xflow->All, TimeHistData, i, UGi));
  }
  yk_calcTimeAverageOutput(ykflow, primal, TimeHistData);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyTimeHistData(TimeHistData));
}

void yk_Init(yk_PrimalSolver *ykflow, Multiverse *multiquation,
	     Cluster *primal_h, Cluster *primal_H, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                 
  //------------------------------------//-------------------------------------
  int i;
  int ierr;
  int nTimeStep = 0;
  char newFace[xf_MAXSTRLEN];
  char jobFile[xf_MAXSTRLEN];
  char eqnReduced[xf_MAXSTRLEN];
  char inputFile[xf_MAXSTRLEN];
  int argcIn = 2;
  char SavePrefix_H[xf_MAXSTRLEN];
  char SavePrefix_h[xf_MAXSTRLEN];
  char *argvIn[] = {"0", "9", NULL};
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver; 
  primal_h->self=(Galaxy *) malloc (sizeof(Galaxy));                   
  primal_H->self = (Galaxy *) malloc (sizeof(Galaxy));
  ykflow->numElemCol = 4;
  //FOR REDUCED ORDER MODELING AND EVERYTHING RELATED WITH YKFLOW (WHAT A NAME)
  //---------------------------------------------------------------------------
  // Define the multiverse/equation parameters EQN
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));
  strcpy(multiquation->equation.nameEqn, jobFile);
  multiquation->equation.numStates = xflow->UG_h->Vector[0]->StateRank;
  sprintf(SavePrefix_H, "%s_Coarse", jobFile);
  sprintf(SavePrefix_h, "%s_Fine", jobFile); 
  //---------------------------------------------------------------------------
  // Initalize Parameters for hyper reduced order modeling INP
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "input", inputFile));
  sprintf(newFace, "%s.inp", inputFile);
  argvIn[1] = newFace;
  //Need to convert the xflow parameters into ykflow parameters initU
  initIs_it(reduced);
  read_input(argcIn, argvIn, multiquation->equation, NULL, NULL, reduced);
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  sprintf(eqnReduced, "%s%s", jobFile, "Reduced");
  strcpy(multiquation->equationReduced.nameEqn, eqnReduced);
  multiquation->equationReduced.numStates = reduced->nBasisFuncs;
  //---------------------------------------------------------------------------
  // Define the Cluster properties here
  //---------------------------------------------------------------------------
  initUtype(primal_h->self);
  primal_h->self->time_method = 2;
  initUtype(primal_H->self); //cehckie
  primal_H->self->time_method = 2;
  //---------------------------------------------------------------------------
  // Interpolation and quad points and basis shit
  //---------------------------------------------------------------------------
  primal_h->self->basis.p = xflow->UG_h->Vector[0]->Order[0];
  primal_H->self->basis.p = xflow->UG->Vector[0]->Order[0];
  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",       
			    &nTimeStep));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",
			       &primal_h->self->time.t_f));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time",     
                             &primal_h->self->time.t_0));  
  primal_h->self->time.globalT_f = primal_h->self->time.t_f;
  primal_h->self->time.globalT_0 = primal_h->self->time.t_0;


  primal_H->self->time.t_f = primal_h->self->time.t_f;
  primal_H->self->time.t_0 = primal_h->self->time.t_0;
  primal_H->self->time.globalT_f = primal_H->self->time.t_f; 
  primal_H->self->time.globalT_0 = primal_H->self->time.t_0;  
  
  primal_h->self->time.dt = (primal_h->self->time.t_f-primal_h->self->time.t_0)
    /nTimeStep;
  
  primal_H->self->time.dt = primal_h->self->time.dt;
  //---------------------------------------------------------------------------
  //Set the total number of elements FIne
  //---------------------------------------------------------------------------
  primal_h->self->space.elem.count = 0;
  for (i=0; i<xflow->All->Mesh->nElemGroup; i++)                               
    primal_h->self->space.elem.count += xflow->All->Mesh->ElemGroup[i].nElem;
  //Set the number of basis nodes
  primal_h->self->basis.nodes = (xflow->UG_h->Vector[0]->GenArray)->r
    /xflow->UG_h->Vector[0]->StateRank;
  //Set the total number of nodes                                              
  primal_h->self->space.node.count = primal_h->self->basis.nodes*
    primal_h->self->space.elem.count;  
  //Set the systemSize
  //tehcnically... thisi s wrong........ (Need to think about this part more)
  primal_h->self->systemSize = (xflow->UG_h->Vector[0]->GenArray)->r*
    primal_h->self->space.elem.count;
  primal_h->self->index = (int *) malloc
    (primal_h->self->systemSize*sizeof(int));
  for (i=0; i<primal_h->self->systemSize; i++)
    primal_h->self->index[i] = i;  
  //clsuterId
  strcpy(primal_h->clusterId, "state"); 
  strcpy(primal_h->self->utypeId, "state");
  createSystem(multiquation->equation, primal_h->self);
  //---------------------------------------------------------------------------
  // Set the Coarse mesh and initialize
  //---------------------------------------------------------------------------
  primal_H->self->space.elem.count = primal_h->self->space.elem.count;
  primal_H->self->basis.nodes =
    (xflow->UG_H->Vector[0]->GenArray)->r/xflow->UG_H->Vector[0]->StateRank;
  //Set the totla number of nodes
  primal_H->self->space.node.count= primal_H->self->basis.nodes
    *primal_H->self->space.elem.count;
  primal_H->self->systemSize = primal_H->self->space.node.count*                        
    multiquation->equation.numStates;  
  strcpy(primal_H->clusterId, "state");
  strcpy(primal_H->self->utypeId, "state");
  createSystem(multiquation->equation, primal_H->self);
  //---------------------------------------------------------------------------
  // Calculate the time average output and concert files
  //---------------------------------------------------------------------------
  yk_xfData2Dat(primal_H->self, xflow->UG_H, SavePrefix_H);
  yk_xfData2Dat(primal_h->self, xflow->UG_h, SavePrefix_h);


  
  yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation, xflow->UG_H, primal_H->self);
  yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation, xflow->UG_h, primal_h->self);
  reduced->actualJbar = primal_H->self->j_bar - primal_h->self->j_bar;
  /*  yk_calcTimeAverageOutput(ykflow, primal_H->self, xflow->TimeHistData_H); */
  /* printf("%g\n", primal_H->self->j_bar); */
 
  /* yk_calcTimeAverageOutput(ykflow, primal_h->self, xflow->TimeHistData_h); */
  /* printf("MARIA %g\n", primal_h->self->j_bar); */
  /* getchar(); */
}

int yk_forwardSimulation(yk_Xflow *xflow){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int ierr;
  int nTimeStep, nTimeStep0, UW;
  int OrderIncrement = 1;
  char jobFile[xf_MAXSTRLEN];
  xf_AdaptInfo AdaptInfo;
  char SavePrefix_H[xf_MAXSTRLEN];    // prefix for saving files 
  char SavePrefix_h[xf_MAXSTRLEN];
  enum xfe_Bool DoRestart;
  real Time0, EndTime0;
  char LogOutput[xf_MAXSTRLEN];
  xf_TimeHistData *TimeHistData=NULL;
  double randomV;
  //---------------------------------------------------------------------------
  // Find all the job and equation related information
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));
  xf_Call(xf_ReadAllFromJobFile(jobFile, xfe_True, &xflow->All));
  xf_Call(xf_LoadEqnSetLibrary(xflow->All->EqnSet->EqnSetLibrary)); //Equation
  xf_Call(xf_EqnSetRegister(xflow->All->EqnSet));
  //Primal State (DOUBLE CHECK THIS PART)
  xf_InitAdaptInfo(&AdaptInfo);
  xf_FillAdaptInfo(xflow->All, xflow->KeyValueArg, &AdaptInfo);
  //---------------------------------------------------------------------------
  // Burn Time to find the new initial conditions
  //---------------------------------------------------------------------------
  // Get the initial condiitons defined in the equations file
  xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart",
  			     &DoRestart));
  xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG));
  //
  srand(rdtsc());                                                          
  //randomV = 0;                                        
  randomV = (double )rand()/(double)RAND_MAX*0.2-0.1; 
  yk_add2VectorGroupScalar(xflow->UG, randomV);
  // Retrieve the values so that you can revert them back for the primal
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time", &Time0));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue,"EndTime",&EndTime0));
  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",
  			    &nTimeStep0));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue,
  			    "UnsteadyWriteInterval",&UW));
  // Set the Burn time characteristics
  nTimeStep = (int) ((AdaptInfo.BurnTime/AdaptInfo.TimeStep) + 1);
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",
  			     Time0 + AdaptInfo.BurnTime));
  xf_Call(xf_SetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",
  			    nTimeStep));
  xf_Call(xf_SetKeyValue(xflow->All->Param->KeyValue, "LogOutput", "None"));
  xf_Call(xf_CreateUniformTimeHistData(xflow->All, AdaptInfo.LogOutput,
  				       &TimeHistData));
  xf_Call(xf_SetKeyValueInt(xflow->All->Param->KeyValue,
  			    "UnsteadyWriteInterval", -1));
  // Run the actual Burn Time
  printf("BURN TIME RUNNING\n");

  xf_Call(xf_ApplyTimeScheme(xflow->All, AdaptInfo.SavePrefix,xfe_False,
  			     &xflow->UG, TimeHistData));
  // Set the values back to what they were before
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",EndTime0));
  xf_Call(xf_SetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",
  			    nTimeStep0));
  xf_Call(xf_SetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_SetKeyValueInt(xflow->All->Param->KeyValue,
  			    "UnsteadyWriteInterval", UW));
  //---------------------------------------------------------------------------
  // Copy Over new initial conditions and run the coarse primal solution
  //---------------------------------------------------------------------------
  ierr =xf_Error(xf_FindSimilarVectorGroup
  		 (xflow->All, xflow->UG, "UG_H", xfe_True, xfe_False,
  		  NULL, &xflow->UG_H, NULL));
  // Set that allocated memory space to the vectorgroup of interest
  ierr = xf_Error(xf_SetVectorGroup(xflow->UG, xfe_Set, xflow->UG_H));
  // Set the initial time
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0));
  // Create the Time History
  xf_Call(xf_CreateUniformTimeHistData
  	  (xflow->All, AdaptInfo.LogOutput, &xflow->TimeHistData_H));
  sprintf(SavePrefix_H, "%s_Coarse", AdaptInfo.SavePrefix);
  // Run the coarse primal solution
  printf("COARSE SOLUTION\n");
  xf_Call(xf_ApplyTimeScheme(xflow->All, SavePrefix_H ,xfe_False, &xflow->UG_H,
  			     xflow->TimeHistData_H));
  //---------------------------------------------------------------------------
  // Copy Over new initial conditions and run fine solution
  //---------------------------------------------------------------------------
  ierr =xf_Error(xf_FindSimilarVectorGroup
              (xflow->All, xflow->UG, "UG_h", xfe_True, xfe_False,
               NULL, &xflow->UG_h, NULL));
  //Set that allocated memory space to the vectorgroup of interest
  ierr = xf_Error(xf_SetVectorGroup(xflow->UG, xfe_Set, xflow->UG_h));
  //Set the initial time
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0));
  // Project the solution into a higher order
  xf_Call(xf_ProjectVectorGroupInPlace_OrderIncrement
  	  (xflow->All->Mesh, NULL, xflow->UG_h, OrderIncrement));
  xf_Call(xf_CreateUniformTimeHistData
          (xflow->All, AdaptInfo.LogOutput, &xflow->TimeHistData_h));
  // Run the coarse primal solution
  sprintf(SavePrefix_h, "%s_Fine", AdaptInfo.SavePrefix);
  printf("FINE SOLUTION\n");
  printf("%s\n", SavePrefix_h);
  xf_Call(xf_ApplyTimeScheme(xflow->All, SavePrefix_h , xfe_False,
  			     &xflow->UG_h, xflow->TimeHistData_h));
  //---------------------------------------------------------------------------
  // Destroy and Free anything tht should be dealt with
  //---------------------------------------------------------------------------
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0));
  xf_Call(xf_DestroyTimeHistData(TimeHistData));
  xf_DestroyAdaptInfo(&AdaptInfo, xfe_False);
  return ierr;
}

void yk_inject(yk_PrimalSolver *ykflow, Multiverse *multiquation,
	       Galaxy *primal_H, Galaxy * primal_h, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                
  //------------------------------------//-------------------------------------
  int ierr;
  int OrderIncrement = 1;  
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  Vec vecObj;
  xf_VectorGroup *UG_h;
  int hrom = reduced->hrom;
  int reducedSolution = reduced->reducedSolution;
  //---------------------------------------------------------------------------
  // Intiialization
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  VecCreate(PETSC_COMM_SELF, &vecObj);
  VecSetSizes(vecObj, primal_h->systemSize, primal_h->systemSize);
  VecSetType(vecObj, VECSEQ);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  ierr =xf_Error(xf_FindSimilarVectorGroup                                
		 (xflow->All, xflow->UG_H, "UG_h", xfe_True, xfe_False,   
		  NULL, &UG_h, NULL));
  yk_array2VectorGroup(xflow->UG_H, primal_H);
  array2Vec(multiquation->equation, primal_H, vecObj);
  ierr = xf_Error(xf_SetVectorGroup(xflow->UG_H, xfe_Set, UG_h)); 
  // this thing destroys data
  //Also i need to make sure to specify time.. that primal->self time and
  //primal_H time match up
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time",
			     primal_H->time.node)); 
  xf_Call(xf_ProjectVectorGroupInPlace_OrderIncrement
          (xflow->All->Mesh, NULL, UG_h, OrderIncrement));
  yk_vector2VecGroup(xflow->All, UG_h, primal_h, vecObj, reduced);
  vec2Array(multiquation->equation, primal_h, vecObj);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  reduced->hrom = hrom;
  reduced->reducedSolution = reducedSolution;
  xf_Call(xf_DestroyVectorGroup(UG_h, xfe_True, xfe_True));  
  VecDestroy(&vecObj);
}


void yk_estimateErrorXflow(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Cluster *primal_h,Cluster *primal_H, Cluster *primal,
			  Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                    
  //------------------------------------//-------------------------------------
  int i, j;                             //initialization for iteration
  int ierr;
  int left, right;
  int numfiles = primal->self->time.count/reduced->nTimeSegs;
  int hrom;
  int _size;
  int IsItreduced;
  int timeSolveNum = 0;
  int systemReducedSize;
  char adjointName[1000];
  double delJ_psi = 0;
  double delJ_dpsi = 0;
  double del_psi = 0;
  double del_dpsi = 0;
  double timeFinal = primal->self->time.globalT_f;
  double timeInitial = primal->self->time.globalT_0;
  double timelength = (timeFinal-timeInitial)/reduced->nTimeSegs;
  Vec residual, residualReduced;
  Vec residualMinv, residualSamples;
  Vec psi, dpsi;
  Vec psifull, psihfinal;
  Universe _eqnOfInterest;
  Universe _eqnOfInterestFull = multiquation->equation;
  Galaxy *_solOfInterest;
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;  
  xf_VectorGroup *Adjoint;
  Psi1->self = (Galaxy *) malloc (sizeof(Galaxy));
  //---------------------------------------------------------------------------
  // Malloc and initialize here
  //---------------------------------------------------------------------------
  hrom = reduced->hrom;
  IsItreduced = reduced->reducedSolution;
  if (reduced->reducedSolution == 1){
    _eqnOfInterest = multiquation->equationReduced;
    _solOfInterest = primal->reduced;
    _size = reduced->nBasisFuncs;
  }else{
    _eqnOfInterest = multiquation->equationLSS;
    _solOfInterest = primal->self;
    _size = primal_h->self->systemSize; 
  }
  systemReducedSize = _eqnOfInterest.numStates; //Size of the LSS system 
  VecCreateSeq(PETSC_COMM_SELF, _size, &residualReduced);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &residual);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &residualMinv);
  VecCreateSeq(PETSC_COMM_SELF, reduced->nSampleNodes, &residualSamples);
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &psi);
  VecCreateSeq(PETSC_COMM_SELF, systemReducedSize, &dpsi);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &psifull);
  VecCreateSeq(PETSC_COMM_SELF, primal_h->self->systemSize, &psihfinal);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  // set to 0 cause want 
  reduced->hrom = 0;
  reduced->reducedSolution = 0;
  // begin error estimation for time segments
  ierr = xf_Error(xf_FindSimilarVectorGroup
		  (xflow->All, xflow->UG_h, "Adjoint", xfe_True, xfe_False,
		   NULL, &Adjoint, NULL));
  for (i=0; i<reduced->nTimeSegs; i++){
    left = numfiles*i;
    right = numfiles*(i+1);
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
      if (j > 0){
	ks_readSolution(_eqnOfInterestFull, primal_H->self, j);
	yk_inject(ykflow, multiquation,primal_H->self,primal_h->self, reduced);
	primal_h->self->time.node = primal_H->self->time.node;
	yk_xflow_totalResidual(ykflow, multiquation, primal_h, residual, NULL,
			       reduced, timeSolveNum);
        vec2Array(_eqnOfInterestFull, primal_h->self, residual);
	yk_MultInvMassMatrixArray(ykflow, primal_h->self,residualMinv,reduced);
	if (IsItreduced == 1)
	  MatMultTranspose(reduced->rOBState, residualMinv, residualReduced);
	else
	  VecCopy(residualMinv, residualReduced);
	//---------------------------------------------------------------------
	// Adjoint Manipulation
	//---------------------------------------------------------------------
	ks_readSolution(_eqnOfInterest, Psi1->self, j);
	array2Vec(_eqnOfInterest, Psi1->self, psi);
	if (IsItreduced == 1)
	  MatMult(reduced->rOBState, psi, psifull);
	else
	  VecCopy(psi, psifull);
	//Multiply by delPsi
        vec2Array(_eqnOfInterestFull, primal_h->self, psifull); 
        yk_array2VectorGroup(Adjoint, primal_h->self);
	xf_Call(xf_ProjectVectorGroupInPlace_OrderIncrement  //project
		(xflow->All->Mesh, NULL, Adjoint, -1));
	xf_Call(xf_ProjectVectorGroupInPlace_OrderIncrement  //inject          
		(xflow->All->Mesh, NULL, Adjoint, 1));
	yk_vector2VecGroup(xflow->All, Adjoint, primal_h->self, psihfinal,
			   reduced);
        VecAXPY(psifull, -1, psihfinal);
        if (IsItreduced == 1)
	  MatMultTranspose(reduced->rOBState, psifull, dpsi);
	else
	  VecCopy(psifull, dpsi);
	//---------------------------------------------------------------------
	// Error Estimation
	//---------------------------------------------------------------------
	VecDot(psi, residualReduced, &del_psi);
	VecDot(dpsi, residualReduced, &del_dpsi);
        if (j==left || j ==right){
	  del_psi *= primal->self->time.dt/2.0; 
          del_dpsi *= primal->self->time.dt/2.0; 
	}else{
	  del_psi *= primal->self->time.dt;       
          del_dpsi *= primal->self->time.dt;
	}
	delJ_psi -= del_psi;
	delJ_dpsi -=del_dpsi;
      }
    }
    destroySystem(_eqnOfInterest, Psi1->self);
  }
  //---------------------------------------------------------------------------
  // Assigned all the error estimates to the reduced structure
  //---------------------------------------------------------------------------
  reduced->delJ_psi = delJ_psi;
  reduced->delJ_dpsi = delJ_dpsi;
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroyVectorGroup(Adjoint, xfe_True, xfe_True)); 
  VecDestroy(&residual);
  VecDestroy(&residualSamples);
  VecDestroy(&residualMinv);
  VecDestroy(&residualReduced);
  VecDestroy(&psi);
  VecDestroy(&dpsi);
  VecDestroy(&psifull);
  VecDestroy(&psihfinal);
  free(Psi1->self);
  free(Psi1);
}

void yk_RunErrorEstChaos(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Cluster *primal_h, Cluster *primal_H,
			 Cluster *primalROM_h, Cluster *primalHROM_h,
			 char *argv[], Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section               
  //------------------------------------//-------------------------------------
  //check this and init about this
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;  
  multiquation->equationLSS.numStates = multiquation->equation.numStates*
    primal_h->self->space.node.count;      
  strcpy(multiquation->equationLSS.nameEqn, "naca");  
  if (reduced->reducedLSS){
    //-------------------------------------------------------------------------
    // Reduced Order Modeling
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|                 Reduced-Order Modeling                 |\n");
    printf("----------------------------------------------------------\n");
    yk_createReducedOrderModel(ykflow, multiquation, primal_h, primalROM_h,
			       reduced);
    if (!reduced->restart)
      yk_runReducedOrderModel(ykflow, multiquation, primal_h, primalROM_h,
    			      reduced);
    yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation, xflow->UG_h,
    					primalROM_h->self);
    //-------------------------------------------------------------------------
    // Hyper-Reduced Order Modeling
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|             Hyper-Reduced-Order Modeling               |\n");
    printf("----------------------------------------------------------\n");
    yk_createHyperReducedOrderModel(ykflow, multiquation, primal_h,
				    primalHROM_h, reduced);
    if (!reduced->restart)
      yk_runHyperReducedOrderModel(ykflow, multiquation, primal_h,
				   primalHROM_h, reduced);
    yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation,xflow->UG_h,
					primalHROM_h->stateFull);  
    //-------------------------------------------------------------------------
    // reduced Least Squares Shaodwing Acted upone the Primal Solution
    //-------------------------------------------------------------------------
    printf("---------------------------------------------------------\n");
    printf("|              HROM-Least Squares Shadowing             |\n");
    printf("---------------------------------------------------------\n");
    reduced->reducedSolution = 1;
    yk_leastSquaresShadowing(ykflow, multiquation,primalHROM_h, argv, reduced);
    yk_estimateErrorXflow(ykflow, multiquation, primal_h, primal_H,
    			  primalHROM_h, reduced);
  }else if (!reduced->reducedLSS){
    printf("---------------------------------------------------------\n");
    printf("|                Least Squares Shadowing                |\n");
    printf("---------------------------------------------------------\n");    
    reduced->reducedSolution = 0;
    yk_leastSquaresShadowing(ykflow, multiquation, primal_h, argv, reduced);
    yk_estimateErrorXflow(ykflow, multiquation, primal_h, primal_H, primal_h,
			  reduced);
  }
}


int yk_finalizeProgramPrintResults(yk_PrimalSolver *ykflow,
				   Multiverse *multiquation, Cluster *primal_h,
				   Cluster *primal_H, Cluster *primalROM_h,
				   Cluster *primalHROM_h, Is_it *reduced){
  int ierr;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;  //polymorphism

  printf("---------------------------------------------------------\n");
  printf("|                        Results                        |\n");
  printf("---------------------------------------------------------\n");
  printf("Actual Error = %0.16e\n", reduced->actualJbar);
  printf("HROM/LSS Estimated Error = %0.16e\n", reduced->delJ_psi);
  printf("HROM/LSS Estimated Error dpsi = %0.16e\n", reduced->delJ_dpsi);
  
  xf_Call(xf_DestroyTimeHistData(xflow->TimeHistData_h));
  xf_Call(xf_DestroyTimeHistData(xflow->TimeHistData_H)); 
  xf_Call(xf_DestroyVectorGroup(xflow->UG_H, xfe_True, xfe_True)); 
  xf_Call(xf_DestroyVectorGroup(xflow->UG_h, xfe_True, xfe_True));
  xf_CloseEqnSetLibrary(&xflow->All);
  xf_Call(xf_DestroyAll(xflow->All));

  xf_Call(xf_DestroyKeyValue(xflow->KeyValueArg));
  //---------------------------------------------------------------------------
  // Free and Unmalloc everything
  //---------------------------------------------------------------------------
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
    delArrayInt(&reduced->reducedMesh.elem);
    delArrayInt(&reduced->reducedMesh.node);
 }
 destroySystem(multiquation->equation, primal_H->self);
 destroySystem(multiquation->equation, primal_h->self);
 free(primal_h->self->index);
 free(primal_h->self);
 free(primal_H->self);
 free(reduced);
 free(primalHROM_h);
 free(primalROM_h);
 free(primal_h);
 free(primal_H);
 free(multiquation);
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
  Multiverse *multiquation = (Multiverse *) malloc (sizeof(Multiverse));
  Cluster *primalROM_h = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primalHROM_h = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primal_h = (Cluster *) malloc (sizeof(Cluster));
  Cluster *primal_H = (Cluster *) malloc (sizeof(Cluster));
  Is_it *reduced =(Is_it*) malloc (sizeof(Is_it));
  yk_PrimalSolver *ykflow = new_Xflow();
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;  //polymorphism
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
  // Search for Attractor/Trajectory Burn Time using xflow only
  //---------------------------------------------------------------------------
  yk_forwardSimulation(xflow);
  //---------------------------------------------------------------------------
  // Initalize everything related to ykflow stuff 
  //---------------------------------------------------------------------------
  yk_Init(ykflow, multiquation, primal_h, primal_H, reduced);
  //---------------------------------------------------------------------------
  // HROM and LSS STUFF HERE
  //---------------------------------------------------------------------------
  
  yk_RunErrorEstChaos(ykflow, multiquation, primal_h, primal_H, primalROM_h,
  		      primalHROM_h, argv, reduced);
  //---------------------------------------------------------------------------
  // Destroy the Petsc and xflow stuff
  //---------------------------------------------------------------------------
  yk_finalizeProgramPrintResults(ykflow, multiquation, primal_h, primal_H,
  				 primalROM_h, primalHROM_h, reduced);
  delete_Xflow(ykflow);
  SlepcFinalize(); 
}
