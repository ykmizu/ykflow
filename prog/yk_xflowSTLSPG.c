
//THIS IS THE RIGHT ONE FEB 6

#include <Python.h>
#include "yk_XFlow.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "/gpfs1/yshimiz/libxlsxwriter/include/xlsxwriter.h"

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
  /* ykflow->Function = yk_xflow_function; */
  /* ykflow->dFunctiondu = yk_xflow_dfunctiondu; */
  ykflow->Residual = yk_xflow_totalResidual;
  ykflow->FomName = yk_fomBuffer;
  ykflow->boundarySeeds = yk_findBoundarySeeds;
  ykflow->minElements = yk_findElementOfTarget;
  ykflow->findMassCoefBDF = yk_getOffDiagonalJacobiansGroup;
  ykflow->adjacentElems = yk_findAdjacentElems;
  ykflow->residualSnapshotFile = yk_createSnapshotResidualFile;
  /* ykflow->aveSpatialOutput = yk_uploadStateCalcSpatialAverageOutput; */
  /* ykflow->spatialdJdU = yk_uploadStateCalcSpatialdJdU; */
  /* ykflow->multiplyByInvMass = yk_MultInvMassMatrixArray; */
  /* ykflow->injectH2h = yk_inject; */
  /* ykflow->anyFunction = yk_anyArrayFunction; */
  ykflow->delete = delete_Xflow;
  return ykflow;
}

void delete_Xflow(yk_PrimalSolver *ykflow){
  free(ykflow->solver);
  free(ykflow);

}

void yk_textBoldBorder(char *text){
  int i;
  int sizechar = strlen(text);
  int spaces = (80-sizechar)/2;
  int extra = (80-sizechar) % 2;
  printf("+======================================="	\
	 "=======================================+\n");
  printf("|");
  for (i=1; i<spaces; i++)
    printf(" ");
  printf("%s", text);
  for (i=0; i<spaces+extra-1; i++)
    printf(" ");
  printf("|\n");
  printf("+======================================="	\
	 "=======================================+\n");
}



void yk_textBorder(char *text){
  int i;
  int sizechar = strlen(text);
  int spaces = (80-sizechar)/2;
  int extra = (80-sizechar) % 2;
  printf("+---------------------------------------"	\
	 "---------------------------------------+\n");
  printf("|");
  for (i=1; i<spaces; i++)
    printf(" ");
  printf("%s", text);
  for (i=0; i<spaces+extra-1; i++)
    printf(" ");
  printf("|\n");
  printf("+---------------------------------------"	\
	 "---------------------------------------+\n");
}



//Also assigns the spatial elems to the xflow framework
//X's are the targets and they go into the sample mesh part of the xflow
//O's are the targets and surroundings that go in the ykflow framework
void yk_findAdjacentElems(yk_PrimalSolver *ykflow, Cluster *primal,
			  PetscInt ** mshOs, PetscInt **temporalSet, Is_it *reduced){
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
  int grpnum;
  int nn = 1208;
  int pelem;
  int egR, eR, faceR;
  int totElem =  primal->self->space.elem.count;
  int *keepTrack = (int *) malloc (totElem*sizeof(int));
  int pie;
  int count = reduced->reducedMesh.elem.count;
  int counter = 0;
  int sysSize = primal->self->systemSize;
  int numBasis = primal->self->basis.nodes;
  int elemSize = numBasis * 4;
  Mat timeIdentity;
  int sysTime = primal->self->time.count;
  int timeCount = reduced->nSampleTime-1;
  MatCreate(PETSC_COMM_SELF, &timeIdentity);
  MatSetSizes(timeIdentity, timeCount, sysTime, timeCount, sysTime);
  MatSetType(timeIdentity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(timeIdentity, 1, NULL);
  //Fill in the identity matrix based on the time indices we want to keep
  for (j=0; j<timeCount; j++)
    MatSetValue(timeIdentity, j, (*temporalSet)[j], 1, INSERT_VALUES);
  MatAssemblyBegin(timeIdentity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(timeIdentity, MAT_FINAL_ASSEMBLY);
  for (grpnum=0; grpnum<xflow->All->Mesh->nElemGroup; grpnum++)
    xflow->All->Mesh->ElemGroup[grpnum].s_nElem = 0;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  //mshOs is important... will work on this part later.
  for (i=0; i<totElem; i++)
    keepTrack[i] = 0;
  for (i=0; i<reduced->reducedMesh.elem.count; i++)
    keepTrack[reduced->reducedMesh.elem.array[i]]=1;

  /* int * sElem = xflow->All->Mesh->ElemGroup[0].sElem; */
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->reducedMesh.elem.count; i++){
    xf_Index2EgrpElem(xflow->All->Mesh, reduced->reducedMesh.elem.array[i],
		      &grpnum, &pelem);
    if (xflow->All->Mesh->ElemGroup[grpnum].s_nElem == 0)
      ierr = xf_Error(xf_Alloc((void **)
			       &xflow->All->Mesh->ElemGroup[grpnum].sElem,
    			       1, sizeof(int)));
    else
      ierr = xf_Error(xf_ReAlloc((void **)
				 &xflow->All->Mesh->ElemGroup[grpnum].sElem,
    				 xflow->All->Mesh->ElemGroup[grpnum].s_nElem+1
    				 , sizeof(int)));
    xflow->All->Mesh->ElemGroup[grpnum].
      sElem[xflow->All->Mesh->ElemGroup[grpnum].s_nElem] = pelem;
    for (j=0; j<xflow->All->Mesh->ElemGroup[grpnum].nFace[pelem];j++){
      xf_NeighborAcrossFace(xflow->All->Mesh, grpnum, pelem, j, &egR, &eR,
			    &faceR);
      if (eR!=-1){
	xf_EgrpElem2Index(xflow->All->Mesh, egR, eR, &pie);

	if (keepTrack[pie] ==0 && eR != -1){
	  keepTrack[pie] = 1;
	  count ++;
	}
      }
    }
    xflow->All->Mesh->ElemGroup[grpnum].s_nElem++;
  }
  xf_NeighborAcrossFace(xflow->All->Mesh, 0, 0, 0,
			&egR, &eR, &faceR);

  //Make the mshOs, the number of columns in the jacobian matrix
  *mshOs = (PetscInt *) malloc (count*sizeof(PetscInt));

  reduced->reducedMesh_Os.i_ElemCom = (int *) malloc (totElem*sizeof(int));
  initArrayInt(&reduced->reducedMesh_Os.elem, count);
  initArrayInt(&reduced->reducedMesh_Os.node, numBasis*count);


  for (i=0; i<totElem; i++)
    reduced->reducedMesh_Os.i_ElemCom[i] = 0;

  MatCreate(PETSC_COMM_SELF, &reduced->spaceZ);
  MatSetSizes(reduced->spaceZ, count*elemSize, sysSize, count*elemSize, sysSize);
  MatSetType(reduced->spaceZ, MATSEQAIJ);
  MatSeqAIJSetPreallocation(reduced->spaceZ, 1, NULL);


  for (i=0 ; i<totElem; i++){
    if (keepTrack[i] > 0){


      reduced->reducedMesh_Os.elem.array[counter]=i;
      reduced->reducedMesh_Os.i_ElemCom[i] = counter;
      for (j=0; j<numBasis; j++){
	reduced->reducedMesh_Os.node.array[counter*numBasis+j]=i*numBasis+j;
      }
      for (j=0; j<elemSize; j++){
	MatSetValue(reduced->spaceZ, counter*elemSize+j, i*elemSize+j, 1, INSERT_VALUES);
      }

    (*mshOs)[counter] = i;
    counter++;
    }
  }


  MatAssemblyBegin(reduced->spaceZ, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(reduced->spaceZ, MAT_FINAL_ASSEMBLY);

  reduced->reducedMesh_Os.systemSize = reduced->reducedMesh_Os.elem.count*12;

  /* PetscInt m, n; */
  yk_kkron(timeIdentity, reduced->spaceZ, &reduced->Z_Os);
  free(keepTrack);
  MatDestroy(&timeIdentity);
}

void yk_findElementOfTarget(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, PetscInt *spatialSet,
			    Is_it* reduced){

  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, k;                          //initialization for iteration
  int count = 0;                        //will give me total # of elems
  int elem;                             //saves the current element to add
  int numBFG;                           //total num of boundary group
  int numBasis = multiquation->equation.numStates*3;
  int elemSize = numBasis*multiquation->equation.numStates;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  xf_BFaceGroup *aBFG = xflow->All->Mesh->BFaceGroup;
  xf_ElemGroup *elemGroup = xflow->All->Mesh->ElemGroup;
  int nSpatialDOF = reduced->nSampleNodes;
  xf_Mesh *entireMesh = xflow->All->Mesh;
  PetscInt numElems = primal->self->space.elem.count;
  PetscInt numNodes = primal->self->space.node.count;
  PetscInt *elemTrack = (PetscInt *) malloc (numNodes*sizeof(PetscInt));
  PetscInt eTarget;
  for (i=0; i<numNodes; i++)
    elemTrack[i] = 0;
  printf("%d\n", entireMesh->nElemGroup);
  printf("Does it work\n");

  xf_ConvTable *Node2Elem;
  printf("Number of nodes in Mesh %d\n", entireMesh->nNode);
  printf("what waht waht %d\n", entireMesh->ElemGroup[0].Node[0][0]);
  printf("what waht waht %d\n", entireMesh->ElemGroup[0].Node[0][1]);
  printf("what waht waht %d\n", entireMesh->ElemGroup[0].Node[0][2]);
  xf_Node2ElemConnectivity(entireMesh, &Node2Elem);
  printf("%d\n", Node2Elem[0].nItem);
  printf("%d\n", Node2Elem[0].Item[0]);
  printf("%d\n", Node2Elem[294].nItem);
  printf("%d\n", Node2Elem[294].Item[0]);
  printf("%d\n", Node2Elem[294].nItem);
  printf("%d\n", Node2Elem[294].Item[1]);
  printf("%d\n", Node2Elem[2171].nItem);
  printf("%d\n", Node2Elem[2171].Item[0]);
  printf("%d\n", Node2Elem[1298].nItem);
  printf("%d\n", Node2Elem[1298].Item[0]);
  printf("%d\n", Node2Elem[1299].nItem);
  printf("%d\n", Node2Elem[1299].Item[0]);

  printf("did it work? i wonder\n");
  printf("num basis? %d \n", numBasis);
  printf("factor factor %d\n", numNodes);

  for(i=0; i<nSpatialDOF; i++)
    printf("%f\n", floor(spatialSet[i]/numBasis));
  for(i=0; i<nSpatialDOF; i++){
    eTarget = spatialSet[i]/numBasis;
    printf("%d\n", elemTrack[eTarget]);
    if (elemTrack[eTarget] == 0){
      elemTrack[eTarget] = 1;
      /* for (j=0; j<Node2Elem[eTarget].nItem; j++) */
      /* 	printf("%d\n", Node2Elem[eTarget].Item[j]); */
      printf("%d\n", eTarget);
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
  //This may need to be changed later............
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
  }else{
    *numSeeds = count*multiquation->equation.numStates;
  }
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
  int r, rN, row, col, actualcol;
  int floorElem = 0;
  int globalElem, localElem;
  int *getMeshIndex = createMeshMap(&pObj->space);
  double value;
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
      if (reduced->hrom==1)
	elemN = reduced->reducedMesh_Os.i_ElemCom[elemN];
      for (m = 0; m<j; m++)
        row +=sr*nBasis;
      for (m = 0; m<getMeshIndex[elemN]; m++)
        col +=sr*nBasis;
      for (m=0; m<nn; m++)
        for (n=0; n<nnN; n++)
          for (p=0; p<sr*sr; p++)
            if (A[(m*nnN+n)*sr*sr+p]!= 0.0){
	      if(!isnan(A[(m*nnN+n)*sr*sr+p])) {
		MatSetValue(matObj, row+m*sr+(p/sr), col+n*sr+(p%sr),
			    A[(m*nnN+n)*sr*sr+p], INSERT_VALUES);
	      }else{
		value = 0.0;
		MatSetValue(matObj, row+m*sr+(p/sr), col+n*sr+(p%sr),
                            value, INSERT_VALUES);
	      }
	    }

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
	if (!isnan(value)) {
	  VecSetValue(vecObj, j*rOrd+k, value, INSERT_VALUES);
	}else{
	  value = 0.0;
	  VecSetValue(vecObj, j*rOrd+k, value, INSERT_VALUES);
	}
      } //else if (A->GenArray[i].iValue != NULL) int value
    }
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

void yk_fomBuffer(yk_PrimalSolver *ykflow, Universe equation, Is_it *reduced, \
		  int p, char *nameOfDir){

  int mach= reduced->paramMeshGrid[p*reduced->numParams]*10000;
  int alfa = round(reduced->paramMeshGrid[p*reduced->numParams+1]*1000);
  int reynolds = reduced->paramMeshGrid[p*reduced->numParams+2];
  sprintf(nameOfDir, "%s/%s_M_%d_A_%d_Re_%d", ykflow->path, equation.nameEqn,
	  mach,  alfa, reynolds);
}


  /* xf_Call(xf_GetTimeStep(TimeHistData, iTime, &TimeScheme, &Time, &TimeStep));   */
  /* //Get the Param_information                                                    */
  /* xf_Call(xf_GetParam_MultiStep(TimeScheme, &nStep, &nStepR, &a, &b));      */

/* int yk_getOffDiagonalJacobians(yk_PrimalSolver *ykflow, Galaxy *primal, */
/* 				Mat *Jacobian, */
/* 				xf_Vector *U, xf_Vector *R, */
/* 				xf_JacobianMatrix *R_U, */
/* 				Is_it *reduced){ */
 /*  //--------------------------------------------------------------------------- */
 /*  // */
 /*  //--------------------------------------------------------------------------- */
 /*  int ierr, egrp, elem; */
 /*  int Order, pOrder, pOrderU, pOrderR, OrderU, OrderR; */
 /*  int i, j, ij, k, sr, sr2, nnU, nnR; */
 /*  real *MM, fac, val, J, *A; */
 /*  xf_Matrix *M; */
 /*  xf_Mesh *Mesh; */
 /*  xf_JacobianMatrix *Mass=NULL; */
 /*  enum xfe_BasisType BasisU, BasisR; */
 /*  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver; */
 /*  xf_All *All = xflow->All; */
 /*  int nStep, nStepR; */
 /*  int iTime = primal->time.node-1; */
 /*  real Time, TimeStep; */
 /*  real *a = NULL; */
 /*  real *b = NULL; */

 /*  xf_TimeHistData *TimeHistData = xflow->TimeHistData; */

 /*  enum xfe_TimeSchemeType TimeScheme; */

 /*  ierr = xf_Error(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, */
 /*                                        xflow->UG, NULL, */
 /*                                        xfe_True, NULL, NULL, &Mass, NULL)); */
 /*  xf_Call(xf_SetZeroJacobian(Mass)); */
 /*  //--------------------------------------------------------------------------- */
 /*  // Implementation */
 /*  //--------------------------------------------------------------------------- */

 /*  Mesh = All->Mesh; */
 /*  sr = U->StateRank; */
 /*  sr2 = sr*sr; */

 /*  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){ */

 /*    BasisU = U->Basis[egrp]; */
 /*    BasisR = ((R==NULL) ? BasisU : R->Basis[egrp]); */
 /*    pOrderU = -1; */
 /*    pOrderR = -1; */
 /*    // loop over elements */
 /*    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){ */
 /*      //get penalization factor */
 /*      // get interpolation order */
 /*      OrderU = xf_InterpOrder(U, egrp, elem); */
 /*      OrderR = ((R==NULL) ? OrderU : xf_InterpOrder(R, egrp, elem)); */

 /*      if ((pOrderU != OrderU) || (pOrderR != OrderR)){ */
 /*        pOrderU = OrderU; */
 /*        pOrderR = OrderR; */

 /*        // find general mass matrix data (if generic, n==1) */
 /*        ierr = xf_Error(xf_FindGenMassMatrixData(All, egrp, BasisR, OrderR, */
 /* 						 BasisU, OrderU, &M)); */
 /*        /\* if (ierr != xf_OK) return ierr; *\/ */
 /*        // nnR (row order) */
 /*        ierr = xf_Error(xf_Order2nNode(BasisR, OrderR, &nnR)); */
 /*        /\* if (ierr != xf_OK) return ierr; *\/ */

 /*        // nnU (column order) */
 /*        ierr = xf_Error(xf_Order2nNode(BasisU, OrderU, &nnU)); */
 /*        /\* if (ierr != xf_OK) return ierr; *\/ */
 /* } */

 /*      ierr = xf_Error(xf_ElemGenMassMatrix(All, egrp, elem, BasisR, OrderR, */
 /* 					   BasisU, OrderU,M,NULL, &MM, &fac)); */
 /*      /\* if (ierr != xf_OK) return ierr; *\/ */
 /*      // R_U += M*(c + 1/dt) */
 /*      if (R_U != NULL){ */
 /* 	A = Mass->Value[egrp][elem][0]; */
 /* 	for (i=0; i<nnR*nnU; i++){ */
 /*          val = MM[i]; */
 /*          for (k=0; k<sr2; k+=(sr+1)) */
 /*            A[i*sr2+k] = val; */
 /*        } // i */
 /*      } */
 /*    } // elem */
 /*  } // egrp */
 /*  //Loop through the number of time steps for the time scheme */

/*   return ierr; */
/* } */
int yk_getOffDiagonalJacobiansGroup(yk_PrimalSolver *ykflow, Galaxy *primal,
				    Mat *Jacobian, Is_it *reduced){
  int ierr;
  enum xfe_Bool Found;
  xf_VectorGroup *RG;
  xf_JacobianMatrix *Mass;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  xf_SolverData *SolverData = NULL;
  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));
  SolverData->DoProcessResidual = xfe_False; // do not process residual
  SolverData->UseArtificialDt = xfe_False;
  SolverData->c = 0;
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG, "residual",
				    xfe_True,
				    xfe_True, NULL, &RG, &Found));
  ierr = xf_Error(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet,
                                        xflow->UG, NULL,
                                        xfe_True, NULL, NULL, &Mass, NULL));
  xf_SetZeroJacobian(Mass);
  xf_Call(xf_AddMassMatrixGroup(xflow->All, 1.0, NULL, xflow->UG, RG, Mass,
				SolverData));
  /* if (*Jacobian[i] != NULL) */
  yk_matrix2D2MatGroup(xflow->All, Mass, primal, *Jacobian, reduced);
  xf_Call(xf_DestroySolverData(SolverData));

  return ierr;
}

/* int yk_getOffDiagonalJacobiansGroup(yk_PrimalSolver *ykflow, Galaxy *primal, */
/* 				    Mat *Jacobian, Is_it *reduced){ */
/*   int ierr; */
/*   enum xfe_VectorRoleType VectorRole; */
/*   enum xfe_Bool Found; */
/*   yk_Xflow *xflow = (yk_Xflow *) ykflow->solver; */
/*   xf_Vector *U, *R; */
/*   xf_JacobianMatrix *R_U = NULL; */
/*   xf_VectorGroup *RG = NULL; */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   // for now, only support mass matrices associated with elements */
/*   VectorRole = xfe_VectorRoleElemState; */

/*   xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG, "residual", */
/* 				    xfe_True, */
/* 				    xfe_True, NULL, &RG, &Found)); */

/*   // pull off desired vectors */
/*   ierr = xf_Error(xf_GetVectorFromGroup(xflow->UG, VectorRole, &U)); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_GetVectorFromGroup(RG, VectorRole, &R)); */

/*   ierr = xf_Error(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, */
/* 					xflow->UG, NULL, */
/* 					xfe_True, NULL, NULL, &R_U, NULL)); */

/*   // call appropriate functions */
/*   switch (VectorRole){ */
/*   case xfe_VectorRoleElemState: */
/*     ierr = yk_getOffDiagonalJacobians(ykflow, primal, Jacobian, U, R, R_U, */
/* 				      reduced); */
/*     if (ierr != xf_OK) return ierr; */
/*     break; */
/*   default: */
/*     return xf_Error(xf_NOT_SUPPORTED); */
/*     break; */
/*   } */
/* } */


void yk_xflow_totalResidual(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Galaxy *state, Vec residual, Mat Jacobian,
			  Is_it *reduced, int timeSolveNum){
  int i, j, k;
  int ierr;
  char jobFile[xf_MAXSTRLEN];
  char title[xf_MAXSTRLEN];
  int nStep, nStepR;

  real Time, TimeStep;
  real *a = NULL;
  real *b = NULL;
  enum xfe_Bool Found, LeanFlag=0;
  xf_VectorGroup **UGi=NULL;
  xf_VectorGroup *SG=NULL;
  xf_VectorGroup *RG = NULL;
  enum xfe_TimeSchemeType TimeScheme;
  xf_JacobianMatrix *R_U;
  int initalNode = round(state->time.t_0/state->time.dt);
  int iTime = state->time.node-initalNode-1;
  double resTime = state->time.node*state->time.dt;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  xf_TimeHistData *TimeHistData = xflow->TimeHistData;
  xf_SolverData *SolverData=NULL;
  xf_DataSet *DataSet = NULL;
  char SavePrefix[xf_MAXSTRLEN];
  xf_Data *D=NULL;
  int ttime;
  char LogOutput[xf_MAXSTRLEN];
  double output;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput",LogOutput));

  /* xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile)); */
  /* sprintf(SavePrefix, "%s", jobFile); */

  xf_Call(xf_CreateSolverData(xflow->All, &SolverData));
  SolverData->DoProcessResidual = xfe_False; // do not process residual
  SolverData->UseArtificialDt = xfe_False;
  /* //Get the time information */
  xf_Call(xf_GetTimeStep(TimeHistData, iTime, &TimeScheme, &Time, &TimeStep));
  /* //Get the Param_information */
  xf_Call(xf_GetParam_MultiStep(TimeScheme, &nStep, &nStepR, &a, &b));

  //Allocate the state for each time that is needed for the discretization

  xf_Call(xf_Alloc( (void **) &UGi, nStep+1, sizeof(xf_VectorGroup *)));
  for (i=0; i<=nStep; i++){
    sprintf(title, "state_%d", i);
    xf_Call(xf_FindSimilarVectorGroup(xflow->All, xflow->UG, title, xfe_True,
  				      xfe_True, NULL, UGi+i, &Found));
    xf_Call(xf_SetZeroVectorGroup(UGi[i])); // zero out if created
  }
  /* //Source vector create */
  sprintf(title, "Residual_%d", 0);
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, UGi[0], title, xfe_True,
				    xfe_True, NULL, &SG, &Found));
  /* //M/dt*(a1*UG^n + a2*UG^{n-1} + ...) + b1*RG^{n} + ... */
  /* //State portion */
  /* //THE GUESS */
  //---------------------------------------------------------------------------
  // Place the guess for the i+1 values (To be used in the residual portion)
  //---------------------------------------------------------------------------
  yk_array2VectorGroup(UGi[0], state);

  //---------------------------------------------------------------------------
  // Create the vector sets for the Ui+1 and Ui parks the non sptial R part
  //---------------------------------------------------------------------------
  xf_Call(xf_CreateDataSet(&DataSet));
  ttime = state->time.node;
  for (i=1; i<=nStep; i++){
    /* printf("%d\n", ttime); */
    /* printf("gala %d\n", i); */
    /* printf("readsolution %d\n", ttime-i); */
    ks_readSolution(multiquation->equation, state, ttime-i);
    yk_array2VectorGroup(UGi[i], state);
    if (i==1){
      xf_Call(xf_VectorGroupMultSet(UGi[i], a[i], xfe_Set, SG));
    }else{
      xf_Call(xf_VectorGroupMultSet(UGi[i], a[i], xfe_Add, SG));
    }
  }
  //---------------------------------------------------------------------------
  // Calculate the Mass Matrix
  //---------------------------------------------------------------------------
/* // now set SG = M/dt * SGhgyt7 */
  xf_Call(xf_MultMassMatrixGroup(xflow->All, 1.0/TimeStep,
				 xfe_VectorRoleElemState,SG));
  // zero out all but ElemState role in SG
  xf_Call(xf_MaskZeroVectorGroup(SG, xfe_VectorRoleElemState));
  SolverData->c = a[0]/TimeStep; // coefficient on M*U^{n+1}
  SolverData->SG = SG; // store source vector group
  //Set the values of the state: meed to edit it

  //Create Jacobian Vector
  xf_Call(xf_FindJacobianMatrix(xflow->All, xflow->All->DataSet, UGi[0], NULL,
  				!LeanFlag, NULL, NULL, &R_U, &Found));

  // locate Residual vector
  /* printf("time to calculate residual at %g\n", Time+primal->self->time.dt); */
  xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time+state->time.dt));

  //---------------------------------------------------------------------------
  // Residual calculuate
  //---------------------------------------------------------------------------
  xf_Call(xf_FindSimilarVectorGroup(xflow->All, UGi[0], "Residual", xfe_False,
  				    xfe_True, NULL, &RG, NULL));
  if (reduced->hrom==1)
    ierr = xf_Error(yk_CalculateResidual(xflow->All, UGi[0], RG, R_U,
  					 SolverData));
  else
    ierr = xf_Error(xf_CalculateResidual(xflow->All, UGi[0], RG, R_U,
  					 SolverData));
  enum xfe_VectorRoleType RnormVector;
  RnormVector = xfe_VectorRoleElemState;
  xf_Call(xf_VectorGroupNorm(RG, 1, RnormVector, &SolverData->RnormL1));
  //---------------------------------------------------------------------------
  // Convert to Petsc Vec format (Residual)
  //---------------------------------------------------------------------------
  yk_vector2VecGroup(xflow->All, RG, state, residual, reduced);
  //yk_vector2VecGroup(primalSolverObj, RG, residual, reduced);
  //Convert to Petsc Mat format (Jacobian)
  if (Jacobian != NULL)
    yk_matrix2D2MatGroup(xflow->All, R_U, state, Jacobian, reduced);


  if (reduced->JFlag == 1){
    for (i=0; i<2; i++)
      xf_CalculateOutput(xflow->All, xflow->All->EqnSet->Outputs->Output[i].Name, UGi[0],
			 &state->J[i], NULL, xfe_Set);
  }
  //---------------------------------------------------------------------------
  // Destory everythiung
  //---------------------------------------------------------------------------
  xf_Call(xf_DestroySolverData(SolverData));
  xf_Call(xf_DestroyDataSet(DataSet));
  xf_Release((void *) UGi);
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
  int leftNode = round(gObj->time.t_0/gObj->time.dt)+1;
  if (leftNode == 1){
    leftNode = 0;
  }

  int rightNode = gObj->time.t_f/gObj->time.dt;
  //CHNAGE THIS CHANGE THIS THAT 0 in the STAT EINPUT ARGV THINKG
  char files[xf_MAXSTRLEN];
  sprintf(files, " %d 1 ", leftNode);
  sprintf(stateInput, "%s_U", jobFile);
  sprintf(numFilesString, "%d", rightNode);
  char * argv[] =  {"/gpfs1/yshimiz/xflow/bin/xf_Data2Text", " -inroot ",
                    stateInput, " -batch", files,  numFilesString, NULL};

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

void definelocalMeshGrid(int paramNumber, Is_it *reduced, double *testparam,
                         int *count);

void yk_Init(yk_PrimalSolver *ykflow, Multiverse *multiquation, Cluster *fom,
	     Is_it *reduced){
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
  char OutputFile[xf_MAXSTRLEN];
  int argcIn = 2;
  char SavePrefix[xf_MAXSTRLEN];
  char *argvIn[] = {"0", "9", NULL};
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  fom->self=(Galaxy *) malloc (sizeof(Galaxy));
  ykflow->numElemCol = 4;
  int numParamSet= 1;
  double value;
  int inter;
  int count=0;
  double *testparam;

  char createParam[xf_MAXSTRLEN];

  //FOR REDUCED ORDER MODELING AND EVERYTHING RELATED WITH YKFLOW (WHAT A NAME)
  //---------------------------------------------------------------------------
  // Define the multiverse/equation parameters EQN/Define jobe file name
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));
  strcpy(multiquation->equation.nameEqn, jobFile);
  multiquation->equation.numStates = xflow->UG->Vector[0]->StateRank;
  sprintf(SavePrefix, "%s", jobFile);
  //---------------------------------------------------------------------------
  // Initalize Parameters for hyper reduced order modeling INP
  //---------------------------------------------------------------------------
  //Not sure if I will use an input file or not
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "input", inputFile));
  sprintf(newFace, "%s.inp", inputFile);
  argvIn[1] = newFace;
  //Need to convert the xflow parameters into ykflow parameters initU
  initIs_it(reduced);
  read_input(argcIn, argvIn, multiquation->equation, NULL, NULL, reduced);

  sprintf(createParam,"python3 createParam.py %g %g %g",
	  reduced->params[0], reduced->params[1], reduced->params[2]);

  system(createParam);
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  sprintf(eqnReduced, "%s%s", jobFile, "Reduced");

  strcpy(multiquation->equationReduced.nameEqn, eqnReduced);
  /* multiquation->equationReduced.numStates = reduced->nBasisFuncs; */
  //---------------------------------------------------------------------------
  // Define the Cluster properties here
  //---------------------------------------------------------------------------
  initUtype(fom->self);
  fom->self->time_method = 2;
  //---------------------------------------------------------------------------
  // Interpolation and quad points and basis shit
  //---------------------------------------------------------------------------
  fom->self->basis.p = xflow->UG->Vector[0]->Order[0];
  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",
                            &nTimeStep));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",
                               &fom->self->time.t_f));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time",
                             &fom->self->time.t_0));
  fom->self->time.globalT_f = fom->self->time.t_f;
  fom->self->time.globalT_0 = fom->self->time.t_0;

  fom->self->time.dt = (fom->self->time.t_f-fom->self->time.t_0)/nTimeStep;
  fom->self->time.count = nTimeStep;
  //---------------------------------------------------------------------------
  //Set the total number of elements FIne
  //---------------------------------------------------------------------------
  fom->self->space.elem.count = 0;
  for (i=0; i<xflow->All->Mesh->nElemGroup; i++){
    fom->self->space.elem.count += xflow->All->Mesh->ElemGroup[i].nElem;
  }//Set the number of basis nodes
  fom->self->basis.nodes = (xflow->UG->Vector[0]->GenArray)->r
    /xflow->UG->Vector[0]->StateRank;
  //Set the total number of nodes
  fom->self->space.node.count = fom->self->basis.nodes*
    fom->self->space.elem.count;
  /* //Set the systemSize */
  /* //tehcnically... thisi s wrong........ (Need to think about this part more) */
  fom->self->systemSize = (xflow->UG->Vector[0]->GenArray)->r*
    fom->self->space.elem.count;
  fom->self->index = (PetscInt *) malloc(fom->self->systemSize*
					 sizeof(PetscInt));

  reduced->nSampleNodes = (int)(reduced->pSampleElems*fom->self->space.elem.count);
  for (i=0; i<fom->self->systemSize; i++)
    fom->self->index[i] = i;
  strcpy(fom->clusterId, "state");
  strcpy(fom->self->utypeId, "state");
  createSystem(multiquation->equation, fom->self);
  //Number of outputs interested in
  fom->self->J = (double *) malloc (2*sizeof(double));
  //---------------------------------------------------------------------------
  // Parmas create the meshgrid signifying the different data
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  // Calculate the time average output and concert files
  //---------------------------------------------------------------------------
  /* printf("%d\n", nTimeStep); */
  /* printf("%d underwear \n", fom->self->time.count); */
  for (i=0; i<fom->self->time.count+1; i++){
    sprintf(OutputFile, "%s_U%d.xfa", SavePrefix, i);
    xf_Call(xf_WriteAllBinary(xflow->All, OutputFile));
  }
  //TURN THIS BACK ON LATER
  yk_xfData2Dat(fom->self, xflow->UG, SavePrefix);



  /* yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation, xflow->UG_H, primal_H->self); */
  /* yk_uploadStateCalcTimeAverageOutput(ykflow, multiquation, xflow->UG_h, primal_h->self); */
  /* reduced->actualJbar = primal_H->self->j_bar - primal_h->self->j_bar; */



  //Calculates the total number of parameter sets
  testparam = (double *) malloc (reduced->numParams*sizeof(double));

  /* for (i=0; i<reduced->numParams; i++){ */
  /*   inter= ((reduced->paramsH[i]*10-reduced->paramsL[i]*10)/ */
  /* 	    reduced->dparams[i])/10+1; */
  /*   numParamSet *= inter; */
  /* } */

  reduced->paramMeshGrid =
    (double *) malloc (reduced->numParams*reduced->numParamSet*sizeof(double));
  definelocalMeshGrid(0, reduced, testparam,  &count);
  free(testparam);

}

void definelocalMeshGrid(int paramNumber, Is_it *reduced, double *testparam,
			 int *count){
  int i;

  double value =reduced->paramsL[paramNumber];
  do{
    testparam[paramNumber] = value;
    if (paramNumber==reduced->numParams-1){
      for(i=0; i<reduced->numParams; i++)
      	reduced->paramMeshGrid[(*count)*reduced->numParams+i] =testparam[i];
      (*count)++;
    }
    if (paramNumber <reduced->numParams-1)
      definelocalMeshGrid(paramNumber+1, reduced, testparam, count);
    //return nextvalue;

    value += reduced->dparams[paramNumber];
  }while ((int)round(value*100) <(int)(reduced->paramsH[paramNumber]*100+reduced->dparams[paramNumber]*100));
}


  //Calculate the meshgrid for each parameter set


/* void createSampleMesh(yk_PrimalSolver *ykflow, xf_All *sampleAll){ */
/*   int ierr; */
/*   char joyk_bFile[xf_MAXSTRLEN]; */
/*   enum xfe_Bool DoRestart=0; */
/*   yk_Xflow *xflow = (yk_Xflow *) ykflow->solver; */
/*   //--------------------------------------------------------------------------- */
/*   // Find all the job and equation related information */
/*   //--------------------------------------------------------------------------- */
/*   xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile)); */
/*   xf_Call(xf_ReadAllFromJobFile(jobFile, xfe_True, &sampleAll)); */

/*   xf_Mesh *Mesh = sampleAll->Mesh; */
/*   int nnode = 3; */
/*   int k; */
/*   int nElem = 1; */
/*    xf_SolverData *SolverData=NULL; */

/*   xf_VectorGroup *sampleU; */
/*   printf("BEFORE IT HAPPENS\n"); */
/*   printf("SAMPLE THING\n"); */
/*   printf("airfoil\n",  xflow->All->Mesh->nBFaceGroup); */
/*   printf("are you a puppy %d\n", sampleAll->Mesh->nElemGroup); */
/*   printf("group %d\n", sampleAll->Mesh->ElemGroup[0].nElem); */
/*   printf("nodes %d\n", sampleAll->Mesh->ElemGroup[0].nNode); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][0]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][1]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][2]); */

/*    printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[1][0]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[1][1]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[1][2]); */

/*   printf("group 2 %d\n", sampleAll->Mesh->ElemGroup[1].nElem); */
/*   printf("nodes 2 %d\n", sampleAll->Mesh->ElemGroup[1].nNode); */

/*   printf("hrm %d\n", sampleAll->Mesh->nIFace); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].nFace[0]); */
/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][0].Group); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][1].Group); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][2].Group); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][0].Number); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][1].Number); */

/*   printf("dsaf %d\n", sampleAll->Mesh->ElemGroup[0].Face[0][2].Number); */
/*   printf("AFTER IT HAPPENS\n"); */
/*   sampleAll->Mesh->nBFaceGroup=0; */
/*   sampleAll->Mesh->nElemGroup = 1; */
/*   sampleAll->Mesh->ElemGroup[0].nElem = 1; */
/*   ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup[0].nFace, nElem, siz\ */
/* eof(int))); */
/*   if (ierr!=xf_OK) return ierr; */
/*   /\* for (elem=0; elem<nElem; elem++) *\/ */
/*   /\*   Mesh->ElemGroup[0].nFace[elem] = nface; *\/ */
/*   Mesh->ElemGroup[0].nNode = nnode; */
/*   /\* ierr = xf_Error(xf_ReVAlloc2((void ***) &Mesh->ElemGroup[0].Face, nElem, *\/ */
/*   /\* 			     Mesh->ElemGroup[0].nFace, sizeof(xf_Face))); *\/ */
/*   if (ierr!=xf_OK) return ierr; */
/*   ierr = xf_Error(xf_ReAlloc2((void ***) &Mesh->ElemGroup[0].Node, nElem, 3, sizeof(int))); */
/*   int elem = 0; */
/*   for (k = 0; k<3; k++) */
/*     Mesh->ElemGroup[0].Node[elem][k] = xflow->All->Mesh->ElemGroup[0].Node[0][k]; */
/*   sampleAll->Mesh->nIFace = 3; */
/*   printf("are you a puppy %d\n", sampleAll->Mesh->nElemGroup); */
/*   printf("group %d\n", sampleAll->Mesh->ElemGroup[0].nElem); */
/*   printf("nodes %d\n", sampleAll->Mesh->ElemGroup[0].nNode); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][0]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][1]); */
/*   printf("gloabal node number %d\n", sampleAll->Mesh->ElemGroup[0].Node[0][2]); */

/*   printf("SIGH %d\n", sampleAll->Mesh->IFace[230].ElemL); */
/*   printf("SIGH %d\n", sampleAll->Mesh->IFace[230].ElemR); */

/*   /\* printf("group 2 %d\n", sampleAll->Mesh->ElemGroup[1].nElem); *\/ */
/*   /\* printf("nodes 2 %d\n", sampleAll->Mesh->ElemGroup[1].nNode); *\/ */
/*     xf_Call(xf_LoadEqnSetLibrary(sampleAll->EqnSet->EqnSetLibrary)); //Equation */
/*   xf_Call(xf_EqnSetRegister(sampleAll->EqnSet)); */
/*   xf_TimeHistData *TimeHistData=NULL; */
/*   char LogOutput[xf_MAXSTRLEN]; */
/*   xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput)); */
/*   xf_Call(xf_FindOrCreatePrimalState(sampleAll, DoRestart, NULL, &sampleU)); */
/*   xf_Call(xf_CreateUniformTimeHistData(sampleAll, LogOutput, &TimeHistData)); */


/*   printf("i wonder i wonder\n"); */

/*   xf_Call(xf_CreateSolverData(sampleAll, &SolverData)); */
/*   SolverData->DoProcessResidual = xfe_False; // do not process residual */
/*   SolverData->UseArtificialDt = xfe_False; */
/*   SolverData->c = 0.; */

/*   xf_VectorGroup *RG = NULL; */
/*   printf("I wanna wish u\n"); */

/*   xf_Call(xf_FindSimilarVectorGroup(sampleAll, sampleU, "RG", */
/* 				    xfe_True, xfe_True, NULL, &RG, NULL)); */

/*   printf("feliz\n"); */
/*   xf_Call(xf_CalculateResidual(sampleAll, sampleU, RG, NULL, SolverData)); */
/*   /\* void yk_vector2Vec(xf_All *All, xf_Vector *A, Galaxy *pObj, Vec vecObj,          *\/ */
/*   /\*                  Is_it *reduced) *\/ */

/*   getchar(); */




/* } */


int yk_forwardSimulation(yk_Xflow *xflow, Cluster *fom){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int ierr;
  int nTimeStep, nTimeStep0, UW;
  real Time0;
  enum xfe_Bool DoRestart;
  char jobFile[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];    // prefix for saving files
  char LogOutput[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  xf_TimeHistData *TimeHistData=xflow->TimeHistData;
  char fName[xf_MAXSTRLEN];
  xf_All *sampleAll;
  char TimeHistFile[xf_MAXSTRLEN];
  clock_t tic;
  clock_t toc;
  //---------------------------------------------------------------------------
  // Find all the job and equation related information
  //---------------------------------------------------------------------------
  yk_textBoldBorder("Run Forward Order Model");
  //---------------------------------------------------------------------------
  // Burn Time to find the new initial conditions
  //---------------------------------------------------------------------------
  // Get the initial condiitons defined in the equations file
  xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart",
                             &DoRestart));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time", &Time0));
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "SavePrefix", SavePrefix));
   xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "SavePrefix", OutputFile));

  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput));
  xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG));
  /* // Create the Time History */
  xf_Call(xf_SetKeyValueInt(xflow->All->Param->KeyValue,
                            "UnsteadyWriteInterval", 1));

  if (fom->self->runFom == 1){
    xf_Call(xf_CreateUniformTimeHistData(xflow->All, LogOutput,
					 &xflow->TimeHistData));
    tic = clock();
    xf_Call(xf_ApplyTimeScheme(xflow->All, SavePrefix ,xfe_False, &xflow->UG,
			       xflow->TimeHistData));
    toc = clock();
    xf_Call(xf_SetKeyValueReal(xflow->All->Param->KeyValue, "Time", Time0));
    for (i=0; i<fom->self->time.count+1; i++){
      sprintf(OutputFile, "%s_U%d.xfa", SavePrefix, i);
      xf_Call(xf_WriteAllBinary(xflow->All, OutputFile));
    }
    yk_xfData2Dat(fom->self, xflow->UG, SavePrefix);

  }else if (fom->self->runFom == 0){
    sprintf(TimeHistFile, "naca_TimeHist.txt");
    ierr = xf_Error(xf_ReadTimeHistData(TimeHistFile, NULL,
                                        &xflow->TimeHistData));

  }
  fom->cpuTime = ((double)(toc-tic))/CLOCKS_PER_SEC;
  //---------------------------------------------------------------------------
  // Convert the fi;es fromd ata to ykflow form and create xfa txt files
  //---------------------------------------------------------------------------

  return ierr;
}

void yk_print2xflowTxt(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Galaxy *sol){
  int i;
  int ierr;
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;
  xf_VectorGroup *UG_approx;
  enum xfe_Bool DoRestart;
  xf_DataSet *DataSet = NULL;
  char Title[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  char temp[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  xf_Data *D;
  int leftNode = round(sol->time.t_0/sol->time.dt);
  int rightNode = round(sol->time.t_f/sol->time.dt);

  //---------------------------------------------------------------------------
  //
  //---------------------------------------------------------------------------
  xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "SavePrefix",
			 SavePrefix));
  sprintf(temp, "%s_U", sol->id);
  for (i=leftNode; i<rightNode+1; i++){
    xf_Call(xf_CreateDataSet(&DataSet));
    xf_Call(xf_DataSetAdd(DataSet, "State", xfe_VectorGroup, xfe_True, (void *)
			  xflow->UG, &D));
    ks_readSolution(multiquation->equation, sol, i);
    yk_array2VectorGroup(xflow->UG, sol);
    sprintf(Title, "%s%d.data", temp, i);
    xf_Call(xf_WriteDataSetBinary(xflow->All->Mesh, DataSet, NULL, Title));
    D->Data= NULL;
    sprintf(OutputFile, "%s%d.xfa", temp, i);
    xf_Call(xf_WriteAllBinary(xflow->All, OutputFile));
    xf_Call(xf_DestroyDataSet(DataSet));
  }

  sprintf(temp, "%s", sol->id);
  yk_xfData2Dat(sol, xflow->UG, temp);
}


void yk_initializeXflow(yk_PrimalSolver *ykflow){
  int ierr;
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  enum xfe_Bool DoRestart;
  char jobFile[xf_MAXSTRLEN];
   //Each Element is surrounded by three other triangles so dRdU has 4 non 0 val
  ykflow->numElemCol = 4;
  //Parse through the input file
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));

  //Read All from job and load and register the equations et
  xf_Call(xf_ReadAllFromJobFile(jobFile, xfe_True, &xflow->All));
  xf_Call(xf_LoadEqnSetLibrary(xflow->All->EqnSet->EqnSetLibrary)); //Equation
  xf_Call(xf_EqnSetRegister(xflow->All->EqnSet));
  //Make sure we have a primal state //initial condition is part of it
  xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart", &DoRestart));
  xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG));
  //Need to convert the xflow parameters into ykflow parameters

}


void yk_createSnapshotResidualFile(yk_PrimalSolver *ykflow,
				    Multiverse *multiquation, Cluster *fom,
				    Cluster *rom, Is_it *reduced){
  //Iterate through the number of parameters
  //For this 2D case lets make the paramset for residuals the same set
  // used for the states for POD snapshot at the beginning
  int ierr;
  int p;
  FILE *fin=NULL;
  int win_i = fom->self->time.window_i;
  char file_out[xf_MAXSTRLEN];
  FILE *fout=NULL;
  int ch;
  char fomBuffer[xf_MAXSTRLEN];
  char cwd[1024];
  char paramCalc[xf_MAXSTRLEN];
  char file_in[xf_MAXSTRLEN];
  char outputFile[xf_MAXSTRLEN];
  enum xfe_Bool DoRestart;
  char TimeHistFile[xf_MAXSTRLEN];
  char jobFile[xf_MAXSTRLEN];
  yk_PrimalSolver *ykflow_p;
  yk_Xflow *xflow_p;  //polymorphism
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;  //polymorphism
  Cluster *rom_p;
  struct stat sb = {0};
  struct stat st = {0};
  int tempSims = reduced->nSims;
  getcwd(cwd, sizeof(cwd));
  printf("%s\n", cwd);
  //---------------------------------------------------------------------------
  // Prepare the ultimate Residual snapshot file
  //---------------------------------------------------------------------------
  reduced->nSims = 1;
  // Create the ultimate residual snapshot file and open it to write into
  sprintf(file_out, "%s/residualSnapshot_%d.dat", cwd, win_i);
  fout = fopen(file_out, "w+");
  if (!fout) {
    perror("fopen output file");
    exit(EXIT_FAILURE);
  }
  //---------------------------------------------------------------------------
  //Iterate through all the oarams used to calcualte the residuals at
  //---------------------------------------------------------------------------
  //name for the Time History Files in the individual data directories
  sprintf(TimeHistFile, "naca_TimeHist.txt");
  //Extract the jobFile name to be assigned for residual creation
  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile));

  //Iterate through the parameter space
  for (p = 0; p< reduced->numParamSet; p++){
    //Create new xflow and ykflow instances
    ykflow_p = new_Xflow();
    xflow_p = (yk_Xflow *) ykflow_p->solver;
    //Create temporary cluster
    rom_p = (Cluster *) malloc (sizeof(Cluster));
    //Retrieve the data set name and enter that directory
    yk_fomBuffer(ykflow, multiquation->equation, reduced, p, fomBuffer);
    printf("Get Data from --->>>> %s\n", fomBuffer);
    if (stat(fomBuffer, &sb) == 0 && S_ISDIR(sb.st_mode)) {
      chdir(fomBuffer);
      //getcwd(cwd, sizeof(cwd));
    } else {
      perror(fomBuffer);
      exit(0);
    }
    //-------------------------------------------------------------------------
    //Initialize the xflow data struc
    //-------------------------------------------------------------------------
    xf_Call(xf_CreateKeyValue(&xflow_p->KeyValueArg));
    ierr = xf_AddKeyValue(xflow_p->KeyValueArg, "job", jobFile, xfe_True);
    strcpy(ykflow_p->path, ykflow->path);
    yk_initializeXflow(ykflow_p);
    xf_Call(xf_ReadTimeHistData(TimeHistFile, NULL,
					&xflow_p->TimeHistData));
    //-------------------------------------------------------------------------
    // Run through parameter space and calcualte ROM for each parameter set
    //-------------------------------------------------------------------------
    //Create directory to save data for each parameter space value
    sprintf(paramCalc, "%s/residualSnapshot_param_%d",  cwd, p);
    if (stat(paramCalc, &st) == -1)
      mkdir(paramCalc, 0700);
    //Copy over initial conditions from the given data to new folders

    if (fom->self->time.window_i == 0){
      sprintf(outputFile, "cp %s_%d.dat %s", fom->self->id, 0,  paramCalc);
      system(outputFile);
    }
    chdir(paramCalc);
    printf("FIRST TWO WRRKS\n");
    //Implement gauss newton solve to find the reduced states
    yk_runReducedOrderModel_ST(ykflow_p, multiquation, fom, rom_p, NULL, reduced);
    //Open the residual file GN solver wrote and copy to the master res file
    sprintf(file_in, "residualSnapshot_%d.dat", win_i);
    fin = fopen(file_in, "r");
    if (!fin) {
      perror("fopen input file");
      fprintf(stderr, "Problem file: %s\n", file_in);
      exit(EXIT_FAILURE);
    }
    while ((ch = fgetc(fin)) != EOF)
      fputc(ch, fout);
    fclose(fin);
    //-------------------------------------------------------------------------
    // Destroy everything for the particualr parameter set p
    //-------------------------------------------------------------------------
    free(rom_p->self->index);
    destroySystem(multiquation->equationReduced, rom_p->reduced);
    destroySystem(multiquation->equation, rom_p->self);
    free(rom_p->reduced);
    free(rom_p->self);
    xf_Call(xf_DestroyTimeHistData(xflow_p->TimeHistData));
    xf_Call(xf_DestroyKeyValue(xflow_p->KeyValueArg));
    xf_CloseEqnSetLibrary(&xflow_p->All);
    xf_Call(xf_DestroyAll(xflow_p->All));
    delete_Xflow(ykflow_p);
    free(rom_p);
  }
  reduced->nSims = tempSims;
  //Close the ultimate residual file, because it is done!
  fclose(fout);
  //Back to the original working directory
  chdir(cwd);
}


void yk_initializeYkflow(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			 Cluster *fom, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int ierr;
  int nTimeStep = 0;
  int argcIn = 2;
  int count = 0;
  int totstbasis;
  double *testparam;
  enum xfe_Bool DoRestart;
  char jobFile[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  char eqnReduced[xf_MAXSTRLEN];
  char inputFile[xf_MAXSTRLEN];
  char newFace[xf_MAXSTRLEN];
  char createParam[xf_MAXSTRLEN];
  char cwd[xf_MAXSTRLEN];
  char *argvIn[] = {"0", "9", NULL};
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  xf_All *All = xflow->All;
  const char* s = getenv("YKFLOWRUNPATH");

  //Each Element is surrounded by three other triangles so dRdU has 4 non 0 val
  /* ykflow->numElemCol = 4; */
  //---------------------------------------------------------------------------
  // Initiaize the forward order model to
  //---------------------------------------------------------------------------
  fom->self=(Galaxy *) malloc (sizeof(Galaxy));
  //---------------------------------------------------------------------------
  // Read through the YKFLOW input file
  //---------------------------------------------------------------------------

  xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "input", inputFile));
  sprintf(newFace, "%s.inp", inputFile);

  argvIn[1] = newFace;
  strcpy(multiquation->equation.nameEqn, inputFile);
  initIs_it(reduced);
  read_input(argcIn, argvIn, multiquation->equation, NULL, NULL, reduced);
  sprintf(createParam,"python3 createParam.py %g %g %g",
          reduced->params[0], reduced->params[1], reduced->params[2]);
  fom->self->runFom = reduced->runFom;
  system(createParam);

  printf("PATH :%s\n",(s!=NULL)? s : "getenv returned NULL");
  strcpy(ykflow->path, s);
  totstbasis = reduced->nTimeSegs*reduced->nSubWindows;


  reduced->final_nBasis_s = (PetscInt *) malloc (totstbasis*sizeof(PetscInt));
  reduced->final_nBasis_st = (PetscInt *) malloc (totstbasis*sizeof(PetscInt));
  reduced->final_Res_nBasis_s = (PetscInt *) malloc (reduced->nTimeSegs*
						     sizeof(PetscInt));
  reduced->final_Res_nBasis_st = (PetscInt *) malloc (reduced->nTimeSegs*
						      sizeof(PetscInt));
  reduced->r_t = (double *) malloc(reduced->nSims*sizeof(double));
  reduced->h_t = (double *) malloc(reduced->nSims*sizeof(double));
  for (i=0; i<reduced->nSims; i++){
    reduced->r_t[i] = 0;
    reduced->h_t[i] = 0;
  }
  //---------------------------------------------------------------------------
  // Read through the XFLOW input file
  //---------------------------------------------------------------------------
  yk_initializeXflow(ykflow);
  //Parse through the input file
  /* xf_Call(xf_GetKeyValue(xflow->KeyValueArg, "job", jobFile)); */

  /* //Read All from job and load and register the equations et */
  /* xf_Call(xf_ReadAllFromJobFile(jobFile, xfe_True, &xflow->All)); */
  /* xf_Call(xf_LoadEqnSetLibrary(xflow->All->EqnSet->EqnSetLibrary)); //Equation */
  /* xf_Call(xf_EqnSetRegister(xflow->All->EqnSet)); */
  /* //Make sure we have a primal state //initial condition is part of it */
  /* xf_Call(xf_GetKeyValueBool(xflow->All->Param->KeyValue, "Restart", &DoRestart)); */
  /* xf_Call(xf_FindOrCreatePrimalState(xflow->All, DoRestart, NULL, &xflow->UG)); */
  /* //Need to convert the xflow parameters into ykflow parameters */
  //---------------------------------------------------------------------------
  // Reduced mutliquation equation parameters
  //---------------------------------------------------------------------------
  sprintf(eqnReduced, "%s%s", inputFile, "Reduced");
  strcpy(multiquation->equationReduced.nameEqn, eqnReduced);
  //---------------------------------------------------------------------------
  // Define the Cluster properties here
  //---------------------------------------------------------------------------
  initUtype(fom->self);
  fom->self->time_method = 2;
  strcpy(fom->self->jobName, inputFile);
  multiquation->equation.numStates = xflow->UG->Vector[0]->StateRank;
  //---------------------------------------------------------------------------
  // Interpolation and quad points and basis shit
  //---------------------------------------------------------------------------
  fom->self->basis.p = xflow->UG->Vector[0]->Order[0];
  xf_Call(xf_GetKeyValueInt(xflow->All->Param->KeyValue, "nTimeStep",
                            &nTimeStep));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "EndTime",
			     &fom->self->time.t_f));
  xf_Call(xf_GetKeyValueReal(xflow->All->Param->KeyValue, "Time",
                             &fom->self->time.t_0));
  fom->self->time.globalT_f = fom->self->time.t_f;
  fom->self->time.globalT_0 = fom->self->time.t_0;

  fom->self->time.dt = (fom->self->time.t_f-fom->self->time.t_0)/nTimeStep;
  fom->self->time.count = nTimeStep;
  fom->self->J = (double *) malloc (2*sizeof(double));
  //---------------------------------------------------------------------------
  //Set the total number of elements FIne
  //---------------------------------------------------------------------------
  fom->self->space.elem.count = 0;
  for (i=0; i<xflow->All->Mesh->nElemGroup; i++){
    fom->self->space.elem.count += xflow->All->Mesh->ElemGroup[i].nElem;
  }//Set the number of basis nodes
  fom->self->basis.nodes = (xflow->UG->Vector[0]->GenArray)->r
    /xflow->UG->Vector[0]->StateRank;
  //Set the total number of nodes
  fom->self->space.node.count = fom->self->basis.nodes*
    fom->self->space.elem.count;

  fom->self->systemSize = (xflow->UG->Vector[0]->GenArray)->r*
    fom->self->space.elem.count;
  fom->self->index = (PetscInt *) malloc(fom->self->systemSize*
                                         sizeof(PetscInt));

  reduced->nSampleNodes = (int)(reduced->pSampleElems*
				fom->self->space.elem.count);
  for (i=0; i<fom->self->systemSize; i++)
    fom->self->index[i] = i;
  strcpy(fom->clusterId, "state");
  strcpy(fom->self->utypeId, "state");
  createSystem(multiquation->equation, fom->self);
  //Create the parameter set
  testparam = (double *) malloc (reduced->numParams*sizeof(double));
  reduced->paramMeshGrid =
    (double *) malloc (reduced->numParams*reduced->numParamSet*sizeof(double));
  definelocalMeshGrid(0, reduced, testparam,  &count);
  free(testparam);


}


void yk_RunErrorEstChaos(yk_PrimalSolver *ykflow, Multiverse *multiquation,
                         Cluster *fom, Cluster *rom, Cluster *hrom,
                         char *argv[], Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             //initialization for iteration
  int j_i;
  int ierr;
  yk_Xflow *xflow  = (yk_Xflow *) ykflow->solver;
  PetscInt sysSize = fom->self->systemSize;
  PetscInt time =fom->self->time.count;
  PetscInt nPerWindow = fom->self->time.count/reduced->nTimeSegs;
  PetscReal normROM;
  PetscReal normHROM;
  PetscReal normFOM;
  PetscReal resNormROM;
  PetscReal resNormHROM;
  PetscScalar val;
  Vec fomVec, romVec, hromVec;
  Mat fomMat, romMat, hromMat, resRomMat, resHromMat;
  enum xfe_Bool DoRestart;
  char windowString[xf_MAXSTRLEN];
  char cwd[xf_MAXSTRLEN];
  Vec delDrag;
  Vec delLift;
  double drag, lift;
  PetscReal normDrag, normLift;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  yk_VecCreateSeq(&fomVec, sysSize);
  yk_VecCreateSeq(&romVec, sysSize);

  yk_VecCreateSeq(&delDrag, time);
  yk_VecCreateSeq(&delLift, time);

  yk_MatCreateSeqDense(&fomMat, sysSize, time);
  yk_MatCreateSeqDense(&romMat, sysSize, time);
  yk_MatCreateSeqDense(&resRomMat, sysSize, time);
  if (reduced->runHrom == 1){
    yk_MatCreateSeqDense(&hromMat, sysSize, time);
    yk_VecCreateSeq(&hromVec, sysSize);
    yk_MatCreateSeqDense(&resHromMat, sysSize, time);
  }
  //---------------------------------------------------------------------------
  // Space-Time Reduced-Order Modeling
  //---------------------------------------------------------------------------
  reduced->JFlag = 0;
  yk_textBoldBorder("ykflow: Windowed Space-Time Model Order Reduction");
  for (i=0; i<reduced->nTimeSegs; i++){ //Number of time windows
    fom->self->time.window_i = i;
    sprintf(windowString, "W i n d o w  %d", i);
    yk_textBorder(windowString);
    //-------------------------------------------------------------------------
    // Reduced Order Modeling
    //-------------------------------------------------------------------------
    yk_textBoldBorder("Windowed ST Least-Squares Petrov-Galerkin");
    yk_textBorder("1. Create the Reduced Order Model");
    yk_createReducedOrderModel_ST(ykflow, multiquation, fom, reduced);
    yk_textBorder("2. Find the Reduced Order Model");
    yk_runReducedOrderModel_ST(ykflow, multiquation, fom, rom, reduced->r_t,
    			       reduced);
    //-------------------------------------------------------------------------
    // Hyper-Reduced Order Modeling
    //-------------------------------------------------------------------------
    if (reduced->runHrom == 1){
      yk_textBoldBorder("Windowed ST Gauss Newton w/ Approx Tensors");
      yk_textBorder("1. Create the hyper reduced order model");
      yk_createHyperReducedOrderModel_ST(ykflow, multiquation, fom, rom,reduced);
      yk_textBorder("2. Find the hyper reduced order model");
      yk_runHyperReducedOrderModel_ST(ykflow, multiquation, fom, hrom,
    				      reduced->h_t, reduced);
    }
    //-------------------------------------------------------------------------
    //Post process here
    //-------------------------------------------------------------------------
    reduced->hrom = 0;
    reduced->JFlag = 1;
    /* if (i==reduced->nTimeSegs-1){ */
    for (j=0; j<nPerWindow; j++){
      j_i = i*nPerWindow+j;
      //STATE
      ks_readSolution(multiquation->equation, fom->self, j_i+1);
      ks_readSolution(multiquation->equation, rom->self, j_i+1);

      array2Vec(multiquation->equation, fom->self, fom->self->space, fomVec);
      array2Vec(multiquation->equation, rom->self, fom->self->space, romVec);
      ks_Vec2MatCol(fomMat, sysSize, fom->self->index, j_i, fomVec,
		    INSERT_VALUES);
      ks_Vec2MatCol(romMat, sysSize, fom->self->index, j_i, romVec,
		    INSERT_VALUES);

      rom->self->time.node = j_i+1;

      yk_xflow_totalResidual(ykflow, multiquation, fom->self, romVec,
			     NULL, reduced, 0);

      yk_xflow_totalResidual(ykflow, multiquation, rom->self, romVec,
			     NULL, reduced, 0);
      drag=rom->self->J[0]-fom->self->J[0];
      lift=rom->self->J[1]-fom->self->J[1];
      VecSetValue(delDrag, j_i, drag, INSERT_VALUES);
      VecSetValue(delLift, j_i, lift, INSERT_VALUES);
      ks_Vec2MatCol(resRomMat, sysSize, fom->self->index, j_i, romVec,
		    INSERT_VALUES);

    /* 	if (reduced->runHrom==1){ */
    /* 	  ks_readSolution(multiquation->equation, hrom->stateFull, j_i+1); */
    /* 	  array2Vec(multiquation->equation, hrom->stateFull, fom->self->space, */
    /* 		    hromVec); */

    /* 	  ks_Vec2MatCol(hromMat, sysSize, fom->self->index, j_i, hromVec, */
    /* 			INSERT_VALUES); */
    /* 	  //RESIDUAL */
    /* 	  hrom->stateFull->time.node = j_i+1; */


    /* 	  yk_xflow_totalResidual(ykflow, multiquation, hrom->stateFull, hromVec, */
    /* 				 NULL, reduced, 0); */
    /* 	  ks_Vec2MatCol(resHromMat, sysSize, fom->self->index, j_i, hromVec, */
    /* 			INSERT_VALUES); */
    /* 	} */
    }

    yk_textBorder("Postprocess Data");
    yk_print2xflowTxt(ykflow, multiquation, rom->self);
    /* if (reduced->runHrom==1) */
    /*   yk_print2xflowTxt(ykflow, multiquation, hrom->stateFull); */
    //-------------------------------------------------------------------------
    // Destroy everything
    //-------------------------------------------------------------------------
    yk_destroyReducedOrderModel_ST(ykflow, multiquation, fom, rom, reduced);
    /* if (reduced->runHrom==1) */
    /*   yk_destroyHyperReducedOrderModel_ST(ykflow, multiquation,fom,hrom,reduced); */
    /* for (j=0; j<xflow->All->Mesh->nElemGroup; j++) */
  /*     xf_Release((void **) xflow->All->Mesh->ElemGroup[j].sElem); */
  }

  yk_textBoldBorder("Numerical Results");

  VecNorm(delDrag, NORM_2, &reduced->normDrag);
  VecNorm(delLift, NORM_2, &reduced->normLift);

  printf("HUMANS %g %g\n", reduced->normDrag, reduced->normLift);

  MatAssemblyBegin(fomMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(fomMat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(romMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(romMat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(resRomMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(resRomMat, MAT_FINAL_ASSEMBLY);

  MatAXPY(romMat, -1, fomMat, SAME_NONZERO_PATTERN);

  MatNorm(fomMat, NORM_FROBENIUS, &normFOM);
  MatNorm(romMat, NORM_FROBENIUS, &normROM);

  MatNorm(resRomMat, NORM_FROBENIUS, &reduced->resNormROM);


  reduced->relnormROM = normROM/normFOM;
  printf("%g\n", reduced->relnormROM);
  /* if (reduced->runHrom == 1){ */
  /*   MatAssemblyBegin(hromMat, MAT_FINAL_ASSEMBLY); */
  /*   MatAssemblyEnd(hromMat, MAT_FINAL_ASSEMBLY); */

  /*   MatAssemblyBegin(resHromMat, MAT_FINAL_ASSEMBLY); */
  /*   MatAssemblyEnd(resHromMat, MAT_FINAL_ASSEMBLY); */
  /*   MatAXPY(hromMat, -1, fomMat, SAME_NONZERO_PATTERN); */
  /*   MatNorm(hromMat, NORM_FROBENIUS, &normHROM); */
  /*   reduced->relnormHROM = normHROM/normFOM; */
  /*   MatNorm(resHromMat, NORM_FROBENIUS, &reduced->resNormHROM); */

  /* } */

  for (i=reduced->nSims-5; i<reduced->nSims; i++){
    reduced->aveR_t+=reduced->r_t[i];
    if (reduced->runHrom==1)
      reduced->aveH_t+=reduced->h_t[i];
  }


  reduced->aveR_t/=5.0;
  printf("Average Time fo find ROM: %g\n", reduced->aveR_t);

  if (reduced->runHrom ==1){
    reduced->aveH_t/=5.0;
    printf("Average Time fo find ROM: %g\n", reduced->aveH_t);
  }

  MatDestroy(&fomMat);
  MatDestroy(&romMat);
  MatDestroy(&resRomMat);
  VecDestroy(&fomVec);
  VecDestroy(&romVec);
  VecDestroy(&delDrag);
  VecDestroy(&delLift);
  if (reduced->hrom == 1){
    VecDestroy(&hromVec);
    MatDestroy(&hromMat);
    MatDestroy(&resHromMat);
  }
}

void yk_write2ExcelFile(Cluster *primal, Is_it *reduced){
  int i, j;
  int index_s = 0, index_t = 0, index_st = 0;
  int index_s_r = 0, index_t_r = 0, index_st_r = 0;
  int no_s=0, no_res=0;
  char str_res[xf_MAXSTRLEN];
  char str_res_st[xf_MAXSTRLEN];
  char str_res_t[xf_MAXSTRLEN];
  char nameOfxlsx[xf_MAXSTRLEN];
  char str_eBasisSpace[xf_MAXSTRLEN];
  char str_eBasisTime[xf_MAXSTRLEN];
  char str_eReBasisSpace[xf_MAXSTRLEN];
  char str_eReBasisTime[xf_MAXSTRLEN];
  char str[xf_MAXSTRLEN];
  char str_st[xf_MAXSTRLEN];
  char str_t[xf_MAXSTRLEN];
  int nSims =reduced->nSims;
  lxw_workbook  *workbook;
  lxw_worksheet *worksheet;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  sprintf(str_eBasisSpace, "%g", reduced->eBasisSpace);
  memmove(&str_eBasisSpace[0],&str_eBasisSpace[2],strlen(str_eBasisSpace)-1);
  sprintf(str_eBasisTime, "%g", reduced->eBasisTime);
  memmove(&str_eBasisTime[0], &str_eBasisTime[2], strlen(str_eBasisTime)-1);
  if (reduced->runHrom == 1){
    sprintf(str_eReBasisSpace, "%g", reduced->eReBasisSpace);
    memmove(&str_eReBasisSpace[0], &str_eReBasisSpace[2],
	    strlen(str_eReBasisSpace)-1);
    sprintf(str_eReBasisTime, "%g", reduced->eReBasisTime);
    memmove(&str_eReBasisTime[0], &str_eReBasisTime[2],
	    strlen(str_eReBasisTime) - 1);
    //set the name of the output file for xlsx for 2D at least
    sprintf(nameOfxlsx, "yk_navier_stokes_nTW_%d_%d_%g_%d_%s_%s_%s_%s.xlsx",
	    reduced->nTimeSegs, reduced->nSubWindows,
	    reduced->pSampleElems*100, reduced->nSampleTime,
	    str_eBasisSpace, str_eBasisTime, str_eReBasisSpace,
	    str_eReBasisTime);
  }else{
    sprintf(nameOfxlsx, "yk_navier_stokes_nTW_%d_%d_%s_%s.xlsx",
	    reduced->nTimeSegs, reduced->nSubWindows,
	    str_eBasisSpace, str_eBasisTime);
  }
  printf("%s\n", nameOfxlsx);
  workbook  = workbook_new(nameOfxlsx);
  worksheet = workbook_add_worksheet(workbook, NULL);
  //Format the strings that will contain the basis data
  //Space data
  for (i=0; i<reduced->nTimeSegs*reduced->nSubWindows; i++){
    index_s += sprintf(&str[index_s], "%d ",reduced->final_nBasis_s[i]);
    index_st += sprintf(&str_st[index_st],"%d ",reduced->final_nBasis_st[i]);
    index_t += sprintf(&str_t[index_t], "[");
    for (j=0; j<reduced->final_nBasis_s[i]; j++){
      index_t+=sprintf(&str_t[index_t],"%d",reduced->final_nBasis_t[no_s+j]);
      if (j<reduced->final_nBasis_s[i]-1)
	index_t +=sprintf(&str_t[index_t], " ");
    }
    no_s+=reduced->final_nBasis_s[i];
    index_t+=sprintf(&str_t[index_t], "]");
  }
  if (reduced->runHrom == 1){
    //Residual Data
    for (i=0 ; i<reduced->nTimeSegs; i++){
      index_s_r +=
	sprintf(&str_res[index_s_r], "%d ", reduced->final_Res_nBasis_s[i]);
      index_st_r +=
	sprintf(&str_res_st[index_st_r], "%d ", reduced->final_Res_nBasis_st[i]);
      index_t_r += sprintf(&str_res_t[index_t_r], "[");
      for (j=0; j<reduced->final_Res_nBasis_s[i]; j++){
	index_t_r += sprintf(&str_res_t[index_t_r], "%d",
			     reduced->final_Res_nBasis_t[no_res+j]);
	if (j<reduced->final_Res_nBasis_s[i]-1)
	  index_t_r +=sprintf(&str_res_t[index_t_r], " ");
      }
      no_res+=reduced->final_Res_nBasis_s[i];
      index_t_r+=sprintf(&str_res_t[index_t_r], "]");
    }
  }
  printf("letitgo\n");
  //---------------------------------------------------------------------------
  //Write to an excel file
  //---------------------------------------------------------------------------
  worksheet_write_string(worksheet, 0, 0, "Filename", NULL);
  worksheet_write_string(worksheet, 1, 0, nameOfxlsx, NULL);
  worksheet_write_string(worksheet, 0, 1, "nTW", NULL);
  worksheet_write_number(worksheet, 1, 1, reduced->nTimeSegs, NULL);
  worksheet_write_string(worksheet, 0, 2, "nSubPerW", NULL);
  worksheet_write_number(worksheet, 1, 2, reduced->nSubWindows, NULL);
  worksheet_write_string(worksheet, 0, 3, "lenWindow", NULL);
  worksheet_write_number(worksheet, 1, 3,
  			 primal->self->time.t_f/reduced->nTimeSegs, NULL);
  worksheet_write_string(worksheet, 0, 4, "lenSubWindow", NULL);
  worksheet_write_number(worksheet, 1, 4, primal->self->time.t_f
			 /reduced->nTimeSegs/
  			 reduced->nSubWindows, NULL);
  worksheet_write_string(worksheet, 0, 5, "ns/es", NULL);
  worksheet_write_number(worksheet, 1, 5, reduced->eBasisSpace, NULL);
  worksheet_write_string(worksheet, 0, 6, "nt/et", NULL);
  worksheet_write_number(worksheet, 1, 6, reduced->eBasisTime, NULL);
  worksheet_write_string(worksheet, 0, 7, "ns", NULL);
  worksheet_write_string(worksheet, 1, 7, str, NULL);
  worksheet_write_string(worksheet, 0, 8, "nt", NULL);
  worksheet_write_string(worksheet, 1, 8, str_t, NULL);
  worksheet_write_string(worksheet, 0, 9, "nst", NULL);
  worksheet_write_string(worksheet, 1, 9, str_st, NULL);
  worksheet_write_string(worksheet, 0, 10, "relative error", NULL);
  worksheet_write_number(worksheet, 1, 10, reduced->relnormROM, NULL);
  worksheet_write_string(worksheet, 0, 11, "residual", NULL);
  worksheet_write_number(worksheet, 1, 11,  reduced->resNormROM, NULL);
  worksheet_write_string(worksheet, 0, 12, "L2 Drag Error", NULL);
  worksheet_write_number(worksheet, 1, 12, reduced->normDrag, NULL);
  worksheet_write_string(worksheet, 0, 13, "L2 Lift Error", NULL);
  worksheet_write_number(worksheet, 1, 13, reduced->normLift, NULL);
  worksheet_write_string(worksheet, 0, 14, "nSim1", NULL);
  worksheet_write_number(worksheet, 1, 14, reduced->r_t[nSims-5], NULL);
  worksheet_write_string(worksheet, 0, 15, "nSim2", NULL);
  printf("letitbgo\n");
  worksheet_write_number(worksheet, 1, 15, reduced->r_t[nSims-4], NULL);
  worksheet_write_string(worksheet, 0, 16, "nSim3", NULL);
  worksheet_write_number(worksheet, 1, 16, reduced->r_t[nSims-3], NULL);
  worksheet_write_string(worksheet, 0, 17, "nSim4", NULL);
  worksheet_write_number(worksheet, 1, 17, reduced->r_t[nSims-2], NULL);
  worksheet_write_string(worksheet, 0, 18, "nSim5", NULL);
  worksheet_write_number(worksheet, 1, 18, reduced->r_t[nSims-1], NULL);
  worksheet_write_string(worksheet, 0, 19, "Average", NULL);
  worksheet_write_number(worksheet, 1, 19, reduced->aveR_t, NULL);
  if (reduced->runHrom == 1){
    worksheet_write_string(worksheet, 0, 20, "ns_r/es_r", NULL);
    worksheet_write_number(worksheet, 1, 20, reduced->eReBasisSpace, NULL);
    worksheet_write_string(worksheet, 0, 21, "nt_r/et_r", NULL);
    worksheet_write_number(worksheet, 1, 21, reduced->eReBasisTime, NULL);
    worksheet_write_string(worksheet, 0, 22, "sample time", NULL);
    worksheet_write_number(worksheet, 1, 22, reduced->nSampleTime, NULL);
    worksheet_write_string(worksheet, 0, 23, "sample space", NULL);
    worksheet_write_number(worksheet, 1, 23, reduced->nSampleNodes, NULL);
    worksheet_write_string(worksheet, 0, 24, "ns_r", NULL);
    worksheet_write_string(worksheet, 1, 24, str_res, NULL);
    worksheet_write_string(worksheet, 0, 25, "nt_r", NULL);
    worksheet_write_string(worksheet, 1, 25, str_res_t, NULL);
    worksheet_write_string(worksheet, 0, 26, "nst_r", NULL);
    worksheet_write_string(worksheet, 1, 26, str_res_st, NULL);
    worksheet_write_string(worksheet, 0, 27, "relative error gnat", NULL);
    worksheet_write_number(worksheet, 1, 27, reduced->relnormHROM, NULL);
    worksheet_write_string(worksheet, 0, 28, "residualgnat", NULL);
    worksheet_write_number(worksheet, 1, 28,  reduced->resNormHROM, NULL);
    worksheet_write_string(worksheet, 0, 29, "nSim1", NULL);
    worksheet_write_number(worksheet, 1, 29, reduced->h_t[nSims-5], NULL);
    worksheet_write_string(worksheet, 0, 30, "nSim2", NULL);
    worksheet_write_number(worksheet, 1, 30, reduced->h_t[nSims-4], NULL);
    worksheet_write_string(worksheet, 0, 31, "nSim3", NULL);
    worksheet_write_number(worksheet, 1, 31, reduced->h_t[nSims-3], NULL);
    worksheet_write_string(worksheet, 0, 32, "nSim4", NULL);
    worksheet_write_number(worksheet, 1, 32, reduced->h_t[nSims-2], NULL);
    worksheet_write_string(worksheet, 0, 33, "nSim5", NULL);
    worksheet_write_number(worksheet, 1, 33, reduced->h_t[nSims-1], NULL);
    worksheet_write_string(worksheet, 0, 34, "Average gnat", NULL);
    worksheet_write_number(worksheet, 1, 34, reduced->aveH_t, NULL);
  }
  workbook_close(workbook);


}

int main(int argc, char** argv) {
//------------------------------------//-------------------------------------
  //Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int ierr;
  int myRank, nProc;
  char *ArgIn[] = {"job", "NULL", ".job file name to read (run parameters)",
                   "input", "NULL", ".inp file name to read (run ROM params)",
                   "nWindow", "10", "starting number of time windows",
                   "\0"};
  Multiverse *multiquation = (Multiverse *) malloc (sizeof(Multiverse));
  Cluster *rom = (Cluster *) malloc (sizeof(Cluster));
  Cluster *fom = (Cluster *) malloc (sizeof(Cluster));
  Cluster *hrom = (Cluster *) malloc (sizeof(Cluster));
  Is_it *reduced =(Is_it*) malloc (sizeof(Is_it));
  yk_PrimalSolver *ykflow = new_Xflow();
  yk_Xflow *xflow = (yk_Xflow *) ykflow->solver;  //polymorphism
  char cwd[xf_MAXSTRLEN];
  //---------------------------------------------------------------------------
  // Intialize everything related to the xflow stuff
  //---------------------------------------------------------------------------

  SlepcInitialize(&argc, &argv, NULL, NULL);
  /* Py_Initialize(); */
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);

  // initialize key-value
  xf_Call(xf_CreateKeyValue(&xflow->KeyValueArg));
  /* // parse arguments */

  ierr = xf_ParseArg(ArgIn, argc, argv, xflow->KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
 /*  //--------------------------------------------------------------------------- */
 /*  // Initialization on ykflow */
 /*  //--------------------------------------------------------------------------- */
  yk_initializeYkflow(ykflow, multiquation, fom, reduced);
 /* //--------------------------------------------------------------------------- */
 /*  // Run the Forward Order Model (FOM) on xflow */
 /*  //--------------------------------------------------------------------------- */
  yk_forwardSimulation(xflow, fom);
 /*  //--------------------------------------------------------------------------- */
 /*  // Run Model Reduction Routines/ WST-LSPG */
 /*  //--------------------------------------------------------------------------- */
  if (reduced->runFomOnly == 0){
    yk_RunErrorEstChaos(ykflow, multiquation, fom, rom, hrom, argv, reduced);
 /*    //------------------------------------------------------------------------- */
 /*    // PostProcess everything and save data */
 /*    //------------------------------------------------------------------------- */
    yk_write2ExcelFile(fom, reduced);
  }
 /*  //--------------------------------------------------------------------------- */
 /*  // Destroy everything */
 /*  //--------------------------------------------------------------------------- */
 /*  /\* printf("CPU TIME\n"); *\/ */
 /*  /\* printf("%0.16e\n", fom->cpuTime); *\/ */
  /* Py_Finalize(); */

  xf_Call(xf_DestroyTimeHistData(xflow->TimeHistData));
  destroySystem(multiquation->equation, fom->self);
  free(fom->self->J);
  free(fom->self->index);
  free(fom->self);
  xf_CloseEqnSetLibrary(&xflow->All);
  xf_Call(xf_DestroyAll(xflow->All));

  xf_Call(xf_DestroyKeyValue(xflow->KeyValueArg));
  free(reduced->r_t);
  free(reduced->h_t);
  free(reduced->final_nBasis_s);
  free(reduced->final_nBasis_t);
  free(reduced->final_nBasis_st);
  free(reduced->final_Res_nBasis_s);
  free(reduced->final_Res_nBasis_t);
  free(reduced->final_Res_nBasis_st);
  free(reduced->win_i);
  free(reduced->params);
  free(reduced->paramsL);
  free(reduced->paramsH);
  free(reduced->dparams);
  free(reduced->paramMeshGrid);
  free(reduced);
  free(rom);
  free(fom);
  free(hrom);
  free(multiquation);
  delete_Xflow(ykflow);

  SlepcFinalize();

  return 0;

}



//SAVE THIS PART IT"S EXTRA BUT NECCESARY


  /* enum xfe_TimeSchemeType TimeScheme; */
  /* real iTime, Time; */
  /* int TimeStep; */

  /* char FileName[xf_MAXLINELEN]; */
  /* char SavePrefix[xf_MAXLINELEN]; */
  /* char LogOutput[xf_MAXLINELEN]; */
  /* xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "SavePrefix", SavePrefix)); */
  /* xf_Call(xf_GetKeyValue(xflow->All->Param->KeyValue, "LogOutput", LogOutput)); */
  /* /\* //--------------------------------------------------------------------------- *\/ */
  /* /\* // Reduced Order Modeling *\/ */
  /* /\* //--------------------------------------------------------------------------- *\/ */
  /* /\* printf("----------------------------------------------------------------\n"); *\/ */
  /* /\* printf("|                    Reduced-Order Modeling                    |\n"); *\/ */
  /* /\* printf("----------------------------------------------------------------\n"); *\/ */
  /* /\* yk_createReducedOrderModel(ykflow, multiquation, fom, rom, reduced); *\/ */
  /* /\* /\\* if (!reduced->restart) *\\/ *\/ */
  /* /\* /\\* yk_runReducedOrderModel(ykflow, multiquation, fom, rom, reduced); *\\/ *\/ */
  /* /\* //--------------------------------------------------------------------------- *\/ */
  /* /\* // Hyper-Reduced Order Modeling *\/ */
  /* /\* //--------------------------------------------------------------------------- *\/ */
  /* /\* printf("----------------------------------------------------------------\n"); *\/ */
  /* /\* printf("|                 Hyper-Reduced-Order Modeling                 |\n"); *\/ */
  /* /\* printf("----------------------------------------------------------------\n"); *\/ */
  /* /\* yk_createHyperReducedOrderModel(ykflow, multiquation, fom, hrom, reduced); *\/ */
  /* /\* yk_runHyperReducedOrderModel(ykflow, multiquation, fom, hrom, reduced); *\/ */
