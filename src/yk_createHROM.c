
#include "yk_createHROM.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <Python.h>
#include "ks_podv2.h"
/*****************************************************************************/
// FUNCTION Definition: yk_createReducedOrderModel_ST
//
// All offline aspects of ROMS is done in this function
// Only deals with primal data
/*****************************************************************************/

void yk_predictInitialCondition(Cluster *primal, Mat *snapshot,Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  //This function will be called for each subwindow.
  int i, j, k;
  int fd;
  char outputS[100];
  FILE* File;
  PetscInt sz, m, n;
  PetscInt we = 0, one = 1;
  PetscInt mu_com;
  PetscInt nParams = reduced->numParamSet;
  PetscInt sizeN = primal->self->systemSize;
  PetscInt snapshotCount = primal->self->time.count;
  PetscInt *value = (PetscInt *) malloc (sizeN*sizeof(PetscInt));
  PetscScalar *avec;
  PetscScalar *snapshot_mu = (PetscScalar *)malloc(sizeN*sizeof(PetscScalar));
  Vec initVec;
  Mat model;
  Mat A;
  Mat workA;
  Mat fullSnapshot;
  PetscViewer viewer, reader;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatGetSize(reduced->ST_rOBState, &m, &n);

  yk_VecCreateSeq(&initVec, n);
  printf("%d %d\n", sizeN*snapshotCount, nParams);
  yk_MatCreateSeqDense(&fullSnapshot, sizeN*snapshotCount, nParams);

  printf("%d %d\n", sizeN*snapshotCount, nParams);

  MatTranspose(reduced->ST_rOBState, MAT_INITIAL_MATRIX, &A);
  MatZeroEntries(A);

  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (k=0; k<snapshotCount; k++){ //time in each subwindow
    for (j=0; j<sizeN; j++)
      value[j] = k*sizeN+j;
    for (j=0; j<nParams; j++){
      mu_com = j*snapshotCount+k;
      MatGetValues(*snapshot, sizeN, primal->self->index, one, &mu_com,
		   snapshot_mu);
      MatSetValues(fullSnapshot, sizeN, value, 1, &j, snapshot_mu,
		   INSERT_VALUES);
    }
  }
  /* MatAssemblyBegin(fullSnapshot, MAT_FINAL_ASSEMBLY); */
  /* MatAssemblyEnd(fullSnapshot, MAT_FINAL_ASSEMBLY); */
  /* //Take the psuedo inverser of the entire basis matrix */
  /* moorePenrosePseudoInv(reduced->ST_rOBState, m, n, &A); */
  /* //Calculate the model in which we will predict the initial conditions */
  /* MatMatMult(A,fullSnapshot,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &model); */

  /* //Write the model to a binary file so that it can be read by python function */
  /* PetscViewerBinaryOpen(PETSC_COMM_SELF,"initialConditions.dat", */
  /* 			FILE_MODE_WRITE,&viewer); */
  /* MatView(model, viewer); */

  /* PetscViewerDestroy(&viewer); */
  Mat modelA;
  MatAssemblyBegin(fullSnapshot, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(fullSnapshot, MAT_FINAL_ASSEMBLY);
  //Take the psuedo inverser of the entire basis matrix
  printf("anymore\n");
  moorePenrosePseudoInvSparse(reduced->ST_rOBState, m, n, &A);
  printf("dsasdasd\n");
  MatMatMult(A,fullSnapshot,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &model);
  /* MatView(modelA, PETSC_VIEWER_STDOUT_SELF); */
  /* printf("IKSDJOIAHDFOISHDOIASD\n"); */
  /* getchar(); */
  /* moorePenrosePseudoInv(reduced->ST_rOBState, m, n, &workA); */
  /* //Calculate the model in which we will predict the initial conditions */
  /* MatMatMult(workA,fullSnapshot,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &model); */
  /* MatView(model, PETSC_VIEWER_STDOUT_SELF); */
  /* printf("adiushaiusdasd\n"); */

  /* MatView(A, PETSC_VIEWER_STDOUT_SELF); */
  /* getchar(); */
  /* MatView(workA, PETSC_VIEWER_STDOUT_SELF); */
  /* getchar(); */

  //Write the model to a binary file so that it can be read by python function
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"initialConditions.dat",
  			FILE_MODE_WRITE,&viewer);
  MatView(model, viewer);

  PetscViewerDestroy(&viewer);
  //---------------------------------------------------------------------------
  // Find the estimated initial conditions guess for the reduced states
  //---------------------------------------------------------------------------
  //Execute the python linear regression function to find the estimated r state
  /* sprintf(outputS, "main([%g, %g, %g, %g, %g, %g])", */
  /* 	  reduced->paramsL[0], reduced->paramsH[0], reduced->dparams[0], */
  /* 	  reduced->paramsL[1],reduced->paramsH[1], reduced->dparams[1] ); */

  /* PyRun_SimpleString("from runs.cases.readLinearRegressMat import main"); */
  /* PyRun_SimpleString(outputS); */

  sprintf(outputS, "python3 ../readLinearRegressMat.py %g %g %g %g %g %g",
  	  reduced->paramsL[0], reduced->paramsH[0], reduced->dparams[0],
  	  reduced->paramsL[1],reduced->paramsH[1], reduced->dparams[1] );
  system(outputS);

  //---------------------------------------------------------------------------
  // Read in the new predicted initial conditions to be used in this window
  //---------------------------------------------------------------------------
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"predictedInit.dat",FILE_MODE_READ,
			&reader);
  PetscViewerBinaryGetDescriptor(reader, &fd);
  //Create vector to put the new intiial conditions in
  VecCreate(PETSC_COMM_SELF, &reduced->init);
  VecSetSizes(reduced->init, PETSC_DECIDE,n);
  VecSetFromOptions(reduced->init);
  VecGetArray(reduced->init, &avec);

  //Read data into vector
  PetscBinaryRead(fd,&sz,1,NULL,PETSC_SCALAR);
  if (sz <=0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_UNEXPECTED,
  		      "Error: Must have array length > 0");
  PetscBinaryRead(fd,avec,sz,NULL,PETSC_SCALAR);
  VecRestoreArray(reduced->init, &avec);
  PetscViewerDestroy(&reader);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&fullSnapshot);
  MatDestroy(&model);
  MatDestroy(&A);
  VecDestroy(&initVec);
  free(value);
  free(snapshot_mu);
}

void yk_createReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				   Multiverse *multiquation, Cluster *primal,
				   Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i,j;                              // initialization for iteration
  PetscInt m;                           // temporary variable for robstate
  double dt = primal->self->time.dt;    // time step
  PetscInt rows, nBasis;                // nBasis = basis in spatial snapshot
  PetscInt snapshotCount;               // total temporal nodes to add
  PetscInt snapshotCount_win;           // total temporal node in window
  PetscInt snapshotCount_sub;           // elem in sub win
  PetscInt snapshotCount_0;             // number of time nodes before hand
  Vec init;                             // local initial conditions
  Mat *snapshot;                         // Snapshot Matrix before POD
  Mat *st_basis_temp;                   // temporal basis
  Mat *st_snapshot;                     // snapshot for the temporal
  Mat *s_mu;                            // same as snapshot but divided on mu
  Mat *subspaceTime;                    // subspace for each time window
  PetscInt win_i =primal->self->time.window_i;
  int info_n;
  PetscInt basis_0;
  char cwd[10000];
  Mat winSnapshot;
/* Mat **total_smu; */
  //---------------------------------------------------------------------------
  // Initialization (Mallocs and stuff)
  //---------------------------------------------------------------------------
  //Total number of tmeporal snapshots in the entire domain (no initial cond)
  snapshotCount = primal->self->time.count;
  /* //Total number of temporal sanpshots in each time window */
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
  //Total number of temporal snapshots in each sub window
  snapshotCount_sub = snapshotCount_win/reduced->nSubWindows;
  // Number of snapshots already finished before this current time window
  //ALlocate the subspaces and its corresponding basis for each window
  /* PetscMalloc1(reduced->nSubWindows, total_smu); */
  PetscMalloc1(reduced->nSubWindows, &snapshot);
  PetscMalloc1(reduced->nSubWindows, &subspaceTime);
  PetscMalloc1(snapshotCount_win, &reduced->ST_rOBState_i);
  snapshotCount_0 = primal->self->time.window_i * snapshotCount_win;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Prepare to find the basis for each subspace
  //Set the count for number of time stuff to deal with

  /* if (getcwd(cwd, sizeof(cwd)) != NULL) { */
  /*   printf("JUNE 18 important one: Current working dir: %s\n", cwd); */
  /* } else { */
  /*   perror("getcwd() error"); */
  /*   return 1; */
  /* } */
  primal->self->time.count  = snapshotCount_sub;
  // For each time window, we need to build the basis based on # of time window


  basis_0 = 0;
  for (i=0; i< win_i*reduced->nSubWindows;i++)
    basis_0+=reduced->final_nBasis_s[i];



  for (i=0; i<reduced->nSubWindows; i++){ //Iterate the number of sub win
    //Set where you want the offline stage to take place in time
    primal->self->time.t_0 = (snapshotCount_0+i*snapshotCount_sub)*dt;

    primal->self->time.t_f = (snapshotCount_0+(i+1)*snapshotCount_sub)*dt;
    //-------------------------------------------------------------------------
    // FOr eac subspace time window calcualte the original POD for space
    //-------------------------------------------------------------------------
    s_mu = yk_createSnapshotState(ykflow, multiquation, primal, &snapshot[i],
                                  reduced->dss, reduced);
    yk_properOrthogonalDecompose(&snapshot[i], primal->self->systemSize,
                                 primal->self->index, &reduced->nBasisFuncs,
                                 reduced->eBasisSpace, &reduced->rOBState);
    //-------------------------------------------------------------------------
    // Now need to calcualte the temmporal POD and et.c TIME
    //-------------------------------------------------------------------------
    MatGetSize(reduced->rOBState, &rows, &nBasis);
    printf("%lld\n", nBasis);
    printf("%d\n", i);
    PetscMalloc1(nBasis, &st_snapshot);
    PetscMalloc1(nBasis, &st_basis_temp);
    yk_createTemporalSnapshotState(ykflow, multiquation, primal,
				   reduced->rOBState, st_snapshot, s_mu,
				   reduced->dss, reduced->numParamSet,reduced);
    yk_properOrthogonalDecomposeGroup(st_snapshot, nBasis, snapshotCount_sub,
                                      primal->self->index,&reduced->nBasisTime,
                                      reduced->eBasisTime, st_basis_temp);
    //-------------------------------------------------------------------------
    // Calculate the ST-HOSVD Tailored space-time basis
    //-------------------------------------------------------------------------
    yk_createSpaceTimeBasis(&reduced->rOBState,st_basis_temp,&subspaceTime[i]);
    //-------------------------------------------------------------------------
    // Destriy Everything that must be within each subspace
    //-------------------------------------------------------------------------
    info_n = win_i*reduced->nSubWindows+i;
    MatGetSize(reduced->rOBState, &m, &reduced->final_nBasis_s[info_n]);
    //st basis stuff
    reduced->final_nBasis_t =
      yk_reallocPetscInt(reduced->final_nBasis_t, basis_0+nBasis);
    for (j=0; j<nBasis; j++){
      MatGetSize(st_basis_temp[j], &m, &reduced->final_nBasis_t[basis_0+j]);
    }
    basis_0+=nBasis;
    MatGetSize(subspaceTime[i], &m, &reduced->final_nBasis_st[info_n]);

    MatDestroy(&reduced->rOBState);
    for (j=0; j<reduced->numParamSet; j++) {MatDestroy(&s_mu[j]);}
    for (j=0; j<nBasis; j++){
      MatDestroy(&st_snapshot[j]);
      MatDestroy(&st_basis_temp[j]);
    }
    PetscFree(st_snapshot);
    PetscFree(st_basis_temp);
    PetscFree(s_mu);
  }
  //---------------------------------------------------------------------------
  //space time basis, and it's broken down version for each dt in subspace
  //---------------------------------------------------------------------------
  yk_formatSpaceTimeBasis(primal, subspaceTime, &reduced->ST_rOBState,
			  reduced->ST_rOBState_i, reduced);
  //-------------------------------------------------------------------------
  // Add the space-time basis to the whole space-time matrix the ultimate one
  //-------------------------------------------------------------------------
  primal->self->time.count = snapshotCount_win;
  primal->self->time.t_0 = primal->self->time.window_i*snapshotCount_win*dt;
  yk_createSnapshotState(ykflow, multiquation, primal, &winSnapshot,
				reduced->dss, reduced);
  printf("porkie\n");
  yk_predictInitialCondition(primal, &winSnapshot, reduced);

  //---------------------------------------------------------------------------
  // Save important ingo to save
  //---------------------------------------------------------------------------
  /* MatGetSize(reduced->rOBState, &m, &n); */
  /* MatGetSize(reduced->ST_rOBState, &m, &n); */
  /* if (getcwd(cwd, sizeof(cwd)) != NULL) { */
  /*   printf("OVER JUNE 18 important one: Current working dir: %s\n", cwd); */
  /* } else { */
  /*   perror("getcwd() error"); */
  /*   return 1; */
  /* } */

  primal->self->time.count = snapshotCount; //rese
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->nSubWindows; i++){ MatDestroy(&subspaceTime[i]);}
  for (i=0; i<reduced->nSubWindows; i++){ MatDestroy(&snapshot[i]);}
  PetscFree(subspaceTime);
  PetscFree(snapshot);
  MatDestroy(&winSnapshot);
}

void yk_runReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Cluster *primalApprox, double *time_spent,
				Is_it *reduced){
  //Need to put back the time spent here and put the dummy variable back
  // int he prog gile.
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  double dt = primal->self->time.dt;    // time step
  PetscInt i;                           // initialization for iteration
  PetscInt leftTimeNode;                // left time node for current t window
  PetscInt innertimesolver = 1;         // time solver, dummy variable
  PetscInt row, nBasis;                 // size of the entire space time
  Universe *_eqnFull = &multiquation->equation;
  Universe *_eqnRcd = &multiquation->equationReduced;
  char approximateID[1000];             // Name given to the approximate state
  char reducedID[1000];                 // Name given to the reduced state
  PetscInt snapshotCount;               // num of snapshot in entire domain
  PetscInt snapshotCount_win;           // nunber of snapshot in entire window
  Vec init;                             // information about the initial cond
  char residualSnapshot[1000];
  clock_t tic, toc;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->nSnapshotsRJ = 0;
  MatGetSize(reduced->ST_rOBState, &row, &nBasis);
  //Create and close the residual data file to fill
  sprintf(residualSnapshot, "%s_%d.dat", "residualSnapshot",
	  primal->self->time.window_i);
  fclose(fopen(residualSnapshot, "w"));
  //Total number of tmeporal snapshots in the entire domain (no initial cond)
  snapshotCount = primal->self->time.count;
  //Total number of temporal sanpshots in each time window
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
 //Total number of temporal snapshots in each sub window
  strcpy(primalApprox->clusterId, primal->clusterId);
  //Set the number of states for the reduced solutions
  multiquation->equationReduced.numStates = nBasis;
  //Set the count to the number of dt in the time window for set up
  primal->self->time.count = snapshotCount_win;
  primal->self->time.t_0 = primal->self->time.window_i*snapshotCount_win*dt;
  primal->self->time.t_f=(primal->self->time.window_i+1)*snapshotCount_win*dt;
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  strcpy(approximateID, "STROM_approximate");
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->self->systemSize = primal->self->systemSize;
  primalApprox->self->index = (PetscInt *) malloc
    (primalApprox->self->systemSize*sizeof(PetscInt));
  for (i=0; i<primalApprox->self->systemSize; i++)
    primalApprox->self->index[i] = i;
  createSystem(*_eqnFull, primalApprox->self);
  reduced->reducedMesh = primalApprox->self->space;
  initArrayInt(&primalApprox->self->time.tNode_i, snapshotCount_win);
  for (i=0; i<snapshotCount_win; i++)
    primalApprox->self->time.tNode_i.array[i] = i;
  initArrayInt(&primalApprox->self->time.time_Os, snapshotCount_win);
  for (i=0; i<snapshotCount_win; i++)
    primalApprox->self->time.time_Os.array[i] = i;
  primalApprox->self->J = (double *) malloc (2*sizeof(double));
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", (*_eqnFull).nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states (ROM) from the reduced solution
  //---------------------------------------------------------------------------
  strcpy(reducedID, "STROM_reduced");
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(*_eqnRcd, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", (*_eqnFull).nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------

  leftTimeNode = round(primal->self->time.t_0/primal->self->time.dt);
  if (leftTimeNode == 0){
    ks_readSolution(*_eqnFull, primal->self, 0);
    ks_copySolutions(*_eqnFull, primal->self, primalApprox->self);
    ks_printSolution(*_eqnFull, primalApprox->self, 0);
  }else
    ks_readSolution(*_eqnFull, primalApprox->self, leftTimeNode);
  //---------------------------------------------------------------------------
  // Calculate the initial conditions for the reduced space states
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->nSims; i++){
    printf("Simulation Number %d\n", i);
    vec2Array(*_eqnRcd, primalApprox->reduced,
	      primalApprox->reduced->space, reduced->init);
    ks_printSolution(*_eqnRcd, primalApprox->reduced,
		     primal->self->time.window_i);
    //-------------------------------------------------------------------------
    // Gauss Newton solver to minimize 1/2 norm(Residua, 2)^2
    //-------------------------------------------------------------------------
    tic = clock();
    yk_gaussNewtonSolve_ST(ykflow, multiquation, primalApprox, reduced,
    			   innertimesolver, &reduced->nSnapshotsRJ);
    //Set the time count back to to the number of dt in the entire time domain
    toc = clock();
    if (time_spent != NULL)
      time_spent[i] += (double)(toc - tic) / CLOCKS_PER_SEC;
  }
  primal->self->time.count = snapshotCount;
  //---------------------------------------------------------------------------
  // Destroy Everything that is only relevant for this current function
  //---------------------------------------------------------------------------
  delArrayInt(&primalApprox->self->time.tNode_i);
  delArrayInt(&primalApprox->self->time.time_Os);

}

void yk_destroyReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				    Multiverse *multiquation, Cluster *primal,
				    Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;                                //initialization for iteration
  int snapshotCount = primal->self->time.count/reduced->nTimeSegs;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Implementation aka destroy everything
  //---------------------------------------------------------------------------
  for (i=0; i<snapshotCount; i++)
    MatDestroy(&reduced->ST_rOBState_i[i]);
  VecDestroy(&reduced->init);
  PetscFree(reduced->ST_rOBState_i);
  MatDestroy(&reduced->ST_rOBState);
  free(primalApprox->self->index);
  destroySystem(_equationReduced, primalApprox->reduced);
  destroySystem(multiquation->equation, primalApprox->self);
  free(primalApprox->reduced);
  free(primalApprox->self->J);
  free(primalApprox->self);
}

void yk_precomputeCoefMat(Mat Z, Is_it *reduced){
  //------------------------------------//------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//------------------------------------
  //Z is the identity matrix in space and time
  PetscInt m, n;
  int col;
  Mat coefInter;
  MatGetSize(reduced->ST_rOBResidual, &m, &n);
  col = n;
  int what = reduced->reducedMesh.elem.count*12*(reduced->nSampleTime-1);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  MatMatMult(Z, reduced->ST_rOBResidual, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
	     &coefInter);

  moorePenrosePseudoInv(coefInter, what, col, &reduced->A);
  /* MatGetSize(reduced->A, &m, &n); */
  /* printf("%d %d\n", m, n); */
  /* getchar(); */
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&coefInter);
}

void yk_createHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow,
					Multiverse *multiquation,
					Cluster *primal, Cluster *primalApprox,
					Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  char residualSnapshots[1000];
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  PetscInt numTime = primal->self->time.count;
  PetscInt systemSize = primal->self->systemSize;
  PetscInt st_systemSize;
  PetscInt *rowIndex = (PetscInt *) malloc (systemSize*sizeof(PetscInt ));
  PetscInt *temporalSet;
  PetscInt *spatialSet;
  PetscInt *mshOs;
  PetscInt snapshotCount;               // total temporal nodes to add
  PetscInt snapshotCount_win;           // total temporal node in window
  PetscInt snapshotCount_sub;           // elem in sub win
  PetscInt snapshotCount_0;             // number of time nodes before hand
  Mat *sRes_mu;
  PetscInt m,n;
  Mat Z, Za, ZZ;
  Mat *st_snapshot;                     // snapshot for the temporal
  Mat *s_mu;                            // same as snapshot but divided on mu
  Mat *st_basis_temp;                    // subspace for each time window
  PetscInt nBasis;
  PetscInt rows;
  Mat subspaceTime;
  Mat rSnapshot;
  Mat Q,R;
  PetscInt colRes;
  int RJLines;
  double dt = primal->self->time.dt;    // time step
  PetscInt info_n;
  PetscInt win_i=primal->self->time.window_i;
  PetscInt basis_0 = 0;
  for (i=0; i< win_i; i++)
    basis_0+=reduced->final_Res_nBasis_s[i];
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Find the residual basis matrix
  //---------------------------------------------------------------------------
  ykflow->residualSnapshotFile(ykflow, multiquation, primal, primalApprox,
			       reduced);

  sprintf(residualSnapshots, "%s_%d.dat", "residualSnapshot",
          primal->self->time.window_i);

  for (i=0; i<systemSize; i++) {rowIndex[i] = i;}
  printf("*750\n");
  //Total number of tmeporal snapshots in the entire domain (no initial cond)
  snapshotCount = primal->self->time.count;
  //Total number of temporal sanpshots in each time window
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
  //Total number of temporal snapshots in each sub window
  /* snapshotCount_sub = snapshotCount_win/reduced->nSubWindows; */
 //Set the number of snapshots for each time subwindow
  // Number of snapshots already finished before this current time window
  snapshotCount_0 = primal->self->time.window_i * snapshotCount_win;
  PetscMalloc1(snapshotCount_win, &reduced->ST_rOBStateBar_i);
  primal->self->time.count  = snapshotCount_win;
  st_systemSize = primal->self->systemSize*snapshotCount_win;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  /* primal->self->time.t_0 = snapshotCount_0; */
  /* primal->self->time.t_f = (snapshotCount_0+snapshotCount_win)*dt; */
  primal->self->time.t_0 = round(primal->self->time.window_i*snapshotCount_win*dt);
  primal->self->time.t_f=round((primal->self->time.window_i+1)*snapshotCount_win*dt);
  //---------------------------------------------------------------------------
  //Readm the snapshot residualmatrix file for manipulation
  //---------------------------------------------------------------------------
  sRes_mu = ks_readMatrix(residualSnapshots, systemSize, snapshotCount_win,
			  &rSnapshot, &RJLines, reduced);
  yk_properOrthogonalDecompose(&rSnapshot, systemSize, rowIndex,
			       &reduced->nBasisFuncsRJ,
			       reduced->eReBasisSpace,&reduced->rOBResidual);
  //-------------------------------------------------------------------------
  // Now need to calcualte the temmporal POD and et.c TIME
  //-------------------------------------------------------------------------

  MatGetSize(rSnapshot, &rows, &colRes);
  MatGetSize(reduced->rOBResidual, &rows, &nBasis);
  PetscMalloc1(nBasis, &st_snapshot);
  PetscMalloc1(nBasis, &st_basis_temp);
  yk_createTemporalSnapshotState(ykflow, multiquation, primal,
				 reduced->rOBResidual, st_snapshot,
				 sRes_mu, reduced->dss, RJLines, reduced);
  MatGetSize(rSnapshot, &m, &n);
  MatGetSize(*st_snapshot, &m, &n);

  yk_properOrthogonalDecomposeGroup(st_snapshot, nBasis, snapshotCount_win,
				    primal->self->index,&reduced->nBasisTime,
				    reduced->eBasisTime, st_basis_temp);

  //---------------------------------------------------------------------------
  // Create the ultimate  space-time residual matrix
  //---------------------------------------------------------------------------
  yk_createSpaceTimeBasis(&reduced->rOBResidual,st_basis_temp, &subspaceTime);
  MatGetSize(reduced->rOBResidual, &m, &n);
  MatGetSize(*st_basis_temp, &m, &n);
  MatGetSize(subspaceTime, &m, &n);

  MatCreateSeqDense(PETSC_COMM_SELF, st_systemSize, n, NULL,
		    &reduced->ST_rOBResidual);
  MatCreateSeqDense(PETSC_COMM_SELF, n, n, NULL, &R);
  qRFactorization(subspaceTime, reduced->ST_rOBResidual, R);
  //---------------------------------------------------------------------------
  // Hyper reductiom initialization
  //---------------------------------------------------------------------------
  //Use greedy to find which time indices to keep
  yk_greedyAlgorithm_temporalSet(primal, &temporalSet, reduced);
  //Use greedy to find which spatial indices to keep
  yk_greedyAlgorithm_spatialSet(primal, &spatialSet, &temporalSet, &Z,reduced);
  minElements(multiquation, primal, spatialSet, temporalSet, &Za, &ZZ,reduced);
  ykflow->adjacentElems(ykflow, primal, &mshOs, &temporalSet, reduced);
  //---------------------------------------------------------------------------
  // Create coefficient for hyper reduction
  //---------------------------------------------------------------------------
  yk_precomputeCoefMat(ZZ, reduced);
  //---------------------------------------------------------------------------
  // Prepare the basis for the hyper reduction
  //---------------------------------------------------------------------------
  MatMatMult(reduced->Z_Os, reduced->ST_rOBState, MAT_INITIAL_MATRIX,
	     PETSC_DEFAULT, &reduced->ST_rOBStateBar);
  for (i=0; i<snapshotCount_win; i++)
    MatMatMult(reduced->spaceZ, reduced->ST_rOBState_i[i], MAT_INITIAL_MATRIX,
	       PETSC_DEFAULT, &reduced->ST_rOBStateBar_i[i]);


  MatGetSize(reduced->rOBResidual, &m, &reduced->final_Res_nBasis_s[win_i]);
  MatGetSize(reduced->ST_rOBResidual, &m, &reduced->final_Res_nBasis_st[win_i]);

  reduced->final_Res_nBasis_t =
    yk_reallocPetscInt(reduced->final_Res_nBasis_t, basis_0+nBasis);
  for (i=0; i<nBasis; i++){
    MatGetSize(st_basis_temp[i], &m, &reduced->final_Res_nBasis_t[basis_0+i]);
  }
  basis_0+=nBasis;
  //-------------------------------------------------------------------------
  // Add the space-time basis to the whole space-time matrix the ultimate one
  //-------------------------------------------------------------------------
  primal->self->time.count = snapshotCount; //reset it

  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  for (i=0; i<nBasis; i++){MatDestroy(&st_basis_temp[i]);}
  for (i=0; i<nBasis; i++){MatDestroy(&st_snapshot[i]);}
  for (i =0 ; i<RJLines; i++){MatDestroy(&sRes_mu[i]);}
  MatDestroy(&subspaceTime);
  PetscFree(sRes_mu);
  PetscFree(st_basis_temp);
  PetscFree(st_snapshot);
  MatDestroy(&Za);
  MatDestroy(&R);
  MatDestroy(&ZZ);
  MatDestroy(&rSnapshot);
  free(mshOs);
  free(rowIndex);
  free(temporalSet);
  free(spatialSet);
  MatDestroy(&Z);
}


void yk_runHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow,
				     Multiverse *multiquation, Cluster *primal,
				     Cluster *primalApprox, double *time_spent,
				     Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  int size_Os = reduced->reducedMesh_Os.node.count;
  int samSysSize = size_Os*_equationFull.numStates;
  double dt = primal->self->time.dt;
  char approximateID[1000];             // Name given to the approximate state
  char reducedID[1000];                 // Name given to the reduced state
  char approximateFinalID[1000];
  PetscInt i, j;
  PetscInt leftTimeNode;
  PetscInt innertimesolver = 1;              //time solver
  PetscInt row, nBasis;
  PetscInt snapshotCount;
  PetscInt snapshotCount_win;
  Vec primal_0;
  Vec primalReducedVec;
  Vec init;
  clock_t tic, toc;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------

  reduced->hrom = 1;
  reduced->nSnapshotsRJ = 0;

  MatGetSize(reduced->ST_rOBState, &row, &nBasis);
  //Total number of tmeporal snapshots in the entire domain (no initial cond)
  snapshotCount = primal->self->time.count;
  //Total number of temporal sanpshots in each time window
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
  //Total number of temporal snapshots in each sub window
  strcpy(primalApprox->clusterId, primal->clusterId);
  /*  printf("words\n"); */
  multiquation->equationReduced.numStates = nBasis;

  primal->self->time.count = snapshotCount_win;
  primal->self->time.t_0 = primal->self->time.window_i * snapshotCount_win*dt;
  primal->self->time.t_f=(primal->self->time.window_i+1)*snapshotCount_win*dt;
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  strcpy(approximateID, "STGNAT_approximate");
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  ks_copyUtype(primal->self, primalApprox->self);
  createSystem(multiquation->equation, primalApprox->self);
  initArrayInt(&primalApprox->self->time.tNode_i, reduced->nSampleTime-1);
  for (i=0; i<primalApprox->self->time.tNode_i.count; i++)
    primalApprox->self->time.tNode_i.array[i] = reduced->reducedTime.array[i];
  initArrayInt(&primalApprox->self->time.time_Os, reduced->reducedTime_Os.count);
  for (i=0; i<primalApprox->self->time.time_Os.count; i++)
    primalApprox->self->time.time_Os.array[i] = reduced->reducedTime_Os.array[i];


  primalApprox->self->systemSize = primal->self->systemSize;
  primalApprox->self->index = (PetscInt *) malloc
    (primalApprox->self->systemSize*sizeof(PetscInt));
  for (i=0; i<primalApprox->self->systemSize; i++)
    primalApprox->self->index[i] = i;
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states (ROM) from the reduced solution
  //---------------------------------------------------------------------------
  strcpy(reducedID, "STGNAT_reduced");
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(multiquation->equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  strcpy(approximateFinalID, "STGNAT_approximateFinal");
  primalApprox->stateFull = (Galaxy *) malloc (sizeof(Galaxy));
  ks_copyUtype(primal->self, primalApprox->stateFull);
  primalApprox->stateFull->systemSize = primal->self->systemSize;
  primalApprox->stateFull->index = (PetscInt *) malloc
    (primalApprox->stateFull->systemSize*sizeof(PetscInt));
  for (i=0; i<primalApprox->stateFull->systemSize; i++)
    primalApprox->stateFull->index[i] = i;
  createSystem(multiquation->equation, primalApprox->stateFull);
  sprintf(primalApprox->stateFull->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->stateFull->basis.p,
          primalApprox->stateFull->time_method, approximateFinalID);
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------
  leftTimeNode = primal->self->time.t_0/dt;
  if (leftTimeNode == 0){
    ks_readSolution(_equationFull, primal->self, leftTimeNode);
    ks_copySolutions(_equationFull, primal->self, primalApprox->self);
    ks_copySolutions(_equationFull, primal->self, primalApprox->stateFull);
    ks_printSolution(_equationFull, primalApprox->self, leftTimeNode);
    ks_printSolution(_equationFull, primalApprox->stateFull, leftTimeNode);
  }else
    ks_readSolution(_equationFull, primalApprox->self, leftTimeNode);



  //---------------------------------------------------------------------------
  // Calculate the initial conditions for the reduced space states
  //---------------------------------------------------------------------------
  for (i=0; i<reduced->nSims; i++){
    vec2Array(multiquation->equationReduced, primalApprox->reduced,
	      primalApprox->reduced->space, reduced->init);
    ks_printSolution(multiquation->equationReduced, primalApprox->reduced,
		     primal->self->time.window_i);
    //--------------------------------------------------------------------------
    // Gauss NEwton to minimize 1/2 norm(Residua, 2)^2
    //--------------------------------------------------------------------------
    tic = clock();
    yk_gaussNewtonSolve_ST(ykflow, multiquation, primalApprox, reduced,
			   innertimesolver, &reduced->nSnapshotsRJ);
    toc = clock();
    time_spent[i] += (double)(toc - tic) / CLOCKS_PER_SEC;
  }
  //Set the time count back to to the number of dt in the entire time domain
  primal->self->time.count = snapshotCount;
  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  VecDestroy(&reduced->init);
}


void yk_destroyHyperReducedOrderModel_ST(yk_PrimalSolver *ykflow,
                                    Multiverse *multiquation, Cluster *primal,
                                    Cluster *primalApprox, Is_it *reduced){
  int i; //initialization for iteration
  int snapshotCount = primal->self->time.count/reduced->nTimeSegs;
  Universe _equationReduced = multiquation->equationReduced;
  int snapshotCount_win;
  int snapshotCount_sub;

  //---------------------------------------------------------------------------
  // Destroy everything
  //---------------------------------------------------------------------------
  MatDestroy(&reduced->A);
  MatDestroy(&reduced->rOBResidual);
  MatDestroy(&reduced->spaceZ);
  MatDestroy(&reduced->Z_Os);
  snapshotCount = primal->self->time.count;
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
  //Total number of temporal sanpshots in each time window
  snapshotCount_win = snapshotCount/reduced->nTimeSegs;
  snapshotCount_sub = snapshotCount_win/reduced->nSubWindows;
  for (i=0; i<snapshotCount_win; i++)
    MatDestroy(&reduced->ST_rOBStateBar_i[i]);
  MatDestroy(&reduced->ST_rOBStateBar);
  MatDestroy(&reduced->rOBResidual);
  MatDestroy(&reduced->ST_rOBResidual);
  destroySystem(_equationReduced, primalApprox->reduced);
  destroySystem(multiquation->equation, primalApprox->self);
  destroySystem(multiquation->equation, primalApprox->stateFull);
  PetscFree(reduced->ST_rOBStateBar_i);
  free(primalApprox->self->index);
  free(primalApprox->stateFull->index);
  delArrayInt(&primalApprox->self->time.tNode_i);
  delArrayInt(&primalApprox->self->time.time_Os);
  delArrayInt(&reduced->reducedTime);
  delArrayInt(&reduced->reducedTime_Os);
  delArrayInt(&reduced->reducedMesh_Os.node);
  delArrayInt(&reduced->reducedMesh_Os.elem);
  free(reduced->reducedMesh_Os.i_ElemCom);
  delArrayInt(&reduced->reducedMesh.elem);
  delArrayInt(&reduced->reducedMesh.node);
  free(primalApprox->reduced);
  free(primalApprox->self);
  free(primalApprox->stateFull);
}


void yk_createHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				     Multiverse *multiquation, Cluster *primal,
				     Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  int seedCount;                        // Number of Boundary Nodes
  int *nodeSet = (int *) malloc (reduced->nSampleNodes*sizeof(int));
  char approximateID[1000];
  char reducedID[1000];
  char realID[1000];
  char *residualSnapshots = "residualSnapshot.dat";
  char *jacobianSnapshots = "jacobianSnapshot.dat";
  Mat sResMat, sJacMat;
  Mat Z;
  PetscInt m,n;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 1;
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->stateFull = (Galaxy *) malloc (sizeof(Galaxy));
  strcpy(approximateID, "HROM_approximate");
  strcpy(reducedID, "HROM_reduced");
  strcpy(realID, "HROM_full");
  strcpy(primalApprox->clusterId, primal->clusterId);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for residuals
  //---------------------------------------------------------------------------
  ks_readMatrix(residualSnapshots, primal->self->systemSize, seedCount, &Z, &m, reduced);
  MatGetSize(sResMat, &m, &n);
  printf("residualsnapshot %d %d\n", m, n);
  getchar();
  yk_properOrthogonalDecompose(&sResMat, primal->self->systemSize,
			       primal->self->index, reduced->nBasisFuncsRJ,
			       reduced->eBasisSpace, &reduced->rOBResidual);
  //---------------------------------------------------------------------------
  // Create reduced model calculation for jacobians
  //---------------------------------------------------------------------------
  ks_readMatrix(jacobianSnapshots, primal->self->systemSize, primal->self->systemSize,&sJacMat, &n, reduced);
  yk_properOrthogonalDecompose(&sJacMat, primal->self->systemSize,
  			       primal->self->index, reduced->nBasisFuncsRJ,
  			       reduced->eBasisSpace, &reduced->rOBJacobian);
  printf("POD DONE\n");
  //---------------------------------------------------------------------------
  // Greedy Algorithm and create of the reduced mesh
  //---------------------------------------------------------------------------
  ykflow->boundarySeeds(ykflow, multiquation, primal, reduced->nSampleNodes,
			&seedCount, nodeSet);
  printf("here\n");
  ks_greedyAlgorithm(seedCount, nodeSet, reduced);
  printf("mayube here\n");
  Mat ZZ;
  minElements(multiquation, primal, nodeSet, nodeSet, &Z, &ZZ,reduced);
  for (i=0; i<reduced->reducedMesh.elem.count; i++)
    printf("%d\n", reduced->reducedMesh.elem.array[i]);
  //---------------------------------------------------------------------------
  // Create the hyper approximate states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->stateFull->systemSize = primal->self->systemSize;
  primalApprox->stateFull->index = (int *) malloc
    (primalApprox->stateFull->systemSize*sizeof(int));
  for (i=0; i<primalApprox->stateFull->systemSize; i++)
    primalApprox->stateFull->index[i] = i;
  findApproxMesh(ykflow, multiquation, primal, primalApprox, reduced);
  assignBasisQuad(primalApprox->self);
  initializeStates(_equationFull, primalApprox->self);
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(_equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Create the hyper real states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->stateFull);
  createSystem(_equationFull, primalApprox->stateFull);
  sprintf(primalApprox->stateFull->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->stateFull->basis.p,
          primalApprox->stateFull->time_method, realID);
  //---------------------------------------------------------------------------
  // Create the corresponding identity matrices from greedy algorithm
  //---------------------------------------------------------------------------
  restrictedIdentity(_equationFull, primal->self, primalApprox->self, nodeSet,
		     reduced);
  MatMatMult(reduced->Zp, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
             &reduced->rOBStateBar);
  MatMatMult(reduced->Zu, reduced->rOBState, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
             &reduced->rOBStateRed);
  //---------------------------------------------------------------------------
  // Destroy and Free things
  //---------------------------------------------------------------------------
  free(nodeSet);
  MatDestroy(&sResMat);
  MatDestroy(&sJacMat);
}

void yk_runHyperReducedOrderModel(yk_PrimalSolver *ykflow,
				  Multiverse *multiquation, Cluster *primal,
				  Cluster *primalApprox,Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int systemApproxSize;                 //Size of the primalApprox vector
  int innertimesolver = 1;              //Only for 1D... gotta figure out this
  Vec primalApprox_0;
  Vec primal_0;
  Vec primalReducedVec;
  Vec primalRealVec;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialization and Mallocing
  //---------------------------------------------------------------------------
  reduced->hrom = 1;
  systemApproxSize =
    _equationFull.numStates*primalApprox->self->space.node.count;
  VecCreate(PETSC_COMM_SELF, &primal_0);
  VecSetSizes(primal_0, primal->self->systemSize, primal->self->systemSize);
  VecSetType(primal_0, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalRealVec);
  VecSetSizes(primalRealVec, primal->self->systemSize,
	      primal->self->systemSize);
  VecSetType(primalRealVec, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, reduced->nBasisFuncs, reduced->nBasisFuncs);
  VecSetType(primalReducedVec, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalApprox_0);
  VecSetSizes(primalApprox_0, systemApproxSize, systemApproxSize);
  VecSetType(primalApprox_0, VECSEQ);
  //---------------------------------------------------------------------------
  // Find A and B for the calculation of residual and jacobian approx
  //---------------------------------------------------------------------------
  ks_offlineMatrix(_equationFull, primal->self, primalApprox->self,
		   reduced);
  //---------------------------------------------------------------------------
  // Print out the initial conditions for completeness (primalApprox);
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  array2Vec(_equationFull, primal->self, primal->self->space, primal_0);
  MatMult(reduced->Zp, primal_0, primalApprox_0);
  vec2Array(_equationFull, primalApprox->self, primalApprox->self->space,
	    primalApprox_0);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // Find the r0 initial solution of reduced solution
  //---------------------------------------------------------------------------


  VecSet(primalReducedVec, 0);
  vec2Array(_equationReduced, primalApprox->reduced,
	    primalApprox->reduced->space, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);


  /* array2Vec(_equationFull, primal->self, primal_0); */
  /* MatMultTranspose(reduced->rOBState, primal_0, primalReducedVec); */
  /* vec2Array(_equationReduced, primalApprox->reduced, primalReducedVec); */
  /* ks_printSolution(_equationReduced, primalApprox->reduced, 0); */

  //---------------------------------------------------------------------------
  // Read Solution initial conditions
  //---------------------------------------------------------------------------
  ks_copySolutions(_equationFull, primal->self, primalApprox->stateFull);
  ks_printSolution(_equationFull, primalApprox->stateFull, 0);
  //---------------------------------------------------------------------------
  // Gauss Newton Solver
  //---------------------------------------------------------------------------
  yk_gaussNewtonSolve(ykflow, multiquation, primalApprox, reduced,
		      innertimesolver, &reduced->nSnapshotsRJ);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  VecDestroy(&primal_0);
  VecDestroy(&primalRealVec);
  VecDestroy(&primalReducedVec);
  VecDestroy(&primalApprox_0);
  MatDestroy(&reduced->A);
  MatDestroy(&reduced->B);
}



//-----------------------------------------------------------------------------
// Initialize and set up everything for the potential to run a ROM
//-----------------------------------------------------------------------------
void yk_createReducedOrderModel(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i;
  char approximateID[1000];             // Name given to the approximate state
  char reducedID[1000];                 // Name given to the reduced state
  Mat snapshot;                         // Snapshot Matrix before POD
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  primalApprox->self = (Galaxy *) malloc (sizeof(Galaxy));
  primalApprox->reduced = (Galaxy *) malloc (sizeof(Galaxy));
  strcpy(approximateID, "ROM_approximate");
  strcpy(reducedID, "ROM_reduced");
  strcpy(primalApprox->clusterId, primal->clusterId);

  //---------------------------------------------------------------------------
  // Create reduced model calculation for states
  //---------------------------------------------------------------------------
  yk_createSnapshotState(ykflow, multiquation, primal,&snapshot, reduced->dss,
			 reduced);
  yk_properOrthogonalDecompose(&snapshot, primal->self->systemSize,
			       primal->self->index, reduced->nBasisFuncs,
			       reduced->eBasisSpace, &reduced->rOBState);
  //---------------------------------------------------------------------------
  // Create the approximated value based on the reduced solution
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->self);
  primalApprox->self->systemSize = primal->self->systemSize;
  primalApprox->self->index = (int *) malloc
    (primalApprox->self->systemSize*sizeof(int));
  for (i=0; i<primalApprox->self->systemSize; i++)
    primalApprox->self->index[i] = i;
  createSystem(multiquation->equation, primalApprox->self);
  reduced->reducedMesh = primalApprox->self->space;
  sprintf(primalApprox->self->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->self->time_method, approximateID);
  //---------------------------------------------------------------------------
  // Create the reduced states
  //---------------------------------------------------------------------------
  ks_copyUtype(primal->self, primalApprox->reduced);
  primalApprox->reduced->basis.p = 0;
  primalApprox->reduced->basis.nodes = 1;
  primalApprox->reduced->space.dx = 0;
  primalApprox->reduced->space.elem.count = 1;
  primalApprox->reduced->space.node.count = 1;
  createSystem(_equationReduced, primalApprox->reduced);
  sprintf(primalApprox->reduced->id, "%s_%s_%d_%d_%s", _equationFull.nameEqn,
          primalApprox->clusterId, primalApprox->self->basis.p,
          primalApprox->reduced->time_method, reducedID);
  //---------------------------------------------------------------------------
  // Destroy Everything
  //---------------------------------------------------------------------------
  MatDestroy(&snapshot);
}

//-----------------------------------------------------------------------------
// Execute to find the reduced order model
//-----------------------------------------------------------------------------
void yk_runReducedOrderModel(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			     Cluster *primal, Cluster *primalApprox,
			     Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int innertimesolver = 1;              //time solver
  Vec primal_0;
  Vec primalReducedVec;
  Universe _equationFull = multiquation->equation;
  Universe _equationReduced = multiquation->equationReduced;
  //---------------------------------------------------------------------------
  // Initialize everything
  //---------------------------------------------------------------------------
  reduced->hrom = 0;
  reduced->nSnapshotsRJ = 0;
  VecCreate(PETSC_COMM_SELF, &primal_0);//Set up the singular vector
  VecSetSizes(primal_0, primal->self->systemSize, primal->self->systemSize);
  VecSetType(primal_0, VECSEQ);
  VecCreate(PETSC_COMM_SELF, &primalReducedVec);
  VecSetSizes(primalReducedVec, reduced->nBasisFuncs, reduced->nBasisFuncs);
  VecSetType(primalReducedVec, VECSEQ);

  fclose(fopen("residualSnapshot.dat", "w"));
  fclose(fopen("jacobianSnapshot.dat", "w"));
  //---------------------------------------------------------------------------
  // Find the initial conditions for the approximated states and print
  //---------------------------------------------------------------------------
  ks_readSolution(_equationFull, primal->self, 0);
  ks_copySolutions(_equationFull, primal->self, primalApprox->self);
  ks_printSolution(_equationFull, primalApprox->self, 0);
  //---------------------------------------------------------------------------
  // Set r0 initial solution of reduced solution to 0
  //---------------------------------------------------------------------------
  /* array2Vec(_equationFull, primal->self, primal_0); */
  /* MatMultTranspose(reduced->rOBState, primal_0, primalReducedVec); */
  VecSet(primalReducedVec, 0);
  vec2Array(_equationReduced, primalApprox->reduced,
	    primalApprox->reduced->space, primalReducedVec);
  ks_printSolution(_equationReduced, primalApprox->reduced, 0);
  //---------------------------------------------------------------------------
  // Gauss NEwton to minimize 1/2 norm(Residua, 2)^2
  //---------------------------------------------------------------------------
  yk_gaussNewtonSolve(ykflow, multiquation, primalApprox, reduced,
		      innertimesolver, &reduced->nSnapshotsRJ);
  VecDestroy(&primal_0);
  VecDestroy(&primalReducedVec);

}


void minElements(Multiverse *multiquation, Cluster *primal, int *nodeSet, int *temporalSet,
		 Mat *spaceIdentity, Mat *ZZ, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             // initialization for iteration
  int elem;                             // elem of interest given nSampleNodes
  int count = 0;
  int totNumElems = 0;                  // saves total number of elements to c
  int numBasis = primal->self->basis.nodes;
  int elemSize = numBasis*multiquation->equation.numStates;

  int *tempNodeUSet = (int *) malloc     // array info on which element to keep
    (primal->self->space.elem.count*sizeof(int));
  int sysTime = primal->self->time.count;
  int timeCount = reduced->nSampleTime-1;
  Mat timeIdentity;
  MatCreate(PETSC_COMM_SELF, &timeIdentity);
  MatSetSizes(timeIdentity, timeCount, sysTime, timeCount, sysTime);
  MatSetType(timeIdentity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(timeIdentity, 1, NULL);
  //Fill in the identity matrix based on the time indices we want to keep
  for (j=0; j<timeCount; j++)
    MatSetValue(timeIdentity, j, temporalSet[j], 1, INSERT_VALUES);
  MatAssemblyBegin(timeIdentity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(timeIdentity, MAT_FINAL_ASSEMBLY);



  //---------------------------------------------------------------------------
  // Implementation given the greedyalgorithm giving me nodes to keep
  //---------------------------------------------------------------------------
  for (i=0; i<primal->self->space.elem.count; i++)
    tempNodeUSet[i] = 0;
  for (i=0; i<reduced->nSampleNodes*12; i++){ //tells me which element included
    elem = floor(nodeSet[i]/(elemSize)); //nodeSet has only samplenodes
    if (tempNodeUSet[elem] == 0){
      tempNodeUSet[elem] = 1;             // include this element
    }
  }
  /* // Count the total of elements in the reduced mesh */
  for (i=0; i<primal->self->space.elem.count; i++)
    totNumElems += tempNodeUSet[i];
  /* //Need to fill in has table thingie in xflow */
  totNumElems = reduced->nSampleNodes;
  //---------------------------------------------------------------------------
  // Create the actual temporary mesh
  //---------------------------------------------------------------------------
  //Take this out later, this is weird
  //Don't need the dx/x thing. Want to be indepedent on the physical mesh trait
  //Make sure this is possible
  reduced->reducedMesh.elem.count = totNumElems;
  //total number of physics nodes
  reduced->reducedMesh.node.count = totNumElems*numBasis;

  initArrayInt(&reduced->reducedMesh.elem, totNumElems);
  initArrayInt(&reduced->reducedMesh.node, numBasis*totNumElems);
  //the array in in elem and nodes arenumbered relative to the resiudal vector


  MatCreate(PETSC_COMM_SELF, spaceIdentity);
  MatSetSizes(*spaceIdentity, numBasis*totNumElems*4, primal->self->systemSize,
	      numBasis*totNumElems*4, primal->self->systemSize);
  MatSetType(*spaceIdentity, MATSEQAIJ);
  MatSeqAIJSetPreallocation(*spaceIdentity, 1, NULL);



  for (i=0; i<primal->self->space.elem.count; i++){
    if (tempNodeUSet[i] == 1){ //If this element exists in reduced mesh
      reduced->reducedMesh.elem.array[count] = i;
      for (j=0; j<numBasis; j++){
  	reduced->reducedMesh.node.array[count*numBasis+j] = i*numBasis+j;
      }
      for (j=0; j<12; j++){
	MatSetValue(*spaceIdentity, count*12+j, i*12+j, 1, INSERT_VALUES);
      }
      count++;
    }

  }

  MatAssemblyBegin(*spaceIdentity, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*spaceIdentity, MAT_FINAL_ASSEMBLY);

  free(tempNodeUSet);

  yk_kkron(timeIdentity, *spaceIdentity, ZZ);

  MatDestroy(&timeIdentity);

}

void findApproxMesh(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *primal, Cluster *primalApprox, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j;                             // initialization for iteration
  int totNumElems = 0;
  int count = 0;
  int elemCount = primal->self->space.elem.count;
  int numBasis = primal->self->basis.nodes;
  int *elemSet = (int *) malloc (elemCount*sizeof(int));
  //---------------------------------------------------------------------------
  // Count the number of total elements in the approximate mesh given reduced m
  //---------------------------------------------------------------------------

  for (i=0; i<elemCount; i++)
    elemSet[i] = 0; //initialize it to empty 0
  /* for (i=0; i<reduced->reducedMesh.elem.count; i++) */
  /*   ykflow->adjacentElems(ykflow, primal, reduced->reducedMesh.elem.array[i], */
  /* 			  elemSet); */
  /* for (i=0; i<reduced->reducedMesh.elem.count; i++) */
  /*   printf(" */
  printf("kitchen\n");

  /* //Count the number of total elements in the approximate mesh */
  /* for (i=0; i<elemCount; i++) */
  /*   totNumElems += elemSet[i]; */
  /* //--------------------------------------------------------------------------- */
  /* // Now create that approximate mesh */
  /* //--------------------------------------------------------------------------- */
  /* initArrayInt(&primalApprox->self->space.elem, totNumElems); */
  /* initArrayInt(&primalApprox->self->space.node, numBasis*totNumElems); */
  /* //Don't need the dx/x thing. Want to be indepdent off the physical mesh */
  /* for (i=0; i<primal->self->space.elem.count; i++){ */
  /*   if (elemSet[i] == 1){ //If it exists */
  /*     primalApprox->self->space.elem.array[count] = i; */
  /*     for (j=0; j<numBasis; j++) */
  /*       primalApprox->self->space.node.array[count*numBasis+j] = i*numBasis+j; */
  /*     count++; */
  /*   } */
  /* } */
  free(elemSet);
}


void yk_findBoundarySeeds1D(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			    Cluster *primal, int nSampNodes, int *numSeeds,
			    int *nodeSet){
  *numSeeds = 2;
  nodeSet[0] = 0;
  nodeSet[multiquation->equation.numStates] =  primal->self->systemSize-1;
}

void yk_findAdjacentElems1D(yk_PrimalSolver *ykflow, Cluster *primal, int elem,
			    int *elemSet){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  //1D
  elemSet[elem]= 1;
  if (elem >0 )
    elemSet[elem-1] = 1;
  if (elem < primal->self->space.elem.count-1)
    elemSet[elem+1] = 1;
}


void createInitialConditions_ST(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal,
				Mat *st_rOBTemp, Vec *estimatedinitST,
				Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int mach = 60;
  int alfa = round(reduced->params[1]);
  /* int alfa = 9; */
  int reynolds = round(reduced->params[2]);
  int timeNode0 = round(primal->self->time.t_0/primal->self->time.dt);
  char nameOfDir[1000];  //Name of the directory to go into , just rounding
  char cwd[1024];
  PetscInt i,j; //initialization for iteration
  PetscInt *index = (PetscInt *) malloc (primal->self->systemSize*
					 sizeof(PetscInt));
  PetscInt row, nBasis;
  PetscInt col;
  PetscInt snapshotCount = primal->self->time.count;
  Vec vecsnapshot; //vector containing the snapshots
  Mat A;

  PetscScalar *state_i = (PetscScalar *) malloc (primal->self->systemSize*
						 sizeof(PetscScalar));
  PetscScalar *state_0 = (PetscScalar *) malloc (primal->self->systemSize*
						 sizeof(PetscScalar));
  char tempS[1000];
  struct stat sb = {0};
  /* if (alfa % 2 == 0) */
  /*   alfa -= 1; */
  //---------------------------------------------------------------------------
  // initialization
  //---------------------------------------------------------------------------

  MatGetSize(*st_rOBTemp, &row, &nBasis);

  //Create the snapshot in vector form
  VecCreate(PETSC_COMM_SELF, &vecsnapshot);
  VecSetSizes(vecsnapshot, row, row);
  VecSetType(vecsnapshot, VECSEQ);

  //Create the reduced states that we are looking for (initial condtions)
  VecCreate(PETSC_COMM_SELF, estimatedinitST);
  VecSetSizes(*estimatedinitST, nBasis, nBasis);
  VecSetType(*estimatedinitST, VECSEQ);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //Find the file name of the folder we want to go into and set it to nameOfDir

  sprintf(nameOfDir, "%s/%s_M_%d_A_%d_Re_%d", ykflow->path,
	  multiquation->equation.nameEqn, 700, 800, reynolds);
  getcwd(cwd, sizeof(cwd));
  printf("CWD ITNERESTED %s\n", cwd);
  printf("%s\n", nameOfDir);
  //Go into the folder now
  if (stat(nameOfDir, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    chdir(nameOfDir);
  } else {
    perror(nameOfDir);
    exit(0);
  }
  printf("%d\n", timeNode0);
  //Read the state solutions for the initial time in that directory
  ks_readSolution(multiquation->equation, primal->self, timeNode0);
  array2Data(multiquation->equation, primal->self, state_0);
  //Add data to the Petsc Vec Type from the snapshots data
  //Iterate through the number of time steps within the time window/simulation
  for (i=0; i<primal->self->time.count; i++){  //creating snapshot matrix
    //Remember to skip the initial condtiions, we already did it
    ks_readSolution(multiquation->equation, primal->self, timeNode0+i+1);
    //Convert the intial condtiion to a double array
    array2Data(multiquation->equation, primal->self, state_i);
    //Take away the initial conditions from the states aka center it man
    for (j=0; j<primal->self->systemSize; j++){
      state_i[j]-=state_0[j];
      //Also set the indices for the vector placement here
      index[j] = i*primal->self->systemSize + j;
    }
    //Now set the centered snapshot to the vector here
    VecSetValues(vecsnapshot, primal->self->systemSize, index, state_i,
		 INSERT_VALUES);
  }
  VecAssemblyBegin(vecsnapshot);
  VecAssemblyEnd(vecsnapshot);

  /* MatGetSize(*st_rOBTemp, &row, &col); */
  /* printf("fact %d %d\n", row, col); */
  /* VecGetSize(vecsnapshot, &row); */
  /* printf("fact %d \n", row); */

  /* MatMultTranspose(*st_rOBTemp, vecsnapshot, *estimatedinitST); */

  moorePenrosePseudoInv(*st_rOBTemp, row, nBasis, &A);

  MatMult(A, vecsnapshot, *estimatedinitST);
  chdir(cwd);
  /* strcpy(primal->self->id, tempS); */
  VecDestroy(&vecsnapshot);
  MatDestroy(&A);
  /* VecView(*estimatedinitST, PETSC_VIEWER_STDOUT_SELF); */
  free(index);
  free(state_0);
  free(state_i);
}
