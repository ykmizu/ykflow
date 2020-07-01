#include "ks_read_solution_file.h"
//-----------------------------------------------------------------------------
// Restrieves the solution file and extracts it to Utype solution.
// Also designates the node number the spatial solution refers to
//
// Yukiko Shimizu
// August 3, 2018
//-----------------------------------------------------------------------------


int ks_checkSolutionExist(Universe eqn, Galaxy *U, int node){
  /* int i, j;             //initialization for iteration */
  /* FILE* solution_file; */
  /* char outputFile[1000]; */
  /* sprinf(outputFile, "%s_%d.dat", U->id, node); */
  /* if( access( outputFile, F_OK ) != -1 ) */
  /*   return 1; */
  /* else */
  /*   return 0; */
  /* return 0; */

}

void ks_printSolution(Universe eqn, Galaxy* U, int node){
  int i, j;             //initialization for iteration
  FILE* solution_file;
  char outputFile[1000];
  char cwd[1000];
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  sprintf(outputFile, "%s_%d.dat", U->id, node);
  /* printf("print %s\n", outputFile); */
  /* if (getcwd(cwd, sizeof(cwd)) != NULL) { */
  /*   printf("important one: Current working dir: %s\n", cwd); */
  /* } else { */
  /*   perror("getcwd() error"); */
  /*   return 1; */
  /* } */
  solution_file = fopen(outputFile, "w+");
  if (solution_file == NULL){
    perror("Error");
    exit(EXIT_FAILURE);
  }
  for (i=0; i<U->space.node.count; i++){
    fprintf(solution_file, "%0.15E", U->solution[0].array[i]);
    for (j=1; j<eqn.numStates; j++)
      fprintf(solution_file, " %0.15E", U->solution[j].array[i]);
    fprintf(solution_file, "\n");
  }
  fclose(solution_file);
}

void ks_printReducedSolution(Universe eqn, Galaxy *U, int node){
  int i;
  FILE* solutionFile;
  char outputFile[1000];
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  sprintf(outputFile, "%s_%d.dat", U->id, node);
  solutionFile = fopen(outputFile, "w+");
  for (i=0; i<U->space.node.count; i++)
    fprintf(solutionFile, "0 %0.16e\n", U->solution->array[i]);
  fclose(solutionFile);
}

void ks_removeSolution(Universe eqn, Galaxy *U, int node){
  char outputFile[1000];
  sprintf(outputFile, "%s_%d.dat", U->id, node);
  remove(outputFile);
}




void ks_readSolution(Universe eqn, Galaxy* U, int node){
  int i;                //initialization for iteration
  FILE* solutionFile;
  char outputFile[1000];
  //char buffer[10000];
  int count = 0;
  char *current_line;
  char *lines = NULL;
  char cwd[1000];
  char *c_lines = NULL;
  int length;
  char *buffer = (char *) malloc (eqn.numStates*23*sizeof(char));
  /* printf("%d\n", eqn.numStates*22); */
  U->time.node = node;
  sprintf(outputFile, "%s_%d.dat", U->id, node);
  /* printf("%s\n", outputFile); */
  /* if (getcwd(cwd, sizeof(cwd)) != NULL) { */
  /*   printf("important one: Current working dir: %s\n", cwd); */
  /* } else { */
  /*   perror("getcwd() error"); */
  /*   return 1; */
  /* } */

  solutionFile = fopen(outputFile, "r");
  if (solutionFile == NULL){
    perror("Error");
    exit(EXIT_FAILURE);
  }else{

    while ((current_line = fgets(buffer, eqn.numStates*23, solutionFile))){
      if (current_line[0] != '\n'){
	lines = strtok(current_line, " ");
	//printf("die die die %s\n", lines);
	U->solution[0].array[count] = atof(lines);
	for (i=1; i<eqn.numStates; i++){
	  c_lines = strtok(NULL, " ");
	  length = strlen(c_lines);
	  if (c_lines[length-1] == ' ' || lines[length-1] == '\n')
	    c_lines[length-1] = '\0';
	  U->solution[i].array[count] = atof(c_lines);
	}
	count ++;
      }
    }
  }
  fclose(solutionFile);
  free(buffer);
}

void yk_readLSSInitialCon(char *filename, double R[]){
  //R is the initial condition for the LSS equations
  int count = 0;
  double initLSS;
  FILE *initialLSSFile;
  initialLSSFile = fopen(filename, "r");
  if (initialLSSFile == NULL){
    perror("Error");
    exit(EXIT_FAILURE);
  }else{
    while(fscanf(initialLSSFile, "%lf\n", &initLSS)!=EOF){
      R[count] = initLSS;
      count++;
    }
    fclose(initialLSSFile);
  }
}

void yk_printLSSInitialCon(int size, double R[], int out_i){
  int i;
  FILE *initialLSSFile;
  char outputFile[1000];
  sprintf(outputFile, "%s_%d.out", "initLSS", out_i);
  initialLSSFile = fopen(outputFile, "w");
  if (initialLSSFile == NULL){
    perror("Error");
    exit(EXIT_FAILURE);
  }else{
    for (i=0; i<size; i++)
      fprintf(initialLSSFile, "%0.16f\n", R[i]);
    fclose(initialLSSFile);
  }
}

Mat * ks_readMatrix(char *filename, int systemSize, int nTime,
		    Mat *resSnapshot, int *RJLinesCount, Is_it *reduced){
  int i, j, k;
  FILE *residualSnapshotFile;
  char *ith_line;
  int count = 0;
  double *snapshotArray;
  Mat *sResidual_mu;
  PetscScalar *bigAssSnapshotArray = NULL;
  PetscInt counter = 0;
  int counter0;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  /* for (i=0; i<reduced->numParamSet; i++){ */
  //Open up the residual snapshot file to view
  residualSnapshotFile = fopen(filename, "r");
  yk_fgets(residualSnapshotFile, systemSize*nTime, &snapshotArray,
	   RJLinesCount);

  fclose(residualSnapshotFile);
  PetscMalloc1(*RJLinesCount, &sResidual_mu);


  MatCreate(PETSC_COMM_SELF, resSnapshot);
  MatSetSizes(*resSnapshot,  systemSize, nTime*(*RJLinesCount), systemSize,
	      nTime*(*RJLinesCount));
  MatSetType(*resSnapshot, MATSEQDENSE);
  MatSetUp(*resSnapshot);


  for (i=0; i<*RJLinesCount; i++){
    MatCreate(PETSC_COMM_SELF, &sResidual_mu[i]);
    MatSetSizes(sResidual_mu[i], systemSize, nTime, systemSize, nTime);
    MatSetType(sResidual_mu[i], MATSEQDENSE);
    MatSetUp(sResidual_mu[i]);

    for (j=0; j<nTime; j++){
      for (k=0; k<systemSize; k++){
	MatSetValue(sResidual_mu[i], k, j, snapshotArray[(i*nTime+j)*systemSize+k],
		    INSERT_VALUES);
	MatSetValue(*resSnapshot, k, i*nTime+j, snapshotArray[(i*nTime+j)*systemSize+k],
		    INSERT_VALUES);
      }
    }
    MatAssemblyBegin(sResidual_mu[i], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sResidual_mu[i], MAT_FINAL_ASSEMBLY);


    /* counter0 = counter; */
    /* counter += systemSize*nTime; */
    /* bigAssSnapshotArray = yk_reallocPetscScalar(bigAssSnapshotArray,counter); */
    /* for (j=counter0; j<counter; j++) */
    /*   bigAssSnapshotArray[j] = snapshotArray[j-counter0]; */
  }
  free(snapshotArray);

  /* for (j=0; j<counter/systemSize; j++){ */
  /*   for (k=0; k<systemSize; k++) */
  /*     MatSetValue(*resSnapshot, k, j, bigAssSnapshotArray[j*systemSize+k], */
  /* 		  INSERT_VALUES); */
  /* } */
  MatAssemblyBegin(*resSnapshot, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*resSnapshot, MAT_FINAL_ASSEMBLY);
  return sResidual_mu;
}

void yk_fgets(FILE *stream, int row, double **snapshotArray, int *count){
  int i;
  int bufferSize = 0;
  int endOfLine = 0;
  char input[1000000] = {0};
  char *targetBuffer = (char *) malloc (sizeof(char));
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  *count = 0;
  targetBuffer[0] = '\0';
  *snapshotArray = (double *) malloc (sizeof(double));
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  while (fgets(input, sizeof(input), stream)){
    bufferSize += strlen(input); // new size
    targetBuffer = (char *)realloc(targetBuffer, (bufferSize+1)*sizeof(char));
    strcat(targetBuffer, input);
    if (strchr(input, '\n') != NULL){
      *snapshotArray=
	(double *)realloc(*snapshotArray,(*count+1)*row*sizeof(double));
      (*snapshotArray)[*count*row] = atof(strtok(targetBuffer, " "));

      for (i=1; i<row; i++){
	(*snapshotArray)[*count*row+i] = atof(strtok(NULL, " "));
      }

      bufferSize = 0;
      targetBuffer[0] = '\0';
      (*count)++;
    }
  }
  free(targetBuffer);
}


/*     //if (strchr(input, "\n") != NULL)


}


void ks_read_solution(Utype* U, int Umode, int node){
  int i;
  double DNU;
  FILE* solution_file;
  char outputFile[1000];

  U->node = node;
  if (Umode == 0){
    sprintf(outputFile, "%s_%d.dat", U->solution.u_state.id, node);
  }else if (Umode == 2){
    sprintf(outputFile, "%s_%d.dat", U->solution.psi_tan.loc_id, node);
  }else if (Umode == 3){
    sprintf(outputFile, "%s_%d.dat", U->solution.psi_lag.loc_id, node);
  }
  solution_file = fopen(outputFile, "r");
  for (i=0; i<U->solution.u_state.size; i++){
    if (Umode == 0){
      fscanf(solution_file, " %lf %lf\n", &DNU, &U->solution.u_state.array[i]);
    }else if (Umode == 2){
      fscanf(solution_file, " %lf %lf\n", &DNU, &U->solution.psi_tan.array[i]);
    }else if (Umode == 3){
      fscanf(solution_file, " %lf %lf\n", &DNU, &U->solution.psi_lag.array[i]);
    }
  }
  fclose(solution_file);
}
/*
void ks_read_solution(Utype* U, int node){
  int i;                             //initialization for iteration
  double DNU;               //temporary variable
  FILE* solution_file;               //File object for the file reader
  char outputFile[1000];
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  U->node = node;
  sprintf(outputFile, "%s_%d.dat", U->solution.u_state.id, node);
  solution_file = fopen(outputFile, "r");
  for (i=0; i<U->solution.u_state.size; i++){
    fscanf(solution_file, " %lf %lf\n", &DNU, &U->solution.u_state.array[i]);
  }
  fclose(solution_file);
}

void ks_read_psi_tan(Utype* U, int seg_num, int node){
  int i;                             //initialization for iteration
  double DNU;               //temporary variable
  FILE* solution_file;               //File object for the file reader
  char outputFile[1000];
  //----------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------
  U->node = node;
  sprintf(outputFile, "%s_%d_%d.dat", U->solution.psi_tan.id, seg_num, node);
  solution_file = fopen(outputFile, "r");
  for (i=0; i<U->solution.u_state.size; i++){
    fscanf(solution_file, " %lf %lf\n", &DNU, &U->solution.psi_tan.array[i]);
  }
  fclose(solution_file);
}
*/
