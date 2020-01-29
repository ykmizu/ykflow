
/*
  =============================================================================
  Name        : Least Squares Shadowing Type III: Checkpoint Design
  Author      : Yukiko Shimizu
  Version     : September 2, 2015
  Description : Control Panel for Least Squares Shadowing Type III
                for the Lorenz Attractor. (LSS TYPE II)
                See the following paper for math:
  Reference   : Least Squares Shadowing for Sensitivity Analysis
                of Turbulent Fluid Flows 16 Jan 2014
  =============================================================================
*/

#include "read_file.h"

Universe read_eqn(int argc, char* argv[]){
  int i;                            //initialization for iteration
  int num_length = 1024;            //define buffer length
  FILE *input_file;                 //input file name
  char buffer[num_length];          //define the buffer
  char *current_line;               //initialize name for line in file
  char *lines = NULL;
  char *lines_temp;
  char *params_line;
  size_t length;
  size_t length_temp;
  Universe eqn;// = (Universe *) malloc (sizeof(Universe));
  memset(buffer, ' ', sizeof(buffer) -1);
  buffer[sizeof(buffer) - 1] = '\0';
  input_file = fopen(argv[1], "r");  //open file
  //Check first to see if the file was successfully opened
  if (input_file == NULL){
    perror("Error");
    exit( EXIT_FAILURE);  //If not exit
  }else{
    while ((current_line = fgets(buffer, sizeof(buffer), input_file))){
      if (current_line[0]!= '\n' && current_line[0]!='%'){
        lines = strtok(current_line, "=");
        if (lines == NULL){
          perror("Error");
          exit( EXIT_FAILURE );
        }
        length = strlen(lines);
        if (lines[length-1] == ' '){
          lines[length-1] = '\0';
        }
        lines_temp = strtok(NULL, "=");
        if (lines_temp == NULL){
          perror("Error");
          exit( EXIT_FAILURE );
        }
        length_temp = strlen(lines_temp);
        if (lines_temp[0] == ' '){
          memmove(lines_temp, lines_temp +1, length_temp);
        }else{
          lines_temp[length_temp-1] = '\0';
        }
        if (strcmp(lines, "eqnName") == 0){
          length = strlen(lines_temp);
          if (lines_temp[length-1] == '\n')
            lines_temp[length-1] = '\0';
          strcpy(eqn.nameEqn, lines_temp);
        }else if(strcmp(lines, "numStates") == 0){
          eqn.numStates = atoi(lines_temp);
        }else if(strcmp(lines, "numParams") == 0){
          eqn.numParams = atoi(lines_temp);
          eqn.c = (double *) malloc (eqn.numParams*sizeof(double));
        }else if(strcmp(lines, "params") == 0){
          params_line = strtok(lines_temp, " ");
          for (i=0; i<eqn.numParams; i++){
            eqn.c[i] = atof(params_line);
            params_line = strtok(NULL, " ");
          }
        }else if(strcmp(lines, "paramIndex") == 0){
	  eqn.paramIndex = atoi(lines_temp);
	}else if(strcmp(lines, "paramNewValue") == 0){
	  eqn.paramValue = atof(lines_temp);
	}
      }
    }

  }
  fclose(input_file);
  return eqn;
}

void read_input(int argc, char *argv[], Universe eqn, Galaxy* u_fine,
                Galaxy* u_coarse, Is_it *reduced){
  int i;
  int inter;
  int num_length = 1024;               //define buffer length
  FILE *input_file;                    //input file name
  char buffer[num_length];             //define the bugger
  char *current_line;                  //initialize name for line in file
  char *lines = NULL;
  char *lines_temp;
  size_t length;
  size_t length_temp;
  char outputfile[1000];
  memset(buffer, ' ', sizeof(buffer) -1);
  buffer[sizeof(buffer) - 1] = '\0';
  //  if (argc > 2)
  //sprintf(outputfile, "%s", argv[2]);
  //else
  sprintf(outputfile, "%s.inp", eqn.nameEqn);
  input_file = fopen(outputfile, "r");  //open file
  //Check first to see if the file was successfully opened
  if (input_file == NULL){
    perror("Error");
    exit( EXIT_FAILURE);  //If not exit
  }else{
    while ((current_line = fgets(buffer, sizeof(buffer), input_file))){
      if (current_line[0]!= '\n' && current_line[0]!='%'){
        lines = strtok(current_line, "=");
        if (lines == NULL){
          perror("Error");
          exit( EXIT_FAILURE );
        }
        length = strlen(lines);
        if (lines[length-1] == ' '){
          lines[length-1] = '\0';
        }
        lines_temp = strtok(NULL, "=");
        if (lines_temp == NULL){
          perror("Error");
          exit( EXIT_FAILURE );
        }
        length_temp = strlen(lines_temp);
        if (lines_temp[0] == ' '){
          memmove(lines_temp, lines_temp +1, length_temp);
        }else{
          lines_temp[length_temp-1] = '\0';
        }
        if (strcmp(lines, "timeMethod") == 0){
          u_coarse->time_method = atoi(lines_temp);
        }else if (strcmp(lines, "burnTimeLength") == 0){
          u_coarse->burnT_f = atof(lines_temp);
        }else if (strcmp(lines, "dt") == 0){
          u_coarse->time.dt = atof(lines_temp);
        }else if (strcmp(lines, "T") == 0){
          u_coarse->time.globalT_f = atof(lines_temp);
          u_coarse->time.t_f = atof(lines_temp);
        }else if (strcmp(lines, "T_0") == 0){
          u_coarse->time.globalT_0 = atof(lines_temp);
          u_coarse->time.t_0 = atof(lines_temp);
        }else if (strcmp(lines, "dx") == 0){
          u_coarse->space.dx = atof(lines_temp);
        }else if (strcmp(lines, "X") == 0){
          u_coarse->space.x_f = atof(lines_temp);
        }else if (strcmp(lines, "X_0") == 0){
          u_coarse->space.x_0 = atof(lines_temp);
        }else if (strcmp(lines, "interpolationNumber") == 0 ){
          u_coarse->basis.p = atoi(lines_temp);
        }else if (strcmp(lines, "newInterpolationNumber") == 0){
          u_fine->basis.p = atoi(lines_temp);
        }else if (strcmp(lines, "newTimeMethod") == 0){
          u_fine->time_method = atoi(lines_temp);
        }else if (strcmp(lines, "numTimeSegments") == 0){
          reduced->nTimeSegs = atoi(lines_temp);
	  reduced->win_i = (PerWindow *) malloc
	    (reduced->nTimeSegs*sizeof(PerWindow));
	}else if (strcmp(lines, "numSubSpacesPerWindow")==0){
	  reduced->nSubWindows = atoi(lines_temp);
	}else if (strcmp(lines, "numSubSpacesPerWindow") == 0){
	  reduced->nSubPWin = atoi(lines_temp);
        }else if (strcmp(lines, "dss") == 0){
          reduced->dss = atoi(lines_temp);
        }else if (strcmp(lines, "numBasisFunctions") == 0){
          reduced->nBasisFuncs = atoi(lines_temp);
	}else if (strcmp(lines, "numBasisTime") == 0){
	  reduced->nBasisTime = atoi(lines_temp);
	}else if (strcmp(lines, "engyBasisFunctions") == 0){
	  reduced->eBasisSpace = atof(lines_temp);
	}else if (strcmp(lines, "engyBasisTime") == 0){
	  reduced->eBasisTime = atof(lines_temp);
        }else if (strcmp(lines, "numSampleNodes") == 0){
          reduced->nSampleNodes = atoi(lines_temp);
        }else if (strcmp(lines, "numBasisFunctionsRJ") == 0){
          reduced->nBasisFuncsRJ = atoi(lines_temp);
        }else if (strcmp(lines, "reducedLSS") == 0){
	  reduced->reducedLSS = atoi(lines_temp);
	}else if (strcmp(lines, "restart") == 0){
	  reduced->restart = atoi(lines_temp);
	}else if (strcmp(lines, "innerGMRESIterationNumber") == 0){
	  reduced->innerStart = atoi(lines_temp);
	}else if (strcmp(lines, "outerGMRESIterationNumber") == 0){
	  reduced->outerStart = atoi(lines_temp);
	}else if (strcmp(lines, "numParams") == 0){
	  reduced->numParams = atoi(lines_temp);
	  reduced->params = (double *) malloc
	    (reduced->numParams*sizeof(double));
	  reduced->paramsL = (double *) malloc
	    (reduced->numParams*sizeof(double));
	  reduced->paramsH = (double *) malloc
	    (reduced->numParams*sizeof(double));
	  reduced->dparams = (double *) malloc
	    (reduced->numParams*sizeof(double));
	}else if (strcmp(lines, "params") == 0){
	  reduced->params[0] = atof(strtok(lines_temp, " "));
	  for (i=1; i<reduced->numParams; i++)
	    reduced->params[i] = atof(strtok(NULL, " "));
	}else if (strcmp(lines, "paramsL") == 0){
	  reduced->paramsL[0] = atof(strtok(lines_temp, " "));
	  for (i=1; i<reduced->numParams; i++)
	    reduced->paramsL[i] = atof(strtok(NULL, " "));
	}else if (strcmp(lines, "paramsH") == 0){
	  reduced->paramsH[0] = atof(strtok(lines_temp, " " ));
	  for (i=1; i<reduced->numParams; i++)
	    reduced->paramsH[i] = atof(strtok(NULL, " "));
	}else if (strcmp(lines, "dparams")== 0 ){
	  reduced->dparams[0] = atof(strtok(lines_temp, " "));
	  for (i=1; i<reduced->numParams; i++)
	    reduced->dparams[i] = atof(strtok(NULL, " "));
	  reduced->numParamSet = 1;
	  for (i=0; i<reduced->numParams; i++){
	    inter= ((reduced->paramsH[i]*10-reduced->paramsL[i]*10)/
		    reduced->dparams[i])/10+1;
	    reduced->numParamSet *= inter;
	  }
	  /* }else if (strcmp(lines, "typeSolutionLSS") ==0){ */
	  /*   miscParameters[6] = atoi(lines_temp); */
	  //}else if (strcmp(lines, "designParameter") ==0){
	  //miscParameters[7] = atoi(lines_temp);
	  //	}else if (strcmp(lines, "newDesignParameter") == 0){
  	  //	  miscParameters[8] = atoi(lines_temp);
	}
      }
    }
  }
  fclose(input_file);
}
