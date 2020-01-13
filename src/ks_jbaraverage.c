#include "ks_jbaraverage.h"
//-----------------------------------------------------------------------------
// Finds the Output Time Average of J // Easy program
//
// Yukiko Shimizu
// April 23, 2016
// Error Estimation for Chaotic Systems
//-----------------------------------------------------------------------------
void ks_jbaraverage(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Galaxy *primal, Is_it *reduced){
  int i, j, k, n;  //initialization for iteration
  FILE* average_file;
  char outfile[10000];
  double average_spatial[primal->time.count+1];
  double instantaneous_dJ;
  double tempValue = 0;
  PetscScalar temp_dJdU;
  PetscInt vec_index;
  primal->j_bar = 0;
  //  VecCreateSeq(PETSC_COMM_SELF, primal->space.node.count, &primal->djdu);
  //---------------------------------------------------------------------------
  // Implementation starts here
  //---------------------------------------------------------------------------
  for (i=0; i<primal->time.count+1; i++) //Initialize here
    average_spatial[i] = 0;
  //---------------------------------------------------------------------------
  // Now calulate the integral of x at all the other time points
  //---------------------------------------------------------------------------
  for (i=0; i<primal->time.count+1; i++){
    ks_readSolution(multiquation->equation, primal, i);
    spatialOutputAverage(ykflow, multiquation, primal, &average_spatial[i], i);
  }
  //---------------------------------------------------------------------------
  // Jbaravaergae here
  //---------------------------------------------------------------------------
  sprintf(outfile, "%s_%s.dat", primal->id, "jbar");
  average_file = fopen(outfile, "w");
  for (i=0; i<primal->time.count; i++){
    primal->j_bar += (primal->time.dt/2.0)*(1.0/(primal->time.globalT_f))
      *(average_spatial[i]+average_spatial[i+1]);
    tempValue += (primal->time.dt/2.0)*
      (average_spatial[i]+average_spatial[i+1]);
    instantaneous_dJ = tempValue*(1.0/(primal->time.dt*(i+1)));
    fprintf(average_file,"%0.16f %0.16f\n", primal->time.dt*(i+1),
	    instantaneous_dJ);
  }
  fclose(average_file);
  //---------------------------------------------------------------------------
  // djdu calculation
  //---------------------------------------------------------------------------
  /* VecSet(primal->djdu, 0); */
  /* if(primal->space.node.count > 1){ */
  /*   for (i=0; i<primal->space.elem.count; i++){ */
  /*     for (j=0; j<multiquation->equation.numStates; j++){ */
  /*       for (k=0; k<primal->basis.p+1; k++){ */
  /*         temp_dJdU = 0; */
  /*         vec_index = (i*multiquation->equation.numStates+j)* */
  /* 	    (primal->basis.p+1)+k; */
  /*         for (n =0; n<primal->quad.n; n++) */
  /*           temp_dJdU += (primal->space.dx* */
  /* 			  phi(primal->basis.p, k, primal->quad.x_i.array[n])* */
  /* 			  primal->quad.w.array[n])/(2*primal->space.x_f); */
  /*         VecSetValue(primal->djdu, vec_index, temp_dJdU, INSERT_VALUES); */
  /*       } */
  /*     } */
  /*   } */
  /* }//addd other things later */
  /* VecAssemblyBegin(primal->djdu); */
  /* VecAssemblyEnd(primal->djdu);        */
}


void spatialOutputAverage(yk_PrimalSolver *ykflow, Multiverse *multiquation,
			  Galaxy *primal, double *aveSpace, int node){
  int i, j;
  double state;
  //This works when there's a space. It Lorenz, you need to edit this more
  //---------------------------------------------------------------------------
  // Initialization
  //--------------------------------------------------------------------------
  *aveSpace = 0;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Average output for a particular time t in space
  //---------------------------------------------------------------------------
  ks_readSolution(multiquation->equation, primal, node);  
  for (i=0; i<primal->space.elem.count; i++){                         
    for (j=0; j<primal->quad.n; j++){                                
      state = U_xi(primal, 0, i, primal->quad.x_i.array[j]);
      *aveSpace += (primal->space.dx/2.0)*state*
	primal->quad.w.array[j]/primal->space.x_f; 
    }                                                                      
  }
}

void ks_dJdU(yk_PrimalSolver *ykflow, Multiverse *multiquation, Galaxy *primal,
	     int time_i, Vec dJdU, Is_it *reduced){
  int i, j, k, n;
  double temp_dJdU;
  PetscInt vec_index;
  Universe _equation = multiquation->equation;
  VecZeroEntries(dJdU); //Make sure to initialize dJdU to all zeros first      
  if(primal->space.node.count > 1){
    for (i=0; i<primal->space.elem.count; i++){
      for (j=0; j<_equation.numStates; j++){
        for (k=0; k<primal->basis.p+1; k++){
          temp_dJdU = 0;
          vec_index = (i*_equation.numStates+j)*
        (primal->basis.p+1)+k;
          for (n =0; n<primal->quad.n; n++)
            temp_dJdU += (primal->space.dx*
              phi(primal->basis.p, k, primal->quad.x_i.array[n])*
              primal->quad.w.array[n])/(2*primal->space.x_f);
          VecSetValue(dJdU, vec_index, temp_dJdU, INSERT_VALUES);
        }
      }
    }
  }
  VecAssemblyBegin(dJdU); //Confused if this should be done once or multiple   
  VecAssemblyEnd(dJdU);   //times                                    
}

void ks_calculatedjdu(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Galaxy *primal, Vec djdu, int timeNode, Is_it *reduced){
  Mesh _meshOfInterest; 
  int systemSize;
  Vec djduT;
  Universe _eqnOfInterest = multiquation->equation; 
  /* if (reduced->reducedSolution == 1){ */
  /*   _meshOfInterest = reduced->reducedMesh; */
  /* }else if (reduced->reducedSolution == 0){ */
  /*   _meshOfInterest = primal->space; */
  /* } */
  systemSize = primal->space.node.count*_eqnOfInterest.numStates;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &djduT);
  ykflow->spatialdJdU(ykflow, multiquation, primal, timeNode, djduT, reduced);
  if (reduced->reducedSolution == 1)
    MatMultTranspose(reduced->rOBState, djduT, djdu);
  else
    VecCopy(djduT, djdu);
  VecDestroy(&djduT);
}


/* void fourierTransformApprox(Universe equation, Utype *primal){ */
/*   //where this function is will act upon the time average outpute J */
/*   //going to approximate the fourier transform just using trapezoidal rule */
/*   int i, j;  //initialization for iteration */
/*   FILE* average_file; */
/*   char outfile[10000]; */
/*   double fourierReal; */
/*   double fourierImg; */
/*   double amplitude[primal->time.count/2]; */
/*   double f[primal->time.count/2]; */
/*   double complex exponential; */
/*   double real; */
/*   double img; */
/*   double frequency = 1/primal->time.t_f; */
/*   double *timeAverageOutput */
/*     = (double *) malloc ((primal->time.count+1)*sizeof(double)); */
/*   double *timet = (double *) malloc ((primal->time.count+1)*sizeof(double)); */
/*   //---------------------------------------------------------------------------- */
/*   // Implementation */
/*   //---------------------------------------------------------------------------- */
/*   sprintf(outfile, "%s_%s.dat", primal->id, "jbar"); */
/*   printf("%s\n", outfile); */
/*   getchar(); */
/*   average_file = fopen(outfile, "r");   */
/*   for (i=0; i<primal->time.count; i++) */
/*     fscanf(average_file, " %lf %lf\n", &timet[i], &timeAverageOutput[i]);  */
/*   fclose(average_file); */

/*   sprintf(outfile, "%s_%s.dat", primal->id, "fourierTransform"); */
/*   average_file = fopen(outfile, "w"); */
/*   for (i=0; i<primal->time.count/2; i++){ */
/*     fourierReal = 0; */
/*     fourierImg = 0; */
/*     for (j=0; j<primal->time.count; j++){ */
/*       real = cos(2.0*M_PI*i*j/primal->time.count); */
/*       img = -sin(2.0*M_PI*i*j/primal->time.count); */
/*       fourierReal += timeAverageOutput[j]*real; */
/*       fourierImg += timeAverageOutput[j]*img; */
/*     } */
/*     f[i] = i*frequency; */
/*     amplitude[i] = sqrt(fourierReal*fourierReal+fourierImg*fourierImg); */
/*     fprintf(average_file,"%0.16f %0.16f\n", f[i], amplitude[i]);  */
/*   } */
/*   fclose(average_file); */
/* }    */
 
