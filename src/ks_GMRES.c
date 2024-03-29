/*                                                                      
==============================================================================
Name        : GMRES specific for Least Squares Shadowing
Author      : Yukiko Shimizu                                             
Version     : March 22, 2015                             
Copyright   : Don't Touch. Hehe                                              
Description : Implementation of LSS Type III: Checkpoint design.          
"Offers a means to reduce the size of the KKT system. Instead    
 of searching for the entire shadow trajectory u(t) at once,      
 this searches for teh shadow trajectory at each checkpoint"      
 =============================================================================
*/
                                                                              
#include "ks_GMRES.h"
//Here for matvec do &name of the function
void ks_GMRES(yk_PrimalSolver *ykflow, Multiverse *multiquation,
	      Cluster *primal, double* R0, Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section
  //------------------------------------//-------------------------------------
  int i, j, n;                          //initialization for iteration
  int krylov_Count=0;
  int systemSize;
  int solverSize;
  int preIt;                            //size of p also; arbitrary really
  int out_i = 0;                        //outer and inner iteration numbers
  int in_i;
  double norm_residual;                 //2nd Norm of the Residual Vector  
  PetscInt hij_0, hij_1;
  PetscScalar hij_v0, hij_v1;
  double *c, *s, *p, *y;                //GMRES constants 
  double gamma;                         //distance between two Hessenberg point
  double tbs;                           //back substitution place holder
  double con = 12;                      //norm of the system residual (error) 
  double *axKrylov;   
  double *R;
  double *krylov_vector;
  char krylov_output[1000];             //Name of the krylov vectors 
  FILE* krylov;                         //Initialize File for krylov to output
  PetscInt *nz;                         //number of nonzero values in hess mat
  Mat hessenberg;                       //Hessenberg matrix to manipulate
  double *hRow;                         //Values at a row in Hessenberg
  double cpu_time_used;
  double cpu_time_used_total = 0;
  clock_t start, end;
  Universe _eqnOfInterest;
  Galaxy *_solOfInterest;
  FILE *krylov_info;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  if (!reduced->restart)
    krylov_info = fopen("krylov.info", "w+");
  else
    krylov_info = fopen("krylov.info", "a");
  if (reduced->reducedSolution == 1){                                        
    _eqnOfInterest = multiquation->equationReduced;
    _solOfInterest = primal->reduced;                                          
  }else{                                                                       
    _eqnOfInterest = multiquation->equation;           
    _solOfInterest = primal->self;                                             
  }
  systemSize = _eqnOfInterest.numStates*_solOfInterest->space.node.count;
  solverSize = systemSize*(2*reduced->nTimeSegs-1);
  axKrylov = (double *) malloc (solverSize*sizeof(double));           
  R = (double *) malloc (solverSize*sizeof(double));               
  krylov_vector = (double *) malloc (solverSize*sizeof(double)); 
  preIt = solverSize * 1.5;
  c = (double *) malloc (preIt*sizeof(double));
  s = (double *) malloc (preIt*sizeof(double));
  p = (double *) malloc ((preIt+1)*sizeof(double));
  y = (double *) malloc (preIt*sizeof(double));
  nz = (PetscInt *) malloc (preIt*sizeof(PetscInt));
  for (i=0; i<preIt; i++)           //designated upper triangle non zero
    nz[i] = i+2;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  do{
    MatCreateSeqAIJ(PETSC_COMM_SELF, preIt, preIt+1, 0, nz, &hessenberg);
    //while (con > pow(10, -6)){           //GMRES outer iterations
    in_i = 0;                          //Initialize inner iteration counter
    //-------------------------------------------------------------------------
    // Calculate initial Residuals R = Ax+b given the guess of R0
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");    
    printf("|             GMRES Solver, R= Ax+b                      |\n"); 
    printf("----------------------------------------------------------\n");  
    primal->self->beta = 1;          //beta = 1 is the R = Ax-b calculati
    ks_adjoint_MATVEC(ykflow, multiquation, primal, R0,  R, reduced);
    //-------------------------------------------------------------------------
    //Calculate the First Krylov Vector from the Output of the MATVEC algorithi
    //-------------------------------------------------------------------------
    for(i=0; i<preIt; i++)             //zero out the p vector
      p[i] = 0;                        //Empty the values
    norm_residual = sqrt(dot(R, R, solverSize)); //compute initial residual nor
    p[0] = norm_residual;              //Initialize right hand side
    sprintf(krylov_output, "%s_%d", "krylov", 0); //create title character
    krylov = fopen(krylov_output, "w");//print this krylov vector to file
    for (i=0; i<solverSize; i++)       //Krylov Vector
      fprintf(krylov, "%0.16f\n", R[i]/norm_residual);  //print to the output
    fclose(krylov);
    //-------------------------------------------------------------------------
    // Beginning of the Arnoldi Formulation
    //-------------------------------------------------------------------------
    do{
      start = clock();
      hRow = (double *) malloc ((in_i+2)*sizeof(double)); //Used to find hess
      //-----------------------------------------------------------------------
    // Matrix-Vector Product of Current Vector R = Ax
      //-----------------------------------------------------------------------
      sprintf(krylov_output, "%s_%d", "krylov", in_i);
      krylov = fopen(krylov_output,"r");
      for (i=0; i<solverSize; i++)
        fscanf(krylov, "%lf\n",  &krylov_vector[i]);
      fclose(krylov);
      primal->self->beta = 0;
      ks_adjoint_MATVEC(ykflow, multiquation, primal, krylov_vector, axKrylov,
			reduced);
      
      //----------------------------------------------------------------------
      // Gramm-Schmidt Orthogonalization
      //----------------------------------------------------------------------
      for (i=0; i<in_i+1; i++){
        sprintf(krylov_output, "%s_%d","krylov", i);
        krylov = fopen(krylov_output,"r");
        for (j=0; j<solverSize; j++)
          fscanf(krylov, "%lf\n",  &krylov_vector[j]);
        fclose(krylov);
        hRow[i] = dot(axKrylov, krylov_vector, solverSize);
        for (j=0; j<solverSize; j++)
          axKrylov[j] -= hRow[i]*krylov_vector[j];
      }
      hRow[in_i+1] = sqrt(dot(axKrylov, axKrylov, solverSize));
      //-----------------------------------------------------------------------
      // Define next Krylov Vector
      //----------------------------------------------------------------------
      sprintf(krylov_output, "%s_%d","krylov", in_i+1);
      krylov = fopen(krylov_output,"w");
      for (j=0; j<solverSize; j++)   //Find the First Krylov Vector
        fprintf(krylov, "%0.16f\n", axKrylov[j]/hRow[in_i+1]);
      fclose(krylov);
      //-----------------------------------------------------------------------
      // Previous Givens rotations on Hessenberg matrix
      //-----------------------------------------------------------------------
      for (i=0; i<in_i; i++){
        hij_v0 = hRow[i];   //Don't want value to be overwritten
        hRow[i] = c[i]*hij_v0 + s[i]*hRow[i+1];
        hRow[i+1] = -s[i]*hij_v0 + c[i]*hRow[i+1];
      }
      //Compute next rotation
      gamma = sqrt(hRow[in_i]*hRow[in_i] + hRow[in_i+1]*hRow[in_i+1]); //dot pro
      c[in_i] = hRow[in_i]/gamma;
      s[in_i] = hRow[in_i+1]/gamma;
      hRow[in_i] = gamma;
      hRow[in_i+1] = 0;
      for (i=0; i<in_i+1; i++)
        MatSetValue(hessenberg, in_i, i, hRow[i], INSERT_VALUES);
      p[in_i+1] = -s[in_i]*p[in_i];       //Given rotation on p
      p[in_i] = c[in_i]*p[in_i];          //Edits the p (rhs) matrix
      con = fabs(p[in_i+1]);
      in_i++;                             //new inner iteration
      free(hRow);
      end = clock();
      cpu_time_used = ((double) (end-start))/CLOCKS_PER_SEC;
      cpu_time_used_total += cpu_time_used;
      printf("----------------------------------------------------------\n");
      printf("|            GMRES Solver, Krylov Iteration %d            |\n",
    	     krylov_Count);
      printf("----------------------------------------------------------\n");
      fprintf(krylov_info, "%d. %s %0.16f %s %0.16f\n", krylov_Count,
    	      "Krylov error", con, "Krylov time in s", cpu_time_used);
      printf("κ_ε = %0.16e, t(s) = %0.16f\n", con, cpu_time_used);
      krylov_Count++;
      if (in_i == reduced->innerStart && reduced->innerStart > 0){
    	printf("why is this so complciated\n");
	break;
      }
    }while (con > pow(10,-6));
    MatAssemblyBegin(hessenberg, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(hessenberg, MAT_FINAL_ASSEMBLY);
    //-------------------------------------------------------------------------
    // Solve the hessenberg matrix using back substitution since LU
    //-------------------------------------------------------------------------
    for (i=0; i<in_i; i++){
      tbs = 0;
      hij_0 = in_i-i-1;
      MatGetValues(hessenberg, 1, &hij_0, 1, &hij_0, &hij_v0);
      for (j=0; j<i; j++){
	hij_1 = in_i-i+j;
	MatGetValues(hessenberg, 1, &hij_1, 1, &hij_0, &hij_v1);
	tbs += -hij_v1*y[in_i-i+j];
      }
      y[in_i-i-1] = (p[in_i-i-1]-tbs)/-hij_v0;
    }
    //-------------------------------------------------------------------------
    // Update and Form Approximate Solution and Recalcute Residuals
    //-------------------------------------------------------------------------
    for (i=0; i<in_i; i++){
      sprintf(krylov_output, "%s_%d","krylov", i);
      krylov = fopen(krylov_output,"r");
      for (j=0; j<solverSize; j++){
	fscanf(krylov, "%lf\n",  &krylov_vector[j]);
	R0[j] += y[i]*krylov_vector[j];
      }
      fclose(krylov);
    }
    MatDestroy(&hessenberg);
    out_i++;
    yk_printLSSInitialCon(solverSize, R0, out_i);
    fprintf(krylov_info, "%s %0.16f %s %0.16f\n", "Final Krylov error", con,
	    "t(s) ", cpu_time_used_total);
    printf("Final Krylov error %g\n", con);
    if (out_i == reduced->outerStart && reduced->outerStart > 0){
      break;
    }
  } while (con > pow(10, -6));
  fclose(krylov_info);
  free(nz);
  free(c);
  free(s);
  free(p);
  free(y);
  free(R);
  free(krylov_vector);
  free(axKrylov);
  MatDestroy(&hessenberg);
}

