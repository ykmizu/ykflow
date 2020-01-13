
#include "ks_dRdU.h"
//-----------------------------------------------------------------------------
// Implementation Of dRdU Via Finite Diferencing
//
// Yukiko Shimizu
// April 23, 2016
// Error Estimation for Chaotic Systems
//-----------------------------------------------------------------------------
void ks_totaldRdU(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		  Cluster *object, Mat dRdU, Is_it *reduced,
		  int innertimesolver){
  int i;
  int systemSize;
  int reducedSize = 0;
  Universe _eqnOfInterest;
  Mat Mij;                        //mass matrix
  double * c = getBDFCoef(innertimesolver);  
  PetscScalar drdu_coef = c[0]/object->self->time.dt;   //Coefficients for BDF
  //---------------------------------------------------------------------------
  // Initialization and Mallocing/ Creating shit
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1)                                          
    _eqnOfInterest = multiquation->equationReduced;
  else if (reduced->reducedSolution == 0 &&                                    
           strncmp(object->self->utypeId, "psi", 3)==0)                       
    _eqnOfInterest =multiquation->equationLSS;                                 
  else                                                                         
    _eqnOfInterest  = multiquation->equation;    
  if (reduced->hrom == 1 && reduced->reducedSolution == 0)                
    reducedSize = reduced->reducedMesh.node.count*_eqnOfInterest.numStates;  
  else if (reduced->hrom == 1 && reduced->reducedSolution == 1)            
    reducedSize = object->self->space.node.count*_eqnOfInterest.numStates;  
  else if (reduced->hrom == 0)                                         
    reducedSize = object->self->space.node.count*_eqnOfInterest.numStates; 
  systemSize =  object->self->space.node.count*_eqnOfInterest.numStates;
  //MatCreateSeqDense(PETSC_COMM_SELF, reducedSize, systemSize, NULL, &Mij);
  MatCreateSeqAIJ(PETSC_COMM_SELF, reducedSize, systemSize,
		  object->self->basis.p+1, PETSC_NULL, &Mij);
  //---------------------------------------------------------------------------
  // Mass Matrix
  //---------------------------------------------------------------------------
  ks_mass(_eqnOfInterest, object->self, Mij, reduced); //Don't need the ide
  ks_spatialdRdU(ykflow, multiquation, object, dRdU, reduced);
  if (reduced->reducedSolution == 1){
    for (i=0; i<systemSize; i++)
      MatSetValue(dRdU, i, i, drdu_coef, ADD_VALUES);
    MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY);
  }else{
    MatAXPY(dRdU, drdu_coef, Mij, SUBSET_NONZERO_PATTERN);
  }
  MatDestroy(&Mij);
  free(c);
}            

//newer version of the ks_spatialdRdU script
void ks_spatialdRdU(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		    Cluster *object, Mat dRdU, Is_it *reduced){
  Universe _equation = multiquation->equation;
  MatZeroEntries(dRdU);
  //---------------------------------------------------------------------------
  // Malloc and Create Vectors and or Matrices here
  //---------------------------------------------------------------------------
  if (strcmp(object->self->utypeId, "state")==0)
    ks_spatialStatedRdU(_equation, object->self, dRdU, reduced);
  else if (strncmp(object->self->utypeId, "psi2", 4) == 0)
    ks_dres_tangent(ykflow, multiquation, object, dRdU, reduced);
  else if (strncmp(object->self->utypeId, "psi1", 4) == 0)
    ks_dres_lagrange(ykflow, multiquation, object, dRdU, reduced);
  MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY);
}

void ks_spatialStatedRdU(Universe equation, Galaxy *object, Mat dRdU,
                         Is_it *reduced){
  if (strncmp(equation.nameEqn, "meks", 4)==0){
    ks_dres_advection(equation, object, dRdU, reduced);
    ks_dres_burgers(equation, object, dRdU, reduced);
    ks_dres_diffusion(equation, object, dRdU, reduced);
    ks_dres_fourth_IPDG(equation, object, dRdU, reduced);
  }
}



/* void ks_spatialdRdUCheck(Multiverse multiquation, Cluster *object, */
/* 			 Is_it *reduced){ */
/*   int i, j, k; */
/*   int alfa = -1; */
/*   int systemSize; */
/*   int col_in; */
/*   int *index; */
/*   double * vals; */
/*   PetscScalar smallstep = pow(10, -8); */
/*   PetscScalar scale = 1.0/smallstep; */
/*   Mat dRdU; */
/*   Vec k1, k2; */
/*   Universe _eqnOfInterest; */
/*   //--------------------------------------------------------------------------- */
/*   // Initialization  */
/*   //--------------------------------------------------------------------------- */
/*   if (reduced->reducedSolution == 1)                                        */
/*     _eqnOfInterest = multiquation.equationReduced;                        */
/*   else if (reduced->reducedSolution == 0 &&                              */
/*            strncmp(object->self->utypeId, "psi", 3)==0)                     */
/*     _eqnOfInterest = multiquation.equationLSS;                              */
/*   else                                                                     */
/*     _eqnOfInterest  = multiquation.equation;     */
/*   systemSize = _eqnOfInterest.numStates*object->self->space.node.count; */
/*   index = (int *) malloc (systemSize*sizeof(int)); */
/*   for (i=0; i<systemSize; i++) */
/*     index[i] = i; */
/*   vals = (double *) malloc (systemSize*sizeof(double)); */
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &k1);                      */
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &k2);   */
/*   MatCreateSeqDense(PETSC_COMM_SELF, systemSize, systemSize, NULL, &dRdU); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation */
/*   //--------------------------------------------------------------------------- */
/*   ks_spatialResidual(multiquation, object, k1, reduced); */
/*   for (i=0; i<object->self->space.elem.count; i++){                        */
/*     for (j=0; j< _eqnOfInterest.numStates; j++){                           */
/*       for (k=0; k<object->self->basis.p+1; k++){                           */
/*         col_in = (i*_eqnOfInterest.numStates+j)*                           */
/* 	  (object->self->basis.p+1)+k;                                        */
/*         object->self->solution[j].array[i*(object->self->basis.p+1)+k] +=  */
/* 	  smallstep;                                                          */
/*         ks_spatialResidual(multiquation, object, k2, reduced); //Calculate  */
/* 	object->self->solution[j].array[i*(object->self->basis.p+1)+k] -=  */
/* 	  smallstep;       */
/* 	VecAXPY(k2, alfa, k1); //Calculates k2-k1   */
/* 	VecScale(k2, scale);   //Calculates (k2-k1)/smallStep  */
/* 	VecGetValues(k2, systemSize, index, vals); //Retrieves nonzero value  */
/* 	//  Now put those nonzero values into the Matrix dRdU */
/* 	MatSetValues(dRdU, systemSize, index, 1, &col_in, vals, INSERT_VALUES); */
/*       } */
/*     } */
/*   } */
/*   MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY); */
/*   MatView(dRdU, PETSC_VIEWER_STDOUT_SELF); */
/*   MatDestroy(&dRdU); */
/*   VecDestroy(&k1); */
/*   VecDestroy(&k2); */
/* } */

/* void ks_spatialdRdU(Multiverse multiquation, Cluster *object, Mat dRdU, */
/*  		    Is_it *reduced){ */
/*    int i, j, k, n;        //initialization for iteration */
/*   int systemSize; */
/*   int count;             //Counter to keep track of indices for row_in */
/*   int jstart, jend;  //gives starting/terminating of nonzero values along the o */
/*   PetscInt row_non_0;    //Number of nonzero values */
/*   PetscInt *row_in;      //Array of row indicies fo one colum */
/*   PetscInt col_in;       // current index value along the entire state vector */
/*   PetscScalar *vals;     // individual values for dRdU used to fill Mat */
/*   PetscScalar alfa = -1; //Used to find difference between two residuals */
/*   PetscScalar smallstep = pow(10, -5); */
/*   PetscScalar scale = 1.0/smallstep; //Used to divide difference by */
/*   Vec k1, k2;            //Create vectors that are used to calculate dRdU */
/*   Universe _eqnOfInterest; */
/*   //--------------------------------------------------------------------------- */
/*   // Malloc and Create Vectors and or Matrices here */
/*   //--------------------------------------------------------------------------- */
/*   if (object->self->space.node.count == 1 )//Firhas diff num of va */
/*     row_non_0 = multiquation.equation.numStates; */
/*   else */
/*     row_non_0 = 3*multiquation.equation.numStates*(object->self->basis.p+1); */
/*   PetscMalloc1(row_non_0, &row_in); */
/*   PetscMalloc1(row_non_0, &vals); */

/*   if (reduced->reducedSolution == 1) */
/*     _eqnOfInterest = multiquation.equationReduced; */
/*   else if (reduced->reducedSolution == 0 && */
/*            strncmp(object->self->utypeId, "psi", 3)==0) */
/*     _eqnOfInterest = multiquation.equationLSS; */
/*   else */
/*     _eqnOfInterest  = multiquation.equation; */

/*   systemSize =  object->self->space.node.count*_eqnOfInterest.numStates; */
  
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &k1); */
/*   VecCreateSeq(PETSC_COMM_SELF, systemSize, &k2); */
/*   //--------------------------------------------------------------------------- */
/*   // Implementation starts here */
/*   //--------------------------------------------------------------------------- */
/*   ks_spatialResidual(multiquation, object, k1, reduced); //Calculate original */
/*   //--------------------------------------------------------------------------- */
/*   // Begin dRdU calculation */
/*   //--------------------------------------------------------------------------- */
/*   for (i=0; i<object->self->space.elem.count; i++){ */
/*     for (j=0; j< _eqnOfInterest.numStates; j++){ */
/*       for (k=0; k<object->self->basis.p+1; k++){ */
/*         col_in = (i*_eqnOfInterest.numStates+j)* */
/* 	  (object->self->basis.p+1)+k; */
/*         count = 0; //Everytime you go to a new column, you need to zero it out */
/*         object->self->solution[j].array[i*(object->self->basis.p+1)+k] += */
/* 	  smallstep; */
/*         ks_spatialResidual(multiquation, object, k2, reduced); //Calculate */
/*         object->self->solution[j].array[i*(object->self->basis.p+1)+k] -= */
/* 	  smallstep; */
/*         //first element, but only looking at a point */
/*         if (i==object->self->space.elem.count-1 && */
/* 	    object->self->space.node.count > 1){ */
/*           for (n=0; n<_eqnOfInterest.numStates* */
/* 		 (object->self->basis.p+1); n++){ */
/*             row_in[count] = n; */
/*             count++; */
/*           } */
/*         } */
/*         if (i!=0){ */
/*           for (n=0; n<_eqnOfInterest.numStates* */
/* 		 (object->self->basis.p+1); n++){ */
/*             row_in[count] = (i-1)*(object->self->basis.p+1)+n; */
/*             count++; */
/*           } */
/*         } */
/*         for (n=0; n<_eqnOfInterest.numStates* */
/* 	       (object->self->basis.p+1); n++){ */
/*           row_in[count] = i*(object->self->basis.p+1)+n; */
/*           count++; */
/*         } */
/*         if (i!=object->self->space.elem.count-1){ */
/*           for (n=0; n<_eqnOfInterest.numStates* */
/* 		 (object->self->basis.p+1); n++){ */
/*             row_in[count] = (i+1)*(object->self->basis.p+1)+n; */
/*             count++; */
/*           } */
/*         } */
/*         if (i==0 && object->self->space.node.count > 1){ */
/*           for (n=0; n<_eqnOfInterest.numStates* */
/* 		 (object->self->basis.p+1); n++){ */
/*             row_in[count] = (object->self->space.elem.count-1)* */
/* 	      (object->self->basis.p+1)+n; */
/*             count++; */
/*           } */
/*         } */
/*         VecAXPY(k2, alfa, k1); //Calculates k2-k1 */
/*         VecScale(k2, scale);   //Calculates (k2-k1)/smallStep */
/* 	if (strncmp(object->self->utypeId, "psi2", 4) == 0){ */
/* 	  VecView(k2, PETSC_VIEWER_STDOUT_SELF); */
/* 	  getchar(); */
/* 	} */
/*         VecGetValues(k2, row_non_0, row_in, vals); //Retrieves nonzero values */
/*         //  Now put those nonzero values into the Matrix dRdU */
/*         MatSetValues(dRdU, row_non_0, row_in, 1, &col_in, vals, INSERT_VALUES); */
/*       } */
/*     } */
/*   } */
/*   //--------------------------------------------------------------------------- */
/*   // Assembly Those Beautiful Matrices */
/*   //--------------------------------------------------------------------------- */
/*   MatAssemblyBegin(dRdU, MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(dRdU, MAT_FINAL_ASSEMBLY); */
/*   //--------------------------------------------------------------------------- */
/*   // Free or Destroy everything */
/*   //--------------------------------------------------------------------------- */
/*   PetscFree(row_in); */
/*   PetscFree(vals); */
/*   VecDestroy(&k1); */
/*   VecDestroy(&k2); */
/* } */
