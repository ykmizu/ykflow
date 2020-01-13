#include "ks_mass.h"

void ks_mass(Universe equation, Galaxy* primal, Mat Mij, Is_it *reduced){
  int i, j, k, n, m;  //initialization for iteration
  int ij;       //index value combination of element number and basis value 
  PetscScalar *vals;
  PetscInt *rows, *cols; //indices refering to their respectful index
  int j_epsilon;
  int *getMeshIndex = createMeshMap(&primal->space);
  Mesh _meshOfInterest;
  PetscMalloc1(pow(primal->basis.p+1,2), &vals);
  PetscMalloc1(pow(primal->basis.p+1,2), &rows);
  PetscMalloc1(pow(primal->basis.p+1,2), &cols);
  if (reduced->hrom==1 && reduced->reducedSolution == 0)               
    _meshOfInterest = reduced->reducedMesh;                                    
  else if (reduced->hrom == 1 && reduced->reducedSolution == 1)      
    _meshOfInterest = primal->space;
  else if (reduced->hrom == 0)
    _meshOfInterest = primal->space;  
  //------------------------------------------------------------------------
  // Implementation
  //------------------------------------------------------------------------
  for (i=0; i<primal->basis.p+1; i++){
    for (j=0; j<primal->basis.p+1; j++)
      vals[j] = 0;
    for (j=0; j<primal->basis.p+1; j++){
      for (n=0; n<primal->quad.n; n++){
        if (primal->space.dx > 0){
          vals[j] += primal->space.dx/2.0*
            phi(primal->basis.p, i, primal->quad.x_i.array[n])*
            phi(primal->basis.p, j, primal->quad.x_i.array[n])*
            primal->quad.w.array[n]; 
        }
        else{
          vals[j] = 1;
        }
      }
    }
    for (j=0; j<_meshOfInterest.elem.count; j++){
      for (m =0; m<equation.numStates; m++){
        ij = (j*equation.numStates+m)*(primal->basis.p+1)+i; //same element 
        for (k=0; k<primal->basis.p+1; k++){
          j_epsilon = getMeshIndex[_meshOfInterest.elem.array[j]];
          cols[k] = (j_epsilon*equation.numStates+m)*(primal->basis.p+1)+k;
        }
        MatSetValues(Mij, 1, &ij, primal->basis.p+1, cols, vals, INSERT_VALUES);
      }
    } 
    
  }
  
  MatAssemblyBegin(Mij, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Mij, MAT_FINAL_ASSEMBLY);
  free(getMeshIndex); 
  PetscFree(vals);
  PetscFree(rows);
  PetscFree(cols);
}                                               



void ks_massReduced(Universe equation, Galaxy* primal, Mat Mij,
		    Is_it *reduced){
  int i, j, k, n, m;  //initialization for iteration                           
  int ij;       //index value combination of element number and basis value    
  PetscScalar *vals;                                                           
  PetscInt *rows, *cols; //indices refering to their respectful index                                                                      
  int *getMeshIndex = createMeshMap(&primal->space);                           
  PetscMalloc1(pow(primal->basis.p+1,2), &vals);                               
  PetscMalloc1(pow(primal->basis.p+1,2), &rows);                               
  PetscMalloc1(pow(primal->basis.p+1,2), &cols);                               
  for (i=0; i<primal->basis.p+1; i++){                                         
    for (j=0; j<primal->basis.p+1; j++)                                        
      vals[j] = 0;                                                             
    for (j=0; j<primal->basis.p+1; j++){                                       
      for (n=0; n<primal->quad.n; n++){                                        
        if (primal->space.dx > 0){                                             
          vals[j] += primal->space.dx/2.0*                                     
            phi(primal->basis.p, i, primal->quad.x_i.array[n])*                
            phi(primal->basis.p, j, primal->quad.x_i.array[n])*                
            primal->quad.w.array[n];                                           
        }                                                                      
        else{                                                                  
          vals[j] = 1;                                                   
        }                                                                      
      }                                                                        
    }                                                                          
    for (j=0; j<reduced->reducedMesh.elem.count; j++){              
      for (m =0; m<equation.numStates; m++){                                  
        ij = (j*equation.numStates+m)*(primal->basis.p+1)+i; //same elment     
        for (k=0; k<primal->basis.p+1; k++)
	  cols[k] = (j*equation.numStates+m)*(primal->basis.p+1)+k;
	MatSetValues(Mij, 1, &ij, primal->basis.p+1, cols, vals, INSERT_VALUES);
      }                                                                       
    }                                                                          
  }                                                                      
  MatAssemblyBegin(Mij, MAT_FINAL_ASSEMBLY);                                   
  MatAssemblyEnd(Mij, MAT_FINAL_ASSEMBLY);                                     
  free(getMeshIndex);                                                          
  PetscFree(vals);                                                             
  PetscFree(rows);                                                             
  PetscFree(cols);                                                             
}     

void ks_massBlock(Universe equation, Galaxy *primal, Mat *MijBlock){
  int i, j, n;                    //initialization for iteration
  int *cols = (int *) malloc ((primal->basis.p+1)*sizeof(int));
  double *vals = (double *) malloc ((primal->basis.p+1)*sizeof(double));
  //---------------------------------------------------------------------------
  // Initialization and Mallocing
  //---------------------------------------------------------------------------
  MatCreateSeqDense(PETSC_COMM_SELF, primal->basis.p+1, primal->basis.p+1,
                    NULL, MijBlock);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<primal->basis.p+1; i++)
    cols[i] = i;
  for (i=0; i<primal->basis.p+1; i++){
    for (j=0; j<primal->basis.p+1; j++){
      vals[j] = 0; 
      for (n=0; n<primal->quad.n; n++){
        if (primal->space.dx > 0){
          vals[j] += primal->space.dx/2.0*
            phi(primal->basis.p, i, primal->quad.x_i.array[n])*
            phi(primal->basis.p, j, primal->quad.x_i.array[n])*
            primal->quad.w.array[n];
          
        }
        else{
          vals[j] = 1;
        }
      }
    }
       MatSetValues(*MijBlock, 1, &i, primal->basis.p+1, cols, vals,INSERT_VALUES);  
   
  }
  
  MatAssemblyBegin(*MijBlock, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*MijBlock, MAT_FINAL_ASSEMBLY);
  free(cols);
  free(vals);
}



void yk_pMass(Universe equation, Galaxy *primal_h, Galaxy* primal_H, Mat *Mij){
  int i, j, n;
  int nElem;
  int row, col;  //position
  PetscScalar value;
  int systemSizeRow = primal_H->space.elem.count*(primal_H->basis.p+1);
  int systemSizeCol = primal_h->space.elem.count*(primal_h->basis.p+1);
  //primal_h is the actual mesh (top)
  //primal_H is what u want to write primal_h as. (bottom)
  MatCreateSeqDense(PETSC_COMM_SELF, primal_H->basis.p+1, primal_h->basis.p+1,
		    NULL, Mij); 
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<primal_H->basis.p+1; i++){
    for (j=0; j<primal_h->basis.p+1; j++){
      for (n=0; n<primal_H->quad.n; n++){
	value = primal_H->space.dx/2.0*
	  phi(primal_H->basis.p, i, primal_H->quad.x_i.array[n])*
	  phi(primal_h->basis.p, j, primal_H->quad.x_i.array[n])*
	  primal_H->quad.w.array[n];
	MatSetValue(*Mij, i, j, value, ADD_VALUES);
      }
    }
  }
  MatAssemblyBegin(*Mij, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Mij, MAT_FINAL_ASSEMBLY);
}

void yk_project_I(Multiverse *multiquation, Galaxy *primal_h,
		  Galaxy *primal_H, Mat *I_proj){
  Mat MijBlock;
  Mat Mij;
  Mat pMij;
  int systemSizeRow = primal_H->space.elem.count*(primal_H->basis.p+1); 
  int systemSizeCol = primal_h->space.elem.count*(primal_h->basis.p+1);  
  /* MatCreateSeqDense(PETSC_COMM_SELF, systemSizeRow, systemSizeCol, */
  /* 		     NULL, I_proj);  */
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  ks_massBlock(multiquation->equation, primal_H, &MijBlock);            
  MatAssemblyBegin(MijBlock, MAT_FINAL_ASSEMBLY);                             
  MatAssemblyEnd(MijBlock, MAT_FINAL_ASSEMBLY);                               
  inverseMatLU(MijBlock, &Mij); 
  
  yk_pMass(multiquation->equation, primal_h, primal_H, &pMij); 
  MatMatMult(Mij, pMij, MAT_INITIAL_MATRIX, PETSC_DEFAULT, I_proj);
  MatDestroy(&MijBlock);
  MatDestroy(&Mij);
  MatDestroy(&pMij);
}


void yk_massH1(Universe equation, Galaxy *primal_c, Galaxy* primal_r,
	       Mat *Mij){
  
  int i, j, n; //initialization for iteration
  PetscScalar value;
  //---------------------------------------------------------------------------
  // Initialization and Mallocing                                              
  //---------------------------------------------------------------------------
  MatCreateSeqDense(PETSC_COMM_SELF, primal_r->basis.p+1, primal_c->basis.p+1,
		    NULL, Mij); 
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<primal_r->basis.p+1; i++){
    for (j=0; j<primal_c->basis.p+1; j++){
      for (n=0; n<primal_r->quad.n; n++){
	value = primal_r->space.dx/2.0*
	  phi(primal_r->basis.p, i, primal_r->quad.x_i.array[n])*
	  phi(primal_c->basis.p, j, primal_r->quad.x_i.array[n])*
	  primal_r->quad.w.array[n]
	  + 2.0/primal_r->space.dx*
	  dphi(primal_r->basis.p, i, primal_r->quad.x_i.array[n])*
	  dphi(primal_c->basis.p, j, primal_r->quad.x_i.array[n])*
	  primal_r->quad.w.array[n];                                         
	MatSetValue(*Mij, i, j, value, ADD_VALUES);                          
      }                                                                      
    }
  }
  MatAssemblyBegin(*Mij, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Mij, MAT_FINAL_ASSEMBLY);
}



void yk_projectH1_I(Multiverse *multiquation, Galaxy *primal_h,            
                  Galaxy *primal_H, Mat *I_proj){        
  Mat MijBlock;                                                          
  Mat Mij;                                                               
  Mat pMij;                                                              
  int systemSizeRow = primal_H->space.elem.count*(primal_H->basis.p+1);  
  int systemSizeCol = primal_h->space.elem.count*(primal_h->basis.p+1);  
  //--------------------------------------------------------------------------
  // Implementation                                                      
  //--------------------------------------------------------------------------
  yk_massH1(multiquation->equation, primal_H, primal_H, &MijBlock);
  inverseMatLU(MijBlock, &Mij);                                          

  yk_massH1(multiquation->equation, primal_h, primal_H, &pMij);          
  MatMatMult(Mij, pMij, MAT_INITIAL_MATRIX, PETSC_DEFAULT, I_proj);      
  MatDestroy(&MijBlock);                                                 
  MatDestroy(&Mij);                                                      
  MatDestroy(&pMij);                                                     
}                         

void yk_Identity(Mat *Imat, int size){
  //Square matrix
  int i;
  //---------------------------------------------------------------------------
  // Initialization
  //---------------------------------------------------------------------------
  MatCreateSeqDense(PETSC_COMM_SELF, size, size, NULL, Imat);  
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<size; i++)
    MatSetValue(*Imat, i, i, 1, INSERT_VALUES);             
  MatAssemblyBegin(*Imat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Imat, MAT_FINAL_ASSEMBLY);
}