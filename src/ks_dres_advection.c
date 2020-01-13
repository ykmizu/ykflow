

#include "ks_dres_advection.h"

void ks_dres_advection(Universe equation, Galaxy *primal, Mat dRdU,
		       Is_it *reduced){
  int i, j, k ,n;                       //initialization for iteration
  int i_epsilon;
  double u_elem;    
  PetscScalar drdu;
  PetscInt row, col;
  int *getMeshIndex = createMeshMap(&primal->space);
  Mesh _meshOfInterest;
  if (reduced->hrom==1)                                       
    _meshOfInterest = reduced->reducedMesh;                                  
  else
    _meshOfInterest = primal->space; 
  //-----------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){ //iterate through all
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<primal->basis.p+1; j++){ //row
      for (k=0; k<primal->basis.p+1; k++){ //column within
        row = i*(primal->basis.p+1)+j; //same element
        col = i_epsilon*(primal->basis.p+1)+k; //same element
        for (n=0; n<primal->quad.n; n++){
          drdu = -equation.c[0]*
            dphi(primal->basis.p, j, primal->quad.x_i.array[n])*
            phi(primal->basis.p, k, primal->quad.x_i.array[n])*
            primal->quad.w.array[n];
	    MatSetValue(dRdU, row, col, drdu, ADD_VALUES);
        }
        if (equation.c[0]>=0){
          if (i_epsilon == 0)
            u_elem = 0;
          //u_elem = phi(primal->basis.p, k, 1);
          else
            u_elem = phi(primal->basis.p, k, 1);
          drdu = equation.c[0]*
            phi(primal->basis.p, k, 1)*phi(primal->basis.p, j, 1);
	  MatSetValue(dRdU, row, col, drdu, ADD_VALUES);
          drdu = -equation.c[0]*u_elem*phi(primal->basis.p, j, -1);   
          if (i_epsilon == 0)
            col = i_epsilon*(primal->basis.p+1)+k; //col = (primal->space.elem.count-1)*(primal->basis.p+1)+k;
          else
            col = (i_epsilon-1)*(primal->basis.p+1)+k;
	  MatSetValue(dRdU, row, col, drdu, ADD_VALUES); 
        }else if (equation.c[0]<0){
          if (i_epsilon == primal->space.elem.count-1)
            u_elem = 0;
            //u_elem = phi(primal->basis.p, k, -1); 
          else
            u_elem = phi(primal->basis.p, k, -1);
          drdu = -equation.c[0]*
            phi(primal->basis.p, k, -1)*phi(primal->basis.p, j, -1);
	  MatSetValue(dRdU, row, col, drdu, ADD_VALUES);
          drdu = equation.c[0]*u_elem*phi(primal->basis.p, j, 1); 
          if (i_epsilon == primal->space.elem.count-1)
            col = i_epsilon*(primal->basis.p+1)+k;
          //col = 0;
          else
            col = (i_epsilon+1)*(primal->basis.p+1)+k;
	  MatSetValue(dRdU, row, col, drdu, ADD_VALUES);
        }
      } 
    }
  }
  free(getMeshIndex); 
}

 
