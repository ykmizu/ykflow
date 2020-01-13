
#include "ks_dres_burgers.h"

void ks_dres_burgers(Universe equation, Galaxy *primal, Mat dRdU,
		     Is_it *reduced){
  int i, j, k, n; //initialization for iteration
  int i_epsilon;
  double state;
  double u1, u2;
  PetscScalar drdu;
  int sign;
  PetscInt row, col, col1, col2;
  int *getMeshIndex = createMeshMap(&primal->space);
  Mesh _meshOfInterest;                                                 
  if (reduced->hrom==1)                                     
    _meshOfInterest = reduced->reducedMesh;                                   
  else
    _meshOfInterest = primal->space; 
  //-----------------------------------------------------------------------
  // I m p l e m e n t a t i o n
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<primal->basis.p+1; j++){
      for (k=0; k<primal->basis.p+1; k++){
        row = i*(primal->basis.p+1)+j;
        col = i_epsilon*(primal->basis.p+1)+k; 
        for (n=0; n<primal->quad.n; n++){
          state = U_xi(primal, 0, i_epsilon, primal->quad.x_i.array[n]);   
          //---------------------------------------------------------------
          // Integral here
          //---------------------------------------------------------------
          drdu = -dphi(primal->basis.p, j, primal->quad.x_i.array[n])*state
            *phi(primal->basis.p, k, primal->quad.x_i.array[n])*
            primal->quad.w.array[n];
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
        }
        for (n=0; n<2; n++){
          if (n==0){ //right
            sign = 1;
            u1 = U_xi(primal, 0, i_epsilon, 1);
            col1 = i_epsilon*(primal->basis.p+1)+k; 
            if (i_epsilon == primal->space.elem.count-1){            //u2 =
              u2 = 0;//U_xi(primal, 0, i_epsilon, -1); //for periodic con
              col2 = i_epsilon*(primal->basis.p+1)+k;  
            }else{
              u2 = U_xi(primal, 0, i_epsilon+1, -1);
              col2 = (i_epsilon+1)*(primal->basis.p+1)+k;  
            }
          }else if (n==1){
            sign = -1;
            if (i_epsilon==0){                 //Calculate the left solu
              u1 = 0;//U_xi(primal, 0, 0 , 1); //u1 = 0
              //col1 = 0;
              col1 = i_epsilon*(primal->basis.p+1)+k;   
            }else{
              u1 = U_xi(primal, 0, i_epsilon-1, 1);
              col1 = (i_epsilon-1)*(primal->basis.p+1)+k;  
            }
            u2 = U_xi(primal, 0, i_epsilon, -1);
            col2 = i_epsilon*(primal->basis.p+1)+k;
          }
          if (u1 <= u2){
            if (u1*u2 < 0)
              drdu = 0;
            else
              if (0.5*u1*u1 <= 0.5*u2*u2){
                drdu = sign*u1*phi(primal->basis.p, k, 1)*
                  phi(primal->basis.p, j, sign);
                col = col1;
              }else if (0.5*u1*u1 > 0.5*u2*u2){
                drdu = sign*u2*phi(primal->basis.p, k, -1)*
                  phi(primal->basis.p, j, sign);
                col = col2;
              }
          }else if(u1 >=u2){
            if (0.5*u1*u1 >= 0.5*u2*u2){
              drdu = sign*u1*phi(primal->basis.p, k, 1)*
                phi(primal->basis.p, j, sign);
              col = col1;
            }else if (0.5*u1*u1 < 0.5*u2*u2){
              drdu = sign*u2*phi(primal->basis.p, k, -1)*
                phi(primal->basis.p, j, sign);
              col = col2;
            }
          }
          MatSetValue(dRdU, row, col, drdu, ADD_VALUES);  
        }
      }
    }
  }
  free(getMeshIndex); 
}

