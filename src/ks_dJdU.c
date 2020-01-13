

#include "ks_dJdU.h"
//-----------------------------------------------------------------------------
// Implementation of the derivative of the output of interest
//
// Yukiko Shimizu
// November 8, 2016
// Error Estimation for Chaotic Systems
//-----------------------------------------------------------------------------

void ks_dJdU(Universe eqn, Galaxy *U, Vec dJdU, int time_i){
  int i, j, k, n;
  double temp_dJdU;
  PetscInt vec_index;
  VecZeroEntries(dJdU); //Make sure to initialize dJdU to all zeros first
  if(U->space.node.count > 1){ //PDE
    for (i=0; i< U->space.elem.count; i++){
      for (j=0; j<eqn.numStates; j++){
        for (k=0; k<U->basis.p+1; k++){
          temp_dJdU = 0;
          vec_index = (i*eqn.numStates+j)*(U->basis.p+1)+k;
          for (n=0; n<U->quad.n; n++){
            if (time_i==0 || time_i==U->time.count){   
              temp_dJdU += -U->time.dt/(2*U->time.globalT_f*U->space.x_f)*
                phi(U->basis.p, k, U->quad.x_i.array[n])*U->quad.w.array[n];
            }else{
              temp_dJdU += -U->time.dt/(U->time.globalT_f*U->space.x_f)*
                phi(U->basis.p, k, U->quad.x_i.array[n])*U->quad.w.array[n];
            }
          }
          VecSetValues(dJdU, 1, &vec_index, &temp_dJdU, ADD_VALUES);
        }
      }
    }
  }else{
    temp_dJdU =  1;
    vec_index = 2;
    VecSetValues(dJdU, 1, &vec_index, &temp_dJdU, INSERT_VALUES);
  }
  VecAssemblyBegin(dJdU); //Confused if this should be done once or multiple
  VecAssemblyEnd(dJdU);   //times   
}
