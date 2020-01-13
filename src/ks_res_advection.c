
#include "ks_res_advection.h"

void ks_res_advection(Universe equation, Galaxy *primal, Vec R,
		      Is_it *reduced){
  int i, j, n;      //initialization for iteration
  int i_epsilon;
  int *getMeshIndex = createMeshMap(&primal->space);
  double element;
  double state;
  double residual;
  double c = equation.c[0];
  PetscInt index;
  Mesh _meshOfInterest;
  if (reduced->hrom==1)
    _meshOfInterest = reduced->reducedMesh;
  else
    _meshOfInterest = primal->space;
  //-----------------------------------------------------------------------
  // Implementation
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];   
    for (j=0; j<primal->basis.p+1; j++){
      index = i*(primal->basis.p+1)+j;
      //-------------------------------------------------------------------
      // Integral term here
      //-------------------------------------------------------------------
      for (n=0; n<primal->quad.n; n++){
        state = U_xi(primal, 0, i_epsilon, primal->quad.x_i.array[n]);
	residual = -c*dphi(primal->basis.p, j, primal->quad.x_i.array[n])*state
          *primal->quad.w.array[n];
	VecSetValue(R, index, residual, ADD_VALUES); 
      }
      //-------------------------------------------------------------------
      // Flux terms
      //-------------------------------------------------------------------
      if (c >= 0 ){ //If advection term is positive/ to the right in ID
        if (i_epsilon==0) //Takes into account of periodic boundary condit
          element= 0;//  U_xi(primal, 0, primal->space.elem.count-1, 1);
        else
          element = U_xi(primal, 0, i_epsilon-1, 1);
        residual = c*U_xi(primal, 0, i_epsilon, 1)*phi(primal->basis.p, j, 1)
	  -c*element*phi(primal->basis.p, j, -1);  
      }else if( c < 0 ){ //If advection term is neagtive/ to the left in 1D
        if (i_epsilon==primal->space.elem.count-1)//Takes into account of 
          element = 0;//element  = U_xi(primal, 0, 0,-1);
        else
          element = U_xi(primal, 0, i_epsilon+1, -1);
        residual =-c*U_xi(primal, 0, i_epsilon, -1)*phi(primal->basis.p, j, -1)
	  + c*element*phi(primal->basis.p, j, 1);
      }
      VecSetValue(R, index, residual, ADD_VALUES);
    }
  }
  free(getMeshIndex); 
}
