
#include "ks_res_burgers.h" 

void ks_res_burgers(Universe equation, Galaxy *primal, Vec R, Is_it *reduced){
  int i, j, n;                      //initialization for iteration
  int i_epsilon;
  double state;                        //Approximation
  double u1, u2;                       //Left and Right solution approximat
  double residual;
  int sign;
  int *getMeshIndex = createMeshMap(&primal->space);
  PetscInt index;
  Mesh _meshOfInterest;
  if (reduced->hrom==1)        
    _meshOfInterest = reduced->reducedMesh;                                    
  else
    _meshOfInterest = primal->space;
  //-----------------------------------------------------------------------
  // Implementation
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){ //Iterate through de
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<primal->basis.p+1; j++){    //Iterate through the basis poi
      index = i*(primal->basis.p+1)+j;
      //-------------------------------------------------------------------
      // Integral term here                                            
      //-------------------------------------------------------------------
      for (n=0; n<primal->quad.n; n++){   //Iterate through the quad points
        state = U_xi(primal, 0, i_epsilon, primal->quad.x_i.array[n]);
	residual = -0.5*dphi(primal->basis.p, j, primal->quad.x_i.array[n])*
	  state*state*primal->quad.w.array[n];
	VecSetValue(R, index, residual, ADD_VALUES);
      }
      //-------------------------------------------------------------------
      // Flux terms                                                        
      //-------------------------------------------------------------------
      for (n=0; n<2; n++){
        if (n==0){ //right
          sign = 1;
          u1 = U_xi(primal, 0, i_epsilon, 1);
          if (i_epsilon == primal->space.elem.count-1)            //u2 = 0;
            u2 = 0;//U_xi(U, 0, i, -1); //for periodic con
          else
            u2 = U_xi(primal, 0, i_epsilon+1, -1);  
        }else if (n==1){
          sign = -1;
          if (i_epsilon==0)                 //Calculate the left solution at el
            u1 = 0; //U_xi(U, 0, 0 , 1); //u1 = 0
          else
            u1 = U_xi(primal, 0, i_epsilon-1, 1); 
          u2 = U_xi(primal, 0, i_epsilon, -1);
        }
        if (u1 <= u2){
          if (u1*u2 < 0)
            residual = 0; 
          else
            if (0.5*u1*u1 <= 0.5*u2*u2)
              residual =sign*0.5*u1*u1*phi(primal->basis.p, j, sign); 
            else if (0.5*u1*u1 > 0.5*u2*u2)
              residual = sign*0.5*u2*u2*phi(primal->basis.p, j, sign);
        }else if(u1 >=u2){
          if (0.5*u1*u1 >= 0.5*u2*u2)
            residual=sign*0.5*u1*u1*phi(primal->basis.p, j, sign); 
          else if (0.5*u1*u1 < 0.5*u2*u2)
            residual= sign*0.5*u2*u2*phi(primal->basis.p, j, sign); 
        }
	VecSetValue(R, index, residual, ADD_VALUES);
      }
    }
  }
  free(getMeshIndex); 
}
