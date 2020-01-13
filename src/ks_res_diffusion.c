
//-------------------------------------------------------------------------
//m = the number of elements you want to find residual values for
//idxm = the global indices of the elements you want values for
//n = the number of basis values you want to get valuesin each element.
//idxn = the globl indices of the basis you want values for in each element
//
// Yukiko Shimizu
// February 16, 2017
//-------------------------------------------------------------------------

#include "ks_res_diffusion.h"

void ks_res_diffusion(Universe equation, Galaxy *primal, Vec R,
                      Is_it *reduced){
  int i, j, n;                      //initialization for iteration
  int i_epsilon;
  int *getMeshIndex = createMeshMap(&primal->space);
  double dstate;                       //Derivative state
  int elempos;                         //index element number
  double delemneg, elemneg;            //State values and derivative state
  double residual;
  PetscInt index;
  double alfa = equation.c[2];              //Magnitude of diffusion
  double eta = 10;                      //For now 10 for p =3 2 for p = 2;
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
        dstate = dU_xi(primal, 0, i_epsilon, primal->quad.x_i.array[n]);      
        residual = alfa*dphi(primal->basis.p, j, primal->quad.x_i.array[n])
	  *dstate*primal->quad.w.array[n];        
	VecSetValue(R, index, residual, ADD_VALUES);
      }
      elempos = i_epsilon;    
      //-------------------------------------------------------------------
      // F i +1/2
      //-------------------------------------------------------------------
      if (i_epsilon==primal->space.elem.count-1){
        elemneg = 0;   //Direchlet Boundary Conditions
        delemneg = 0;
      }else{
        elemneg = U_xi(primal, 0, i_epsilon+1, -1);
        delemneg = 0.5*(dU_xi(primal, 0, elempos, 1) +
                        dU_xi(primal, 0, i_epsilon+1, -1))-
          eta/primal->space.dx*(U_xi(primal, 0, elempos, 1)-elemneg);
      }
      residual = -alfa*(1.0/primal->space.dx)*U_xi(primal, 0, elempos, 1)
        *dphi(primal->basis.p, j, 1)
	+alfa*(1.0/primal->space.dx)*elemneg*dphi(primal->basis.p, j, 1)
	-alfa*phi(primal->basis.p, j, 1)*delemneg;
      VecSetValue(R, index, residual, ADD_VALUES);
      //-------------------------------------------------------------------
      // F i -1/2
      //-------------------------------------------------------------------
      if (i_epsilon==0){
        elemneg = 0;    //Direchlet Boundary Conditions
        delemneg = 0;
      }else{
        elemneg = U_xi(primal, 0, i_epsilon-1, 1);
        delemneg = 0.5*(dU_xi(primal, 0, elempos, -1) +
			dU_xi(primal, 0, i_epsilon-1, 1))-
          eta/primal->space.dx*(-U_xi(primal, 0, elempos, -1)+elemneg);
      }
      residual = alfa*(1.0/primal->space.dx)*U_xi(primal, 0, elempos, -1)
	*dphi(primal->basis.p, j, -1)
	-alfa*(1.0/primal->space.dx)*elemneg*dphi(primal->basis.p, j, -1)
	+alfa*phi(primal->basis.p, j, -1)*delemneg;
      VecSetValue(R, index, residual, ADD_VALUES);
    }
  }
  free(getMeshIndex); 
}

