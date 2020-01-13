
//-----------------------------------------------------------------------------
//m = the number of elements you want to find residual values for
//idxm = the global indices of the elements you want values for
//n = the number of basis values you want to get valuesin each element.
//idxn = the globl indices of the basis you want values for in each element
//
// Yukiko Shimizu
// February 16, 2017
//------------------------------------------------------------------------------

#include "ks_res_fourth_IPDG.h"

void ks_res_fourth_IPDG(Universe equation, Galaxy *primal, Vec R,
                        Is_it *reduced){
  int i, j, n;                         //initialization for iteration
  int i_epsilon;
  int *getMeshIndex = createMeshMap(&primal->space);
  double d2state;                      //2nd partial derivative of the stat
  double sigma=1*pow(primal->basis.p, 6)/pow(primal->space.dx,3); //ft i
  double tau = 1*pow(primal->basis.p, 2)/pow(primal->space.dx,1);   //se
  int vplus;                           //basis function on left side
  int vneg;                            //basis function on right side
  double c = equation.c[3];                 //Fourth order magnitude 
  double elemplus;                     //element information from LHS
  double delemplus;                    //first derivative LHS
  double d2elemplus;                   //second derivative LHS
  double d3elemplus;                   //third derivative LHS
  double elemneg;                      //element information from RHS
  double delemneg;                     //first derivative RHS 
  double d2elemneg;                    //second derivative RHS
  double d3elemneg;                    //third derivative RHS
  double residual;
  Mesh _meshOfInterest;                                                    
  PetscInt index;
  if (reduced->hrom==1)                                      
    _meshOfInterest = reduced->reducedMesh;
  else
    _meshOfInterest = primal->space; 
  //-----------------------------------------------------------------------
  // Implementation
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){   //Iterate through a
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<primal->basis.p+1; j++){    //Iterate through all basis fu
      index = i*(primal->basis.p+1)+j;  
      //-------------------------------------------------------------------
      // Integral term Conditions
      //-------------------------------------------------------------------
      for (n=0; n<primal->quad.n; n++){  //Interate through quadrature poi
        d2state = d2U_xi(primal, 0, i_epsilon, primal->quad.x_i.array[n]);
        residual = c*primal->space.dx*(2.0/primal->space.dx)*
	  (2.0/primal->space.dx)*d2state*
	  d2phi(primal->basis.p, j, primal->quad.x_i.array[n])*
          primal->quad.w.array[n];
	VecSetValue(R, index, residual, ADD_VALUES); 
      }
    }
    //---------------------------------------------------------------------
    // Flux terms
    //---------------------------------------------------------------------
    for (n=0; n<2; n++){               //Iterate through both left and r
      if (n==0){                       //Left first
        vplus = 0;
        vneg = 1;
        elemneg = U_xi(primal, 0, i_epsilon, -1);
        delemneg = dU_xi(primal, 0, i_epsilon, -1);
        d2elemneg = d2U_xi(primal, 0, i_epsilon, -1);
        d3elemneg = d3U_xi(primal, 0, i_epsilon, -1);
        if (i_epsilon == 0){            //If operating on the left boundary
          d3elemplus = d3U_xi(primal, 0, i_epsilon, -1); //Direchlet Bound
          d2elemplus = d2U_xi(primal, 0, i_epsilon, -1);
          delemplus = 0;
          elemplus = 0;
        }else{
          d3elemplus = d3U_xi(primal, 0, i_epsilon-1, 1);
          d2elemplus = d2U_xi(primal, 0, i_epsilon-1, 1);
          delemplus = dU_xi(primal, 0, i_epsilon-1, 1);
          elemplus = U_xi(primal, 0, i_epsilon-1, 1);
        }
      }else if(n == 1){                 //Right next
        vplus = 1;
        vneg = 0;
        elemplus = U_xi(primal, 0, i_epsilon, 1);
        delemplus = dU_xi(primal, 0, i_epsilon, 1);
        d2elemplus = d2U_xi(primal, 0, i_epsilon, 1);
        d3elemplus = d3U_xi(primal, 0, i_epsilon, 1);
        if (i_epsilon == primal->space.elem.count-1){ //If operating on the
          d3elemneg = d3U_xi(primal, 0, i_epsilon, 1); //Direchlet Bound
          d2elemneg = d2U_xi(primal, 0, i_epsilon, 1);
          delemneg = 0;
          elemneg = 0;
        }else{
          d3elemneg = d3U_xi(primal, 0, i_epsilon+1, -1);
          d2elemneg = d2U_xi(primal, 0, i_epsilon+1, -1);
          delemneg = dU_xi(primal, 0, i_epsilon+1, -1);
          elemneg = U_xi(primal, 0, i_epsilon+1, -1);
        }
      }
      for (j=0; j<primal->basis.p+1; j++){ //Calculate the boundary value
	index = i*(primal->basis.p+1)+j;   
	residual =
	  //Number 1
	  c*0.5*(d3elemplus + d3elemneg)*
          (phi(primal->basis.p, j, 1)*vplus-phi(primal->basis.p, j, -1)*vneg)
	  //Number 2
          +c*0.5*(2.0/primal->space.dx)*(2.0/primal->space.dx)*
          (2.0/primal->space.dx)*
          (d3phi(primal->basis.p, j, 1)+ d3phi(primal->basis.p, j, -1))*
          (elemplus - elemneg)
	  //Number 3
          -c*0.5*(2.0/primal->space.dx)*(d2elemplus + d2elemneg)*
          (dphi(primal->basis.p, j, 1)*vplus -
           dphi(primal->basis.p, j, -1)*vneg)
	  //Number 4
	  -c*0.5*(2.0/primal->space.dx)*(2.0/primal->space.dx)*
          (d2phi(primal->basis.p, j, 1)*vplus + d2phi(primal->basis.p, j, -1)*
           vneg)*(delemplus - delemneg)
	  //Number 5
	  +c*sigma*(elemplus - elemneg)*
	  (phi(primal->basis.p, j, 1)*vplus - phi(primal->basis.p, j, -1)
	   *vneg)
	  //Number 6
	  +c*tau*(2.0/primal->space.dx)*(delemplus - delemneg)*
          (dphi(primal->basis.p, j, 1)*vplus -
           dphi(primal->basis.p, j, -1)*vneg);
	 VecSetValue(R, index, residual, ADD_VALUES);  
      }
    }
  }
  free(getMeshIndex); 
}
