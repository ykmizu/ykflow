
#include "ks_dres_fourth_IPDG.h"

void ks_dres_fourth_IPDG(Universe equation, Galaxy *primal, Mat dRdU,
			 Is_it *reduced){
  int i, j, k, n;   //initialization for iteration
  int i_epsilon;
  double sigma= 1*pow(primal->basis.p, 6)/pow(primal->space.dx,3); //first i
  double tau = 1*pow(primal->basis.p, 2)/pow(primal->space.dx,1);   //second
  int vplus; //existence of basis function on left side
  int vneg;  //existence of basis function on right side
  double elemplus; //element to extract information from left hand side
  double delemplus;
  double d2elemplus;
  double d3elemplus;
  double elemneg;  //element to extract information from right hand side
  double delemneg;
  double d2elemneg;
  double d3elemneg;
  double c = equation.c[3];
  PetscScalar drdu;
  PetscInt row, col, col2;
  int *getMeshIndex = createMeshMap(&primal->space);
  Mesh _meshOfInterest;                                                   
  if (reduced->hrom==1)                                      
    _meshOfInterest = reduced->reducedMesh;
  else                                                                    
    _meshOfInterest = primal->space; 
  //-----------------------------------------------------------------------
  // Implementation
  //-----------------------------------------------------------------------
  for (i=0; i<_meshOfInterest.elem.count; i++){  //Iterate through al
    i_epsilon = getMeshIndex[_meshOfInterest.elem.array[i]];
    for (j=0; j<primal->basis.p+1; j++){  //Iterate through all the basis func
      for (k=0; k<primal->basis.p+1; k++){
        row = i*(primal->basis.p+1)+j;
        col = i_epsilon*(primal->basis.p+1)+k;  
        for (n=0; n<primal->quad.n; n++){ //Interate through quadrature poi
          //---------------------------------------------------------------
          // Integral term here
          //---------------------------------------------------------------
          drdu = equation.c[3]*
            primal->space.dx*(2.0/primal->space.dx)*(2.0/primal->space.dx)*
            d2phi(primal->basis.p, k, primal->quad.x_i.array[n])*
            (2.0/primal->space.dx)*(2.0/primal->space.dx)*
            d2phi(primal->basis.p, j, primal->quad.x_i.array[n])*
            primal->quad.w.array[n];
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
        }
      }
    }
    for (k=0; k<primal->basis.p+1; k++){
      for (n=0; n< 2; n++){
        if (n==0){// k+1/2
            vplus = 0;
            vneg = 1;
            elemneg = phi(primal->basis.p, k, -1);
            delemneg = dphi(primal->basis.p, k, -1)*(2.0/primal->space.dx);
            d2elemneg = d2phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
              (2.0/primal->space.dx);
            d3elemneg = d3phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
              (2.0/primal->space.dx)*(2.0/primal->space.dx);
            col = i_epsilon*(primal->basis.p+1)+k;
            if (i_epsilon==0){
              d3elemplus = d3phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx)*(2.0/primal->space.dx);
              d2elemplus = d2phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx);
              delemplus = 0;
              elemplus = 0;
              col2 = i_epsilon*(primal->basis.p+1)+k;
            }else{
              d3elemplus = d3phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx)*(2.0/primal->space.dx);
              d2elemplus = d2phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx);
              delemplus = dphi(primal->basis.p, k, 1)*(2.0/primal->space.dx);
              elemplus = phi(primal->basis.p, k, 1);
              col2 = (i_epsilon-1)*(primal->basis.p+1)+k;
            }
          }else if(n==1){
            vplus = 1;
            vneg = 0;
            elemplus = phi(primal->basis.p, k, 1);
            delemplus = dphi(primal->basis.p, k, 1)*(2.0/primal->space.dx);
            d2elemplus = d2phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
              (2.0/primal->space.dx);
            d3elemplus = d3phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
              (2.0/primal->space.dx)*(2.0/primal->space.dx);
            col2 = i_epsilon*(primal->basis.p+1)+k;
            if (i_epsilon==primal->space.elem.count-1){
              d3elemneg = d3phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx)*(2.0/primal->space.dx); //t Boundary Con
              d2elemneg = d2phi(primal->basis.p, k, 1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx);
              delemneg = 0;
              elemneg = 0;
              col = i_epsilon*(primal->basis.p+1)+k;
            }else{
              d3elemneg = d3phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx)*(2.0/primal->space.dx);
              d2elemneg = d2phi(primal->basis.p, k, -1)*(2.0/primal->space.dx)*
                (2.0/primal->space.dx);
              delemneg = dphi(primal->basis.p, k, -1)*(2.0/primal->space.dx);
              elemneg = phi(primal->basis.p, k, -1);
              col = (i_epsilon+1)*(primal->basis.p+1)+k;
            }
          }
        for (j=0; j<primal->basis.p+1; j++){
          row = i*(primal->basis.p+1)+j;
          //Number 1
          drdu = c*0.5*d3elemplus*
            (phi(primal->basis.p, j, 1)*vplus-phi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu = c*0.5*d3elemneg*
            (phi(primal->basis.p, j, 1)*vplus-phi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
          
          //Number 2
          drdu =elemplus*c*0.5*
            (2.0/primal->space.dx)*(2.0/primal->space.dx)*(2.0/primal->space.dx)
            *(d3phi(primal->basis.p,j, 1) + d3phi(primal->basis.p, j, -1));
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu= -elemneg*c*0.5*
            (2.0/primal->space.dx)*(2.0/primal->space.dx)*(2.0/primal->space.dx)
            *(d3phi(primal->basis.p,j, 1) + d3phi(primal->basis.p, j, -1));
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
          
          /* //Number 3 */
          drdu = -c*0.5*(2.0/primal->space.dx)*d2elemplus*
            (dphi(primal->basis.p, j, 1)*vplus-
             dphi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu = -c*0.5*(2.0/primal->space.dx)*d2elemneg*
            (dphi(primal->basis.p, j, 1)*vplus-
             dphi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
          
          //Number 4
          drdu = -c*0.5*(2.0/primal->space.dx)*(2.0/primal->space.dx)*
            (d2phi(primal->basis.p, j, 1)*vplus +
             d2phi(primal->basis.p, j, -1)*vneg)*delemplus;
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu = +c*0.5*(2.0/primal->space.dx)*(2.0/primal->space.dx)*
            (d2phi(primal->basis.p, j, 1)*vplus +
             d2phi(primal->basis.p, j, -1)*vneg)*delemneg;
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
          
          //Number 5
          drdu = c*sigma*elemplus*
            (phi(primal->basis.p, j, 1)*vplus-phi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu = -c*sigma*elemneg*
            (phi(primal->basis.p, j, 1)*vplus-phi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
          
          //Number 6
          drdu = c*tau*(2.0/primal->space.dx)*delemplus*
            (dphi(primal->basis.p, j, 1)*vplus-
             dphi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col2, &drdu, ADD_VALUES);
          drdu = -c*tau*(2.0/primal->space.dx)*delemneg*
            (dphi(primal->basis.p, j, 1)*vplus
             - dphi(primal->basis.p, j, -1)*vneg);
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);

        }
      }
    }
  }
  free(getMeshIndex);
}
