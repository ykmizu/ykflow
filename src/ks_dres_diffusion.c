
#include "ks_dres_diffusion.h"

void ks_dres_diffusion(Universe equation, Galaxy *primal, Mat dRdU,
		       Is_it *reduced){
  int i, j, k, n; //initialization for iteration
  int i_epsilon;
  double delemnegPrt1, delemnegPrt2, elemneg;  
  double eta = 10; //For now 10 for p =3 2 for p = 2;
  PetscInt row, col;
  PetscScalar drdu;
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
      for (k=0; k<primal->basis.p+1;k++){
        row = i*(primal->basis.p+1)+j;  
        col = i_epsilon*(primal->basis.p+1)+k; 
        for (n=0; n<primal->quad.n; n++){
          //---------------------------------------------------------------
          // Integral term her
          //---------------------------------------------------------------
          drdu = equation.c[2]*
            dphi(primal->basis.p, j, primal->quad.x_i.array[n])*
            (2.0/primal->space.dx)*
            dphi(primal->basis.p, k,  primal->quad.x_i.array[n])*
            primal->quad.w.array[n];
          MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
       }
        //-----------------------------------------------------------------
        // F i +1/2
        //-----------------------------------------------------------------
        if (i_epsilon==primal->space.elem.count-1){
          elemneg = 0;   //Direchlet Boundary Conditions
          delemnegPrt1 = 0;
          delemnegPrt2 = 0;
        }else{
          elemneg = phi(primal->basis.p, k, -1); //i+1
          delemnegPrt1 = +(2.0/primal->space.dx)*0.5*dphi(primal->basis.p, k, 1)
            -eta/primal->space.dx*phi(primal->basis.p, k, 1);
          delemnegPrt2 = (2.0/primal->space.dx)*0.5*dphi(primal->basis.p, k, -1)
            +eta/primal->space.dx*elemneg;
        }
        col = i_epsilon*(primal->basis.p+1)+k;
        drdu = -equation.c[2]*(1.0/primal->space.dx)*
          phi(primal->basis.p, k, 1)*dphi(primal->basis.p, j, 1)
          -equation.c[2]*phi(primal->basis.p, j, 1)*delemnegPrt1;
        MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
        if (i_epsilon==primal->space.elem.count-1)
          col = i_epsilon*(primal->basis.p+1)+k;
        else
          col = (i_epsilon+1)*(primal->basis.p+1)+k;
        drdu = equation.c[2]*(1.0/primal->space.dx)*elemneg*
          dphi(primal->basis.p, j, 1)-equation.c[2]*phi(primal->basis.p, j, 1)*
          delemnegPrt2;
        MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
        //---------------------------------------------------------------------
        // F i -1/2
        //---------------------------------------------------------------------
        if (i_epsilon==0){
          elemneg = 0;    //Direchlet Boundary Conditions
          delemnegPrt1 = 0;
          delemnegPrt2 = 0;
        }else{
          elemneg = phi(primal->basis.p, k, 1); //i-1
          delemnegPrt1 = (2.0/primal->space.dx)*0.5*dphi(primal->basis.p, k,-1)
            +eta/primal->space.dx*phi(primal->basis.p,k,-1);
          delemnegPrt2 = (2.0/primal->space.dx)*0.5*dphi(primal->basis.p, k, 1)
            -eta/primal->space.dx*elemneg; //i-1
          }
        col = i_epsilon*(primal->basis.p+1)+k;
        drdu = +equation.c[2]*(1.0/primal->space.dx)*
          phi(primal->basis.p,k,-1)*dphi(primal->basis.p,j,-1)
          +equation.c[2]*phi(primal->basis.p,j, -1)*delemnegPrt1;
        
        MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
        if (i_epsilon==0)
          col = i_epsilon*(primal->basis.p+1)+k;
        else
          col= (i_epsilon-1)*(primal->basis.p+1)+k;
        drdu = -equation.c[2]*(1.0/primal->space.dx)*elemneg*
          dphi(primal->basis.p, j, -1)
          +equation.c[2]*phi(primal->basis.p, j, -1)*delemnegPrt2;
        MatSetValues(dRdU, 1, &row, 1, &col, &drdu, ADD_VALUES);
      } 
    }
  }
  free(getMeshIndex); 
}

