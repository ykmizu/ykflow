/*
 ==========================================================================
 Name        : LorenzError.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ==========================================================================
 */

#include "u_dg.h"

double U_xi(Galaxy *U, int comU, int element, double xi){
  int i;                             //initialization for iteration
  double sol_xi=0;                   //solution approximation zeroed out
  for (i=0; i<U->basis.p+1; i++){
    sol_xi += U->solution[comU].array[i+element*(U->basis.p+1)]*
      phi(U->basis.p, i, xi);
  }
  return sol_xi;
}

double dU_xi(Galaxy *U, int comU, int element, double xi){
  int i;
  double sol_xi=0;
  for (i=0; i<U->basis.p+1; i++){
    sol_xi+=(2.0/U->space.dx)*U->solution[comU].array[i+element*(U->basis.p+1)]*
      dphi(U->basis.p, i, xi);
  }
  return sol_xi;
}  

double d2U_xi(Galaxy *U, int comU, int element, double xi){
  int i;   //initialization for interpolation purposes
  double sol_xi = 0;
  for (i=0; i<U->basis.p+1; i++){
    sol_xi += (2.0/U->space.dx)*(2.0/U->space.dx)*
      U->solution[comU].array[i+element*(U->basis.p+1)]*
      d2phi(U->basis.p, i, xi);
  }
  return sol_xi;
}

double d3U_xi(Galaxy *U, int comU, int element, double xi){
  int i;   //initialization for iteration purposes
  double sol_xi = 0;
  for (i=0; i<U->basis.p+1; i++){
    sol_xi+= (2.0/U->space.dx)*(2.0/U->space.dx)*(2.0/U->space.dx)*
      U->solution[comU].array[i+element*(U->basis.p+1)]*
      d3phi(U->basis.p, i, xi);
  }
  return sol_xi;
}

double d4U_xi(Galaxy *U, int comU, int element, double xi){
  int i;
  double sol_xi = 0;
  for (i=0; i<U->basis.p+1; i++){
    sol_xi+= (2.0/U->space.dx)*(2.0/U->space.dx)*(2.0/U->space.dx)*
      (2.0/U->space.dx)*U->solution[comU].array[i+element*(U->basis.p+1)]*
      d4phi(U->basis.p, i, xi);
  }
  return sol_xi;
}



                           
     
       
