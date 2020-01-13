
//Need a function to pull off the coefficients for the timesolver to use
//for later
#include "getBDFCoef.h"
//initialization for the time solver before you actually go and do it
double* getBDFCoef(int timeSolveNum){
  double * c = (double *) malloc ((timeSolveNum+1)*sizeof(double));
  //----------------------------------------------------------------------------
  // Extract the coefficients for the backward time solver
  //----------------------------------------------------------------------------
  if (timeSolveNum == 1){    //assigns coefficients based on method of choice
    c[0] = 1; c[1] = -1;
  }else if(timeSolveNum == 2){   //BDF 2
    c[0] = 3.0/2.0; c[1] = -2.0; c[2] = 1.0/2.0;
  }else if (timeSolveNum == 3){  //BDF 3
    c[0] = 11.0/6.0; c[1] = -3.0; c[2] = 3.0/2.0; c[3] = -1.0/3.0;
  }
  return c;
}


