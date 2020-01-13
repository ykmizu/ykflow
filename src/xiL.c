/*
 ============================================================================
 Name        : basisNquad.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
 */
#include "xiL.h"

double xiL(int r, int k){
  //================Initialization of Parameters and Variables=================
  double dxi;  
  if (r>0){
    dxi = 2.0/r;   //Makes constant distant between each phi node
  }else{
    dxi = 2.0;
  }
  double xi_;
  //===========================Implementation==================================
  xi_=-1.0+k*dxi;  //Calculates the all the xi nodes
  
  return xi_;
}
