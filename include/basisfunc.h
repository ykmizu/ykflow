/*
 ============================================================================
 Name        : basisNquad.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
 */

#ifndef BASISFUNC_H_
#define BASISFUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "xiL.h"

double phi(int r, int k, double xi);
double dphi(int r, int k, double xi);
double d2phi(int r, int k, double xi);
double d3phi(int r, int k, double xi);
double d4phi(int r, int k, double xi);

#endif
