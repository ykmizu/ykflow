/*
 ============================================================================
 Name        : residual.h
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
 */

#ifndef RESIDUAL_H_
#define RESIDUAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct_def.h"
#include "petsc.h"

void residual(Utype *U, int t, double *mR);
#endif
