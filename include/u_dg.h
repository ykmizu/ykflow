/*
 ============================================================================
 Name        : LorenzError.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
 */

#ifndef U_DG_H_
#define U_DG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct_def.h"
#include "basisfunc.h"

double U_xi(Galaxy *U, int comU, int element, double xi);

double dU_xi(Galaxy *U, int comU, int element, double xi);

double d2U_xi(Galaxy *U, int comU, int element, double xi);

double d3U_xi(Galaxy *U, int comU, int element, double xi);

double d4U_xi(Galaxy *U, int comU, int element, double xi);     
#endif
