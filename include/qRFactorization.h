

#ifndef QRFACTORIZATION_H_
#define QRFACTORIZATION_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "struct_def.h"
#include "tools.h"
#include "ks_copyUtype.h"
#include "petsc.h"
#include "ks_read_solution_file.h"
#include "ks_Residual.h"
#include "ks_dRdU.h"
#include "ks_podv2.h"  

void linearLeastSquares(Mat A, Vec b, Vec x);
void qRFactorization(Mat A, Mat Q, Mat R);

#endif
