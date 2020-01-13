
#ifndef YK_LEASTSQUARESSHADOWING_H_
#define YK_LEASTSQUARESSHADOWING_H_

#include <stdio.h>                                                             
#include <stdlib.h>                                                            
#include "math.h"                                                              
#include "struct_def.h"                                                        
#include "petsc.h"                                                             
#include "injection.h"                                                         
#include "ks_GMRES.h"                                                          
#include "ks_adjoint_MATVEC.h"                                                 
#include "petsc.h"                                                             
#include <stdarg.h>                                                            
#include "yk_createHROM.h"

double yk_leastSquaresShadowing(yk_PrimalSolver *ykflow,
				Multiverse *multiquation, Cluster *primal, 
                                char *argv[], Is_it *reduced);

#endif
