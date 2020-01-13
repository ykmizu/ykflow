
#ifndef KS_ADJOINT_UNSTEADY_H_
#define KS_ADJOINT_UNSTEADY_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "struct_def.h"                                                        
#include "ks_dRdU.h"                                                           
#include "ks_mass.h"                                                           
#include "ks_read_solution_file.h"                                             
#include "petsc.h"                                                             
#include "ks_copyUtype.h"                                                      
#include "tools.h"                                                             
#include "ks_function.h"                                                       
#include "ks_jbaraverage.h"

void ks_adjoint_Unsteady(Multiverse multiquation, Cluster *primal);

#endif
