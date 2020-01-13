
#ifndef YK_CALCULATEERROR_H_
#define YK_CALCULATEERROR_H_

#include <stdio.h>                                                           
#include <stdlib.h>                                                            
#include "math.h"                                                              
#include "struct_def.h"                                                        
#include "ks_Residual.h"                                                       
#include "petsc.h"                                                             
#include "tools.h"                                                             
#include "ks_ApplyTimeScheme.h"                                                
#include "ks_read_solution_file.h"                                             
#include "injection.h"                                                         
#include "ks_jbaraverage.h"                                                    
#include "ks_adjoint_Unsteady.h"                                               
#include "yk_solverCFD.h"

void yk_estimateError(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *primal_h, Cluster *primal_H, Cluster *primal,
		      Is_it *reduced);

#endif
