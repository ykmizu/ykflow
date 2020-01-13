
#ifndef KS_GREEDYALGORITHM_H_
#define KS_GREEDYALGORITHM_H_

#include <stdio.h>                                                             
#include <stdlib.h>                                                            
#include "math.h"                                                              
#include "struct_def.h"                                                        
#include "tools.h"                                                             
#include "petsc.h"                                                             
#include "ks_read_solution_file.h"                                             
#include "ks_Residual.h"                                                       
#include "ks_dRdU.h"                                                           
#include "ks_podv2.h"                                                          
#include "ks_gaussNewtonSolve.h"                                               
#include <time.h>                                                              
#include "system_initialization.h"                                             
#include "qRFactorization.h"                                                   
#include "ks_offlineMatrix.h"                                                  
#include "ks_copyUtype.h" 

void ks_greedyAlgorithm(int seedCount, PetscInt *nodeSet, Is_it *reduced);

void partialSortInt(PetscInt i[], int first, int last);
#endif