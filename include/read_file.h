
/*                                                                
 ============================================================================  
 Name        : LorenzError.c                                                       
 Author      : Yukiko Shimizu                                       
 Version     : 1.0                                   
 Copyright   : Don't Touch. Hehe                                          
 Description : Error Estimation for the Lorenz System                
 ============================================================================   
*/



#ifndef READ_FILE_H_
#define READ_FILE_H_

#include <stdio.h>
#include <stdlib.h>
#include "struct_def.h"
#include <string.h> //Required for the strtok() function      

Universe read_eqn(int argc, char* argv[]);
void read_input(int argc, char* argv[], Universe eqn, Galaxy* u_coarse,
               Galaxy* u_fine, Is_it *reduced);

#endif

