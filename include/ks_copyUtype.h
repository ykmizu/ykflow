

#ifndef KS_COPYUTYPE_H_
#define KS_COPYUTYPE_H_

#include <stdio.h>
#include <stdlib.h>
#include "struct_def.h"
#include "system_initialization.h"
#include <stdarg.h>
#include <string.h>

void ks_copyUtype(Galaxy *U, Galaxy *Ucopy);

void ks_copySolutions(Universe equation, Galaxy *U, Galaxy *Ucopy);

void ks_setSolutionZero(Universe equation, Galaxy *U);
#endif
