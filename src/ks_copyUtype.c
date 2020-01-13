
#include "ks_copyUtype.h"

void ks_copyUtype(Galaxy *U, Galaxy *Ucopy){
  strcpy(Ucopy->utypeId, U->utypeId);
  Ucopy->burnT_f = U->burnT_f;
  Ucopy->time_method = U->time_method;
  Ucopy->time.dt = U->time.dt;
  Ucopy->time.globalT_0 = U->time.globalT_0;
  Ucopy->time.globalT_f = U->time.globalT_f;
  Ucopy->time.t_0 = U->time.t_0;
  Ucopy->time.t_f = U->time.t_f;
  Ucopy->time.count = U->time.count;
  Ucopy->space.elem.count = U->space.elem.count;
  Ucopy->space.node.count = U->space.node.count;
  Ucopy->space.dx = U->space.dx;
  Ucopy->space.x_0 = U->space.x_0;
  Ucopy->space.x_f = U->space.x_f;
  Ucopy->basis.p = U->basis.p;
  Ucopy->basis.nodes = U->basis.nodes;
  Ucopy->quad.n = U->quad.n;
  
}

void ks_copySolutions(Universe equation, Galaxy *U, Galaxy *Ucopy){
  int i, j;
  Ucopy->time.node = U->time.node;
  for (i=0; i<equation.numStates; i++)
    for (j=0; j<U->space.node.count; j++)
      Ucopy->solution[i].array[j] = U->solution[i].array[j];
}

void ks_setSolutionZero(Universe equation, Galaxy *U){
  int i, j;
  for (i=0; i<equation.numStates; i++)
    for (j=0; j<U->space.node.count; j++)
      U->solution[i].array[j] = 0;    
}

