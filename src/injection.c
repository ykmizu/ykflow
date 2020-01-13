
#include "injection.h"

void inject(yk_PrimalSolver *ykflow, Multiverse *multiquation, Galaxy *primal,
	    Galaxy *primalFine, Is_it *reduced){
  //Just takes into account of higher order interpolation
  int i, j, k;  //initialization for iteration
  int count=0; //placement
  double xi_;
  primalFine->time.node = primal->time.node;
  for (i=0; i<primal->space.elem.count; i++){    
    for (k=0; k<primalFine->basis.p + 1; k++){
      xi_ = xiL(primalFine->basis.p , k);     
      for (j=0; j<multiquation->equation.numStates; j++){
        primalFine->solution[j].array[count] = U_xi(primal, j, i, xi_);
        count++;
	
      }
    }
  }
}

/* void injectAll(Universe equation, Galaxy *primal, Galaxy *primalFine){ */
/*   int i; */
/*   for (i=0; i<primal->time.count+1; i++){ */
/*     ks_readSolution(equation, primal, i); */
/*     inject(equation, primal, primalFine); */
/*     ks_printSolution(equation, primalFine, i); */
/*   } */
/* } */
