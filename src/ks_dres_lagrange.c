
#include "ks_dres_lagrange.h"

void ks_dres_lagrange(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		      Cluster *Psi1, Mat dRdU, Is_it *reduced){
  int timeNode = Psi1->self->time.node;
  ks_dfunctiondu(ykflow, multiquation, Psi1, reduced, dRdU, timeNode); 
  MatScale(dRdU, -1);
  MatHermitianTranspose(dRdU, MAT_INPLACE_MATRIX, &dRdU);
}
