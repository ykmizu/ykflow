#include "yk_ResidualDGTime.h"

void yk_ResidualDGTime(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Cluster *primal, Vec residual, Is_it *reduced,
		       int timeSolveNum){

//timeSolverNum here will just be NULL for now...


  int r;  //time integration interpolation order
  Universe _equation = multiquation->equation;
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  for (i=0; i<r+1; i++){ //There's goign to r+1 residuals value
    for (j=0; j<primal->quad.n; n++){
      R = ks_spatialStateResidual(_equation, primal->self, residual, reduced);
      state = U_xi(primal, 0, , primal->quad.x_i.array[n]);
      residual = -dphi(primal->basis.r, j, primal->quad.x_i.array[j])*state*
	primal->quad.w.array[j] + primal->self->time.dt/2.0*
	phi(primal->basis.r,j,primal->quad.x_i.array[j])*R*
	primal->quad.w.array[j];
      VecSetValue(R, index, residual, ADD_VALUES);
    }
    //-------------------------------------------------------------------------
    // Flux Time terms Upwind
    //-------------------------------------------------------------------------
    timePlus = U_xi(primal->self, elem-1, 1)*phi(primal->basis.r, i, 1)-
      U_xi(primal->self, elem, 1)*phi(primal->basis.r, i, -1);
    VecSetValue(timePlus, index, residual, ADD_VALUES);
  }
}
