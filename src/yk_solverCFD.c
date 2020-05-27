
#include "yk_solverCFD.h"


yk_PrimalSolver* new_Ykflow(){
  //allocationg memory
  yk_PrimalSolver *ykflow =
    (yk_PrimalSolver *) malloc (sizeof(yk_PrimalSolver));
  //point to itself to desinate the base class object
  ykflow->solver = ykflow;
  ykflow->Function = yk_ykflow_function;
  ykflow->dFunctiondu = yk_ykflow_dfunctiondu;
  ykflow->Residual = yk_ykflow_totalResidual;
  //  ykflow->FomName = yk_fomBuffer;
  ykflow->boundarySeeds = yk_findBoundarySeeds1D;
  /* ykflow->adjacentElems = yk_findAdjacentElems1D; */
  ykflow->aveSpatialOutput = spatialOutputAverage;
  ykflow->spatialdJdU = ks_dJdU;
  ykflow->anyFunction = yk_any1D_function;
  /* ykflow->multiplyByInvMass = yk_MultInvMass; */
  ykflow->injectH2h = inject;
  ykflow->delete = delete_Ykflow;
  return ykflow;
}

void delete_Ykflow(yk_PrimalSolver *ykflow){
  free(ykflow);
}
