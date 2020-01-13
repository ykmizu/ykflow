

#include "ks_adjoint_MATVEC.h"

void ks_adjoint_MATVEC(yk_PrimalSolver *ykflow, Multiverse *multiquation,
		       Cluster *primal, double *R, double *Rnew,
		       Is_it *reduced){
  //------------------------------------//-------------------------------------
  // Variables here                     // Comments section                  
  //------------------------------------//-------------------------------------
  int i, j;                             // Intialization for iteration
  double wCon;
  double aveSpace;
  int left, right;
  int numFiles = primal->self->time.count/reduced->nTimeSegs;
  double timeFinal = primal->self->time.globalT_f;
  double timeInitial = primal->self->time.globalT_0;
  int systemSize;
  int *systemIndex;
  double denom;
  double timelength = (timeFinal-timeInitial)/reduced->nTimeSegs;
  double *vpArray;
  double *wpArray;
  Vec vp;
  Vec f;
  Vec wp;
  Universe _eqnOfInterest;  //Current universe to look at given Reduced info
  Galaxy *_solOfInterest; 
  Galaxy *_stateOfInterest;
  Cluster *Psi2 = (Cluster *) malloc (sizeof(Cluster));
  Cluster *Psi1 = (Cluster *) malloc (sizeof(Cluster));
  //---------------------------------------------------------------------------
  //Initialization
  //---------------------------------------------------------------------------
  if (reduced->reducedSolution == 1){
    _eqnOfInterest = multiquation->equationReduced;
    _solOfInterest = primal->reduced;
  }else{
    _eqnOfInterest = multiquation->equationLSS; //primal can be 
    _solOfInterest = primal->self;
  }
  if (reduced->hrom == 1)
    _stateOfInterest = primal->stateFull;
  else if (reduced->hrom == 0)
    _stateOfInterest = primal->self;
  //  systemSize = _eqnOfInterest.numStates*_solOfInterest->space.node.count;
  systemSize = _eqnOfInterest.numStates; //since dx = 0 for the LSS equations
  systemIndex = (int *) malloc (systemSize*sizeof(int));
  for (i=0; i<systemSize; i++)
    systemIndex[i] = i;
  vpArray = (double *) malloc (systemSize*sizeof(double));
  wpArray = (double *) malloc (systemSize*sizeof(double));;
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &vp);                
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &f);                 
  VecCreateSeq(PETSC_COMM_SELF, systemSize, &wp);
  //---------------------------------------------------------------------------
  // Implementation
  //---------------------------------------------------------------------------
  createAdjointCluster(ykflow, multiquation, primal, Psi2, reduced);
  createAdjointCluster(ykflow, multiquation, primal, Psi1, reduced);
  for (i=0; i<reduced->nTimeSegs; i++){
    //-------------------------------------------------------------------------
    // Time Segment initialization
    //-------------------------------------------------------------------------
    aveSpace =0;
    left = numFiles*i;     //Refers to the left node of the time segment
    right = numFiles*(i+1); //Refers to the right node of the time segment
    ks_copyUtype(_solOfInterest, Psi2->self); //initialize the adjoints
    Psi2->self->time.t_0 = timelength*i; //Refers tot he time t of the left no
    Psi2->self->time.t_f = timelength*(i+1);
    Psi2->self->basis.p = 0;
    Psi2->self->basis.nodes = 1;
    Psi2->self->space.elem.count = 1;
    Psi2->self->space.node.count = 1;
    sprintf(Psi2->clusterId, "psi2_seg_%d", i);
    sprintf(Psi2->self->utypeId, "psi2_seg_%d", i);
    createSystem(_eqnOfInterest, Psi2->self);
    ks_copyUtype(_solOfInterest, Psi1->self); //initialize the adjoints
    Psi1->self->time.t_0 = timelength*i; //Refers to the time t of the left
    Psi1->self->time.t_f = timelength*(i+1); //Refers to the time t of the rig
    Psi1->self->basis.p = 0;
    Psi1->self->basis.nodes = 1;
    Psi1->self->space.elem.count = 1;
    Psi1->self->space.node.count = 1;
    sprintf(Psi1->clusterId, "psi1_seg_%d", i);
    sprintf(Psi1->self->utypeId, "psi1_seg_%d", i);
    createSystem(_eqnOfInterest, Psi1->self);
    //-------------------------------------------------------------------------
    // Tangent Solve
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|              Tangent Solver: Segment %d                 |\n", i);
    printf("----------------------------------------------------------\n");
    getSegArray(_eqnOfInterest, Psi2->self, R, i, vp);
    ks_function(ykflow, multiquation,Psi2, reduced, f, left); //left +1 or lef
    projection(_eqnOfInterest, f, vp);
    vec2Array(_eqnOfInterest, Psi2->self, vp);
    ks_printSolution(_eqnOfInterest, Psi2->self, left);
    ks_ApplyTimeScheme(ykflow, multiquation, Psi2, reduced);
    //-------------------------------------------------------------------------
    // Lagrange Solve
    //-------------------------------------------------------------------------
    printf("----------------------------------------------------------\n");
    printf("|             Lagrange Solver: Segment %d                 |\n", i);
    printf("----------------------------------------------------------\n");
    if (i == reduced->nTimeSegs-1)
      VecZeroEntries(wp);
    else
      getSegArray(_eqnOfInterest, Psi1->self, R, reduced->nTimeSegs+i, wp);
    ks_function(ykflow, multiquation, Psi1, reduced, f, right);
    projection(_eqnOfInterest, f, wp);
    ykflow->aveSpatialOutput(ykflow, multiquation, _stateOfInterest, &aveSpace,
			     right);
    VecNorm(f, NORM_2, &denom);
    wCon = -(-primal->self->beta*(aveSpace-_stateOfInterest->j_bar)
	     /(primal->self->time.globalT_f*denom*denom));
    VecAXPY(wp, wCon, f);
    vec2Array(_eqnOfInterest, Psi1->self, wp);
    ks_printSolution(_eqnOfInterest, Psi1->self, right);
    ks_ApplyTimeScheme(ykflow, multiquation, Psi1, reduced);
    //-------------------------------------------------------------------------
    // Residual Calculation for Rw
    //-------------------------------------------------------------------------
    if (i<reduced->nTimeSegs - 1){
      ks_readSolution(_eqnOfInterest, Psi2->self, right);
      array2Vec(_eqnOfInterest, Psi2->self, vp);
      ks_function(ykflow, multiquation, Psi2, reduced, f, right);
      projection(_eqnOfInterest, f, vp);
      VecGetValues(vp, systemSize, systemIndex, vpArray);
      for (j=0; j<systemSize; j++)
        Rnew[(reduced->nTimeSegs+i)*systemSize+j]=
	  R[(i+1)*systemSize+j] - vpArray[j];
    }
    //-------------------------------------------------------------------------
    // Residual Calculation for Rv
    //-------------------------------------------------------------------------
    ks_readSolution(_eqnOfInterest, Psi1->self, left);
    array2Vec(_eqnOfInterest, Psi1->self, wp);
    ks_function(ykflow, multiquation, Psi1, reduced,  f, left);
    projection(_eqnOfInterest, f, wp);
        VecGetValues(wp, systemSize, systemIndex, wpArray);
    for (j=0; j<systemSize; j++)
      if (i==0)
        Rnew[i*systemSize+j] = -wpArray[j];
      else
        Rnew[i*systemSize+j] =
	  R[(reduced->nTimeSegs+i-1)*systemSize+j]-wpArray[j];
    destroySystem(_eqnOfInterest, Psi2->self);
    destroySystem(_eqnOfInterest, Psi1->self);
  }
  //---------------------------------------------------------------------------
  // Destroy and Clear out memory
  //---------------------------------------------------------------------------
  destroyAdjointCluster(ykflow, multiquation, primal, Psi2, reduced);
  destroyAdjointCluster(ykflow, multiquation, primal, Psi1, reduced);
  VecDestroy(&vp);
  VecDestroy(&wp);
  VecDestroy(&f);
  free(systemIndex);
  free(vpArray);
  free(wpArray);
  free(Psi2);
  free(Psi1);
}
