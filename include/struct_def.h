
#ifndef STRUCT_DEF_H_
#define STRUCT_DEF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "petsc.h"
//-----------------------------------------------------------------------------
// Struct_def.h Contains all the tools for organization
//
// Yukiko Shimizu
// August 9, 2016
// Error Estimatiofor Chaotic Systems
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//Creates the Array structure so that I can call and size an array on the fly
//-----------------------------------------------------------------------------
typedef struct{
  int size;               //Size of the data
  double *array;          //Contains the actual data
}Array;

typedef struct{
  int count;
  int *array;
}intArray;


typedef struct{
  int xlength;            //Number of the rows for the data
  int ylength;            //Number of the columns for the data
  double **array2;        //Contains the actual data for 2 dimensional arrays
}Array2D;

void initArray(Array *a, int size);         //initialize data for Array
void initArrayInt(intArray *a, int size);
void initArray2D(Array2D *a, int m, int n); //initialize data for 2D Array
void delArray(Array *a);                    //delete data for Array
void delArrayInt(intArray *a);
void delArray2D(Array2D *a);                //delete data for Array

//-----------------------------------------------------------------------------
// Quadrature and Basis stuff
//-----------------------------------------------------------------------------
typedef struct {
  int n;                  //number of quad points to evaluate integral over
  Array x_i;              //the corresponding spatial value for each quad point
  Array w;                //the weights refering to each quadrature point
}GuassianQuad;

typedef struct {
  int p;                  //spatial interpolation number
  int nodes;
  char pId[1000];         //corresponding name for the designated p given
  Array2D phi;            //basis functions
  Array2D dphi;           //derivative of basis (Do I still use this function)
}BasisFunctions;

//-----------------------------------------------------------------------------
// Solution organization
//--------- -------------------------------------------------------------------

typedef struct{
  double dx;               //Spatial step for entire simulation (1 dim so far)
  Array x;
  intArray elem;        //Contains information for the elements
  intArray node;           //Contains information for the nodes
  double x_0;              //Initial left boundary spatial value
  double x_f;              //Final right boundary spatial value
}Mesh;

typedef struct{
  double dt;               //Time step for entire simulation
  int count;               //Number of time step elements in entire simulation
  int node;                //Refering to the time node the spatial is on
  double globalT_0;        //Initial time for the entire simulation
  double globalT_f;        //Final time for the entire simulation
  double t_0;              //Local initial time for the currenttime window
  double t_f;              //Local final time for the the current time window
}Time;

//-----------------------------------------------------------------------------
// General Program Solution organization: Contains solution for 1 time node
//-----------------------------------------------------------------------------
typedef struct { //the command module
  char utypeId[1000];       //Id name for Utype such as coarse vs fine
  char time_method_name[1000];
  int time_method;          //Type of time solver to use
  int typeSolution;         //full order, reduced, hyper-reduced
  double burnT_f;           //Time length to burn to get on the attractor
  double j_bar;             //Statistical time average output
  Vec djdu;                 //djdu part
  double error_estimate;    //The output of interest
  Time time;                //Contains structure for time variables/data
  Mesh space;               //Contains structure for space variables/data
  BasisFunctions basis;     //Basis functions neccesary for DG methods
  GuassianQuad quad;        //Quadrature to solve integrals
  char id[1000];            //Name for this set of solutions
  Array* solution;
  int beta;                //if pde or ode requires any constant as part of eq
  Mat Mij;
  int systemSize;
  PetscInt *index;
}Galaxy; //Galaxies


typedef struct{
  double actualJbar;
  int typeLSS;
  double delJ_psi;
  double delJ_dpsi;
  double delJ_dpsi_H1;
  double tradDelJ;
  int dss;
  int hrom; //indiciates if we are dealing with a rom or hrom solution
  int reducedSolution;    //1 is yes 0 is no
  int nTimeSegs;
  int restart;
  int reducedLSS;
  int innerStart;
  int outerStart;
  int numParamSet;
  int numParams;
  double *params;
  double *paramsL;
  double *paramsH;
  double *dparams;
  double *paramMeshGrid;
  int nSampleNodes;  //number of sample nodes
  int nBasisFuncs;
  int nBasisTime;
  int nBasisFuncsRJ;
  int nSnapshotsRJ;
  Mesh reducedMesh;       //sample nodes+ nodes within the same element
  Mat rOBState;           //if 1, then rOBState, A, and B should be filled up
  Mat ST_rOBState;
  Mat *ST_rOBState_i;
  Mat rOBStateRed;        //number of rows is the sixe of reducedMesh
  Mat rOBStateBar;        //Number of rows is the size of primalApprox
  Mat rOBResidual;
  Mat rOBJacobian;
  Mat Z;
  Mat Zp;
  Mat Zu;
  Mat Zmult;
  Mat A;   //Offline Matrix
  Mat B;   //Ofline Matrix
}Is_it;

void initIs_it(Is_it *reduced);

void initUtype(Galaxy *primal);


typedef struct{
  char clusterId[1000];
  Galaxy *self;
  Galaxy *reduced;
  Galaxy *state;
  Galaxy *stateFull; //full state when found from rOBState
  Galaxy *adjointFull;
}Cluster;

typedef struct {
  char nameEqn[10000];
  double *c;
  int numStates;
  int numParams;
  int paramIndex;
  double paramValue;
}Universe;

typedef struct{
  Universe equation;
  Universe equationLSS;
  Universe equationReduced;
}Multiverse;

#endif