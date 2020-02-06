/*------------------------------------------------------------------*/
/* XFLOW: A discontinuous Galerkin finite element software library. */
/*                                                                  */
/*                    Copyright 2008-2015                           */
/*                 The University of Michigan                       */
/*                    All rights reserved                           */
/*      Contact:  Krzysztof J. Fidkowski, kfid@alum.mit.edu         */
/*                                                                  */
/*                    Copyright  2007-2008                          */
/*                   Krzysztof J. Fidkowski                         */
/*                                                                  */
/* This library is intended to be useful but is distributed without */
/* any warranty, not even merchantability or fitness for a          */
/* particular purpose.  It is free software: you can redistribute   */
/* it and/or modify it under the terms of the GNU Lesser General    */
/* Public License (LGPLv3).                                         */
/*                                                                  */
/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free       */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.        */
/*------------------------------------------------------------------*/



#ifndef _xf_ResidualNewDG_h
#define _xf_ResidualNewDG_h 1

/*
  FILE:  xf_ResidualNewDG.h

  This file contains headers for NewDG-specific residual-calculation
  functions

*/

#include "xf_MeshStruct.h"
#include "xf_DataStruct.h"
#include "xf_AllStruct.h"



extern int
yk_CalculateResidual_NewDG(xf_All *All, xf_VectorGroup *UG, xf_VectorGroup *RG,
                        xf_JacobianMatrix *R_U, xf_SolverData *SolverData);


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateResidual_NewDG
extern int
xf_CalculateResidual_NewDG(xf_All *All, xf_VectorGroup *UG, xf_VectorGroup *RG,
                        xf_JacobianMatrix *R_U, xf_SolverData *SolverData);
/*

PURPOSE:

  NewDG-specific version of residual calculation.  See full description
  in xf_Residual.h.

INPUTS:

  All: All structure
  UG : state vector group

OUTPUTS:

  R_U : Jacobian matrix
  RG : Residual vector group

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototpye: xf_BoundaryFluxIntegral_NewDG
extern int
xf_BoundaryFluxIntegral_NewDG(xf_All *All, xf_VectorGroup *UG,
                           xf_OutputEvalData *OutputEval,
                           int nBFG, int *BFGs);
/*

PURPOSE:

  NewDG-specific version of the boundary flux calculation function.  See
  xf_Residual.h for a full description.

INPUTS:

  All: All structure
  UG : state vector group
  OutputEval : output evaluation structure (see xf_SolverStruct.h).
  nBFG : number of boundary face groups over which to integrate
  BFGs : numbers of boundary face groups for the integration

OUTPUTS:

  OutputEval: variables in this structure are calculated

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ElemConnectivity_NewDG
extern int
xf_ElemConnectivity_NewDG(xf_All *All, xf_VectorGroup *UG, xf_Vector *R);
/*
PURPOSE:

  Calculates a real-valued directed connectivity magnitudes for
  interior faces of the discretization.

  DG: These magnitudes are Frobenius norms of the face flux linearizations.

  This may be called with a parallelized mesh.

INPUTS:

  All : All structure
  UG : State vector group

OUTPUTS:

  R : Vector of real-valued weights for each face

RETURNS:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PingResidual_NewDG
extern int
xf_PingResidual_NewDG(xf_All *All, xf_VectorGroup *UG, real eps, real tol);
/*
PURPOSE:

  Pings NewDG residual Jacobian matrix (not Schur complement)

INPUTS:

  All: All structure
  UG: state vector group
  ep  : epsilon used to perturb state
  tol : if > 0, this function will exit with error if
        the ping fails by tol*ep^2

OUTPUTS:

  pings output to stdout

RETURN:

  Error code
*/


#endif // end ifndef _xf_ResidualNewDG_h
