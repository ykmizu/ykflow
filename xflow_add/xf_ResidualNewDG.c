
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



/*
  FILE:  xf_ResidualNewDG.c

  "Clean" DG implementation.  This file will eventually replace
  xf_ResidualDG.c.

*/
#include <stdlib.h>
#include "xf_String.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_MeshTools.h"
#include "xf_Solver.h"
#include "xf_ResidualCommon.h"
#include "xf_ResidualEqnSet.h"
#include "xf_ResidualStab.h"
#include "xf_State.h"
#include "xf_MPI.h"


/* Storage for temporary vectors and matrices (work space) */
typedef struct
{
  int Asize;
  real *A;

  int Psize;
  int *P;
}
xf_DGWorkSpace;


/* DG storage structure used to assist in the calculation of residual
   and Jacobian contributions. */
typedef struct
{
  /*-------------------------*/
  /*   General information   */
  /*-------------------------*/

  // pointer to vector group state
  xf_VectorGroup *UG;

  // pointer to vector group residual
  xf_VectorGroup *RG;

  /*------------------------------*/
  /*  Residual storage structure  */
  /*------------------------------*/

  xf_ResidualStorage *RS;


  /*-------------------------------*/
  /*  State and Auxiliary vectors  */
  /*-------------------------------*/

  // shortcuts to state vector
  xf_Vector *U;

  /*----------------------*/
  /*   Residual vectors   */
  /*----------------------*/

  // shortcuts to vector
  xf_Vector *R;

  /*---------------------*/
  /*  Jacobian Matrices  */
  /*---------------------*/

  xf_JacobianMatrix *R_U;   // pointer to Jacobian matrix

  xf_GenMatrix RE_UE;       // self residual w.r.t self state
  xf_GenMatrix RF_UF[2][2]; // left/right residual w.r.t left/right state

  /*-------------*/
  /*  Workspace  */
  /*-------------*/

  // Temporary storage vectors and matrices
  xf_DGWorkSpace WS;

}
xf_DGStorage;



/******************************************************************/
//   FUNCTION Definition: xf_InitDGStorage
static void
xf_InitDGStorage(xf_DGStorage *DS)
{
/*
  PURPOSE: Initializes DG storage structure
  INPUTS : DS = structure to initialize
  OUTPUTS: DS, with data initialized to 0/NULL
  RETURN : None
*/
  int i;

  DS->UG = NULL;
  DS->RG = NULL;
  DS->RS = NULL;

  // States
  DS->U = NULL;

  // Residuals
  DS->R  = NULL;

  // initialize Jacobian matrices
  DS->R_U = NULL;
  xf_InitGenMatrix(&DS->RE_UE);
  for (i=0; i<4; i++) xf_InitGenMatrix(&DS->RF_UF[i/2][i%2]);

  // Workspace
  DS->WS.Asize = 0;
  DS->WS.A     = NULL;
  DS->WS.Psize = 0;
  DS->WS.P     = NULL;

}


/******************************************************************/
//   FUNCTION Definition: xf_AllocateDGStorage
static int
xf_AllocateDGStorage(xf_All *All, xf_VectorGroup *UG, xf_VectorGroup *RG,
                     xf_SolverData *SolverData, xf_JacobianMatrix *R_U,
                     xf_DGStorage **pDS)
{
/*
PURPOSE:

  Allocates DG storage structure, including room for self

INPUTS:

  All        : all structure
  UG         : state vector group
  RG         : residual vector group
  SolverData : solver data structure

OUTPUTS:

  (*pDS)     : allocated DG storage structure

RETURN:

  Error code

*/

  int ierr;
  int i;
  xf_DGStorage *DS;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // allocate self
  xf_Call(xf_Alloc( (void **) pDS, 1, sizeof(xf_DGStorage)));
  DS = (*pDS);

  // Initialize all entries to NULL;
  xf_InitDGStorage(DS);

  // state vector group
  DS->UG = UG;

  // residual vector group
  DS->RG = RG;

  // Allocate residual storage
  xf_Call(xf_AllocateResidualStorage(All, SolverData, &DS->RS));

  // set pointer to state vector
  xf_Call(xf_GetVectorFromGroup(UG, xfe_VectorRoleElemState , &DS->U));

  // set pointer to residual vector
  if (RG != NULL)
    xf_Call(xf_GetVectorFromGroup(RG, xfe_VectorRoleElemState, &DS->R));

  // Jacobian matrix pointer
  DS->R_U = R_U;

  // Local Jacobian matrices will just be pointers, so do not resize
  DS->RE_UE.DoResize = xfe_False;
  for (i=0; i<4; i++) DS->RF_UF[i/2][i%2].DoResize = xfe_False;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyDGStorage
static int
xf_DestroyDGStorage(xf_DGStorage *DS)
{
/*
  PURPOSE: Destroys DG structure, including self
  INPUTS : DS = DG structure
  OUTPUTS: None
  RETURN : Error code
*/

  int ierr;
  int i;

  // destroy Jacobian matrix meta-data
  xf_DestroyGenMatrix(&DS->RE_UE, xfe_False);
  for (i=0; i<4; i++) xf_DestroyGenMatrix(&DS->RF_UF[i/2][i%2], xfe_False);

  // Destroy residual storage
  xf_Call(xf_DestroyResidualStorage(DS->RS));

  // destroy workspace
  xf_Release( (void *) DS->WS.A);
  xf_Release( (void *) DS->WS.P);

  // destroy self
  xf_Release( (void *) DS);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PrepareElemIntDG
static int
xf_PrepareElemIntDG(xf_All *All, int egrp, int elem, xf_DGStorage *DS)
{
/*
PURPOSE:

  Prepares DG storage structure DS for calculating residuals/Jacobians
  associated with element egrp, elem.

INPUTS:

  All        : all structure
  egrp, elem : element in question
  DS         : DG storage structure

OUTPUTS: None
RETURN : Error code

*/

  int i;
  xf_ResidualStorage *RS;
  xf_GenMatrix **GM;

  // pointer to residual storage structure
  RS = DS->RS;

  // reset minimum dt (for artificial time stepping)
  xf_ResetArtificialDtInfo(&RS->ArtificialDt);

  // make local data look stale
  xf_ForceReCalcResidualStorage(RS);

  // point vectors in RS to appropriate states/residuals
  RS->VU[xfe_ElemStateU] = DS->U;
  RS->VR[xfe_ElemStateU] = DS->R;

  // point residuals in RS to correct location
  for (i=0; i<xfe_ElemStateLast; i++) RS->R[i] = NULL; // first set all to NULL
  RS->R[xfe_ElemStateU] = DS->R->GenArray[egrp].rValue[elem];

  // point linearization request in RS to R_U
  GM = RS->R_V->GM[xfe_ElemStateU];
  if ((DS->R_U != NULL) && (DS->R_U->Value != NULL)){
    RS->NeedLin = xfe_True;
    DS->RE_UE.A = DS->R_U->Value[egrp][elem][0];
    GM[xfe_ElemStateU] = &DS->RE_UE;
  }
  else{
    RS->NeedLin = xfe_False;
    GM[xfe_ElemStateU] = NULL;
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElemDG
static int
xf_CalculateResidualElemDG(xf_All *All, int egrp, int elem, xf_DGStorage *DS,
                           xf_SolverData *SolverData)
{
/*
PURPOSE:

  Calculates an element's contribution to the DG residuals and
  Jacobians (if requested).

INPUTS:

  All        : all structure
  egrp, elem : element number
  DS         : DG storage structure
  SolverData : solver data structure

OUTPUTS:

  Residuals and Jacobians are updated

RETURN:

  Error code
*/
  int ierr;

  // Prepare DG storage for element integration
  xf_Call(xf_PrepareElemIntDG(All, egrp, elem, DS));

  // RU{i} += - int_k gradphi{i,d} * (F+G){d}  (implied sum)
  xf_Call(xf_AddResTermElemInt(All, SolverData, egrp, elem, DS->RS, xfe_ElemIntGPhiDotFlux, -1.0));

  // RU{i} += int_k phi{i} * s  (s = source)
  xf_Call(xf_AddResTermElemInt(All, SolverData, egrp, elem, DS->RS, xfe_ElemIntPhiSource, 1.0));

  if ((DS->RS->IncludeUnsteady) || (DS->RS->ArtificialDt.Active)){
    // RU{i} += (1/dt + 1/dtartificial) * int_k phi{i} * u   (u = element interior state)
    xf_Call(xf_AddResTermElemInt(All, SolverData, egrp, elem, DS->RS, xfe_ElemIntPhiStateUnsteady, 1.0));
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PrepareFaceIntDG
static int
xf_PrepareFaceIntDG(xf_All *All, int FGroup, int FNumber, xf_DGStorage *DS)
{
/*
PURPOSE:

  Prepares DG storage structure DS for calculating residuals/Jacobians
  associated with face FGroup, FNumber.

INPUTS:

  All             : all structure
  SolverData      : solver data structure
  FGroup, FNumber : face in question
  DS              : DG storage structure

OUTPUTS: None
RETURN : Error code

*/


  int ierr, f, i;
  xf_ResidualStorage *RS;
  enum xfe_Bool OnBC = xfe_False;
  enum xfe_Bool OnBCL = xfe_False, OnBCR = xfe_False;
  enum xfe_FaceIntVecType FIV[2] = {xfe_FaceIntVecLeft, xfe_FaceIntVecRight};
  real *RF_UF[2][2];
  real *RL, *RR;
  real *RL_UL, *RL_UR, *RR_UL, *RR_UR;
  xf_FaceInfo FI;
  xf_GenMatrix ***GM;

  // information about the face
  xf_Call(xf_GetFaceInfo(All->Mesh, FGroup, FNumber, &FI));
  OnBC = (FI.FaceType == xfe_FaceBoundary);

  // pointer to residual storage structure
  RS = DS->RS;

  // reset minimum dt (for artificial time stepping)
  xf_ResetArtificialDtInfo(&RS->ArtificialDt);

  // make local data look stale
  xf_ForceReCalcResidualStorage(RS);

  // reset state/residual pointers to NULL
  for (f=0; f<xfe_FaceIntVecLast; f++)  // loop over face integration groups
    for (i=0; i<xf_FaceIntVec2nState(f); i++)
      RS->VUF[f][i] = RS->VRF[f][i] = NULL;

  // flags indicating whether L or R elems are on BC (i.e. do not exist)
  OnBCL = (FI.egrpL < 0);
  OnBCR = ((FI.egrpR < 0) || OnBC);

  // point to correct state/residual vectors
  if (!OnBCL){
    RS->VUF[xfe_FaceIntVecLeft][xfe_ElemStateU] = DS->U;
    RS->VRF[xfe_FaceIntVecLeft][xfe_ElemStateU] = DS->R;
  }
  if (!OnBCR){
    RS->VUF[xfe_FaceIntVecRight][xfe_ElemStateU] = DS->U;
    RS->VRF[xfe_FaceIntVecRight][xfe_ElemStateU] = DS->R;
  }

  // point residuals to correct locations
  RL = RR = NULL;
  if (DS->R != NULL){
    RL = (OnBCL) ? NULL : DS->R->GenArray[FI.egrpL].rValue[FI.elemL];
    RR = (OnBCR) ? NULL : DS->R->GenArray[FI.egrpR].rValue[FI.elemR];
  }
  RS->RF[xfe_FaceIntVecLeft ][xfe_ElemStateU] = RL;
  RS->RF[xfe_FaceIntVecRight][xfe_ElemStateU] = RR;

  // point to residual linearizations
  if ((DS->R_U != NULL) && (DS->R_U->Value != NULL)){
    RS->NeedLin = xfe_True;
    RL_UL = RL_UR = RR_UL = RR_UR = NULL;
    if (!OnBCL) RL_UL = DS->R_U->Value[FI.egrpL][FI.elemL][0];
    if (!OnBCR) RR_UR = DS->R_U->Value[FI.egrpR][FI.elemR][0];
    if ((!OnBCL) && (!OnBCR)){
      RL_UR = DS->R_U->Value[FI.egrpL][FI.elemL][1+FI.faceL];
      RR_UL = DS->R_U->Value[FI.egrpR][FI.elemR][1+FI.faceR];
    }

    // more convenient for access of linearizations
    RF_UF[0][0] = RL_UL;
    RF_UF[0][1] = RL_UR;
    RF_UF[1][0] = RR_UL;
    RF_UF[1][1] = RR_UR;

    // dependence of left/right elem residuals on left/right elem states
    for (i=0; i<4; i++){
      GM = RS->RF_VF[FIV[i/2]][FIV[i%2]]->GM;
      GM[xfe_ElemStateU][xfe_ElemStateU] = NULL;
    }
    if ((DS->R_U != NULL) && (DS->R_U->Value != NULL)){
      for (i=0; i< (OnBC ? 1 : 4); i++){
        GM = RS->RF_VF[FIV[i/2]][FIV[i%2]]->GM;
        DS->RF_UF[i/2][i%2].A = RF_UF[i/2][i%2];
        GM[xfe_ElemStateU][xfe_ElemStateU] = (RF_UF[i/2][i%2] == NULL) ? NULL : &DS->RF_UF[i/2][i%2];
      }
    }
  }
  else{
    RS->NeedLin = xfe_False;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualFaceDG
static int
xf_CalculateResidualFaceDG(xf_All *All, int FGroup, int FNumber, xf_DGStorage *DS,
                           xf_SolverData *SolverData)
{
/*
PURPOSE:

  Calculates an element's contribution to the DG residuals and
  Jacobians (if requested).

INPUTS:

  All        : all structure
  egrp, elem : element number
  DS         : DG storage structure
  SolverData : solver data structure

OUTPUTS:

  Residuals and Jacobians are updated

RETURN:

  Error code

*/

  int ierr;

  // Reset structures and prepare for face integration
  xf_Call(xf_PrepareFaceIntDG(All, FGroup, FNumber, DS));

  // RF{i} += int_face phi{i,d} * (F+G){d}  (implied sum)
  xf_Call(xf_AddResTermFaceInt(All, SolverData, FGroup, FNumber, DS->RS,
                               xfe_FaceIntPhiTraceFluxHatLR, 1.0));

  return xf_OK;
}



int yk_CalculateResidual_NewDG(xf_All *All, xf_VectorGroup *UG,
			       xf_VectorGroup *RG, xf_JacobianMatrix *R_U,
			       xf_SolverData *SolverData){
  int ierr;
  int i_nface;
  int i;
  int egrp, elem;
  int ifgrp, iface, iiface;
  int FGroup;
  int nIFaceRegular, nIFaceTot;
  enum xfe_Bool Physical;
  xf_VectorGroup *SG0, *SG;
  xf_DGStorage *DS=NULL;
  xf_Mesh *Mesh;
  int *faceHash = (int *) malloc (All->Mesh->nIFace*sizeof(int));

  Mesh = All->Mesh;
  // begin halo exchange of vector group UG
  xf_Call(xf_HaloExchangeVectorGroupBegin(UG, -1));

  // source vectors
  SG0 = (SolverData != NULL) ? SolverData->SG0 : NULL;
  SG = (SolverData != NULL) ? SolverData->SG : NULL;

  // initialize residual to zero or SG if present
  xf_Call(xf_SetZeroVectorGroup(RG));
  if (SG0 != NULL) xf_Call(xf_SetVectorGroup(SG0, xfe_Add, RG));
  if (SG  != NULL) xf_Call(xf_SetVectorGroup(SG,  xfe_Add, RG));

  // initialize residual Jacobian, R_U, to zero
  if ((R_U != NULL) && (R_U->Value != NULL)){
    xf_Call(xf_SetZeroJacobian(R_U));
  }

  // determine memory requirements and allocate local storage for linearization
  xf_Call(xf_AllocateDGStorage(All, UG, RG, SolverData, R_U, &DS));

  // physical error flag check
  Physical = xfe_True;

  // do we need artificial stabilization?
  if (SolverData != NULL){
    ierr = xf_CalculateStabilization(All, DS->U, DS->RS->NeedLin, SolverData);
    if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
    else if (ierr != xf_OK) return xf_Error(ierr);
    if (SolverData->StabData.StabType != xfe_StabilizationNone)
      DS->RS->StabData = &SolverData->StabData;
  }
  //Now need to loop through the elements this time
  //Lets jsut keep this loop here for now

  nIFaceRegular = (Mesh->ParallelInfo != NULL) ? Mesh->ParallelInfo->nIFaceRegular : -1;
  nIFaceTot = Mesh->nIFace;

  for (i =0 ; i<nIFaceTot; i++)
    faceHash[i] = 0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){;
      /* //Add the main element contributions from weak form */
      ierr = xf_CalculateResidualElemDG(All, egrp, elem, DS, SolverData);
      if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
      else if (ierr != xf_OK) return xf_Error(ierr);
      //Add the fluxes contributions from the boundary,
      //Check if the element has a boundary face
        // residuals on self boundary faces

      for (i_nface = 0; i_nface<Mesh->ElemGroup[egrp].nFace[elem]; i_nface++){
      	ifgrp = Mesh->ElemGroup[egrp].Face[elem][i_nface].Group;
      	iface = Mesh->ElemGroup[egrp].Face[elem][i_nface].Number;


	if (ifgrp >0 ){
	  if (strncmp(Mesh->BFaceGroup[ifgrp-1].Title, "ZeroMeasure", 11) == 0) continue;
          xf_GetFaceGroup(Mesh, xfe_FaceBoundary, ifgrp-1, xfe_False, &FGroup);
	  ierr = xf_CalculateResidualFaceDG(All, FGroup, iface, DS, SolverData);
          if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
          else if (ierr != xf_OK) return xf_Error(ierr);
	}else if (ifgrp == 0 && faceHash[iface]==0){ //Boundary groups
	  faceHash[iface] = 1;
	  if (iface == nIFaceRegular) // halo exchange must be completed at this point
	    xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
	  // continue if neither L nor R elems are on current proc
	  if (!xf_FullLinkageIFace(Mesh, iface)) continue;
	  xf_GetFaceGroup(Mesh, xfe_FaceInterior, -1, xfe_False, &FGroup);
	  ierr = xf_CalculateResidualFaceDG(All, FGroup, iface, DS, SolverData);
	  if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
	  else if (ierr != xf_OK) return xf_Error(ierr);
	}

	// Physical check
      	xf_Call(xf_MPI_Allreduce((int *) &Physical, 1, xfe_SizeInt, xfe_MPI_MIN));
      	if (!Physical){
      	  xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
      	  return xf_Error(xf_NON_PHYSICAL);
      	}

      }
    }
  }
      // Physical check
  xf_Call(xf_MPI_Allreduce((int *) &Physical, 1, xfe_SizeInt, xfe_MPI_MIN));
  if (!Physical){
    xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
    return xf_Error(xf_NON_PHYSICAL);
  }
  free(faceHash);
  // destroy local storage
  xf_Call(xf_DestroyDGStorage(DS));

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidual_NewDG
int
xf_CalculateResidual_NewDG(xf_All *All, xf_VectorGroup *UG, xf_VectorGroup *RG,
                           xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  int ierr;
  int egrp, elem;
  int ibfgrp, ibface, iiface;
  int FGroup;
  int nIFaceRegular, nIFaceTot;
  enum xfe_Bool Physical;
  xf_VectorGroup *SG0, *SG;
  xf_DGStorage *DS=NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  // begin halo exchange of vector group UG
  xf_Call(xf_HaloExchangeVectorGroupBegin(UG, -1));

  // source vectors
  SG0 = (SolverData != NULL) ? SolverData->SG0 : NULL;
  SG = (SolverData != NULL) ? SolverData->SG : NULL;

  // initialize residual to zero or SG if present
  xf_Call(xf_SetZeroVectorGroup(RG));
  if (SG0 != NULL) xf_Call(xf_SetVectorGroup(SG0, xfe_Add, RG));
  if (SG  != NULL) xf_Call(xf_SetVectorGroup(SG,  xfe_Add, RG));

  // initialize residual Jacobian, R_U, to zero
  if ((R_U != NULL) && (R_U->Value != NULL)){
    xf_Call(xf_SetZeroJacobian(R_U));
  }

  // determine memory requirements and allocate local storage for linearization
  xf_Call(xf_AllocateDGStorage(All, UG, RG, SolverData, R_U, &DS));

  // physical error flag check
  Physical = xfe_True;

  // do we need artificial stabilization?
  if (SolverData != NULL){
    ierr = xf_CalculateStabilization(All, DS->U, DS->RS->NeedLin, SolverData);
    if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
    else if (ierr != xf_OK) return xf_Error(ierr);
    if (SolverData->StabData.StabType != xfe_StabilizationNone)
      DS->RS->StabData = &SolverData->StabData;
  }

  // residuals on self elements
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ierr = xf_CalculateResidualElemDG(All, egrp, elem, DS, SolverData);
      if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
      else if (ierr != xf_OK) return xf_Error(ierr);
    } // elem
  } // egrp

  // residuals on self boundary faces
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    // skip zero-measure boundary groups
    if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) == 0) continue;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      xf_GetFaceGroup(Mesh, xfe_FaceBoundary, ibfgrp, xfe_False, &FGroup);
      ierr = xf_CalculateResidualFaceDG(All, FGroup, ibface, DS, SolverData);
      if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
      else if (ierr != xf_OK) return xf_Error(ierr);
    }
  } // ibfgrp

  // Physical check
  xf_Call(xf_MPI_Allreduce((int *) &Physical, 1, xfe_SizeInt, xfe_MPI_MIN));
  if (!Physical){
    xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
    return xf_Error(xf_NON_PHYSICAL);
  }

  /* residuals on interior faces */

  // regular faces are ones for which L and R elems are on current proc
  // total faces > self faces if some faces are next to halos (in parallel)
  nIFaceRegular = (Mesh->ParallelInfo != NULL) ? Mesh->ParallelInfo->nIFaceRegular : -1;
  nIFaceTot = Mesh->nIFace;
  if (Mesh->ParallelInfo != NULL) nIFaceTot += Mesh->ParallelInfo->nIFaceHalo;

  for (iiface=0; iiface < nIFaceTot; iiface++){
    if (iiface == nIFaceRegular){ // halo exchange must be completed at this point
      xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
    }
    // continue if neither L nor R elems are on current proc
    if (!xf_FullLinkageIFace(Mesh, iiface)) continue;

    xf_GetFaceGroup(Mesh, xfe_FaceInterior, -1, xfe_False, &FGroup);
    ierr = xf_CalculateResidualFaceDG(All, FGroup, iiface, DS, SolverData);
    if (ierr == xf_NON_PHYSICAL) Physical = xfe_False;
    else if (ierr != xf_OK) return xf_Error(ierr);
  } // iiface

    // Physical check
  xf_Call(xf_MPI_Allreduce((int *) &Physical, 1, xfe_SizeInt, xfe_MPI_MIN));
  if (!Physical){
    xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
    return xf_Error(xf_NON_PHYSICAL);
  }


  // destroy local storage
  xf_Call(xf_DestroyDGStorage(DS));

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BoundaryFluxIntegral_NewDG
int
xf_BoundaryFluxIntegral_NewDG(xf_All *All, xf_VectorGroup *UG,
                              xf_OutputEvalData *OutputEval,
                              int nBFG, int *BFGs)
{
  int ierr;
  int ibfg, ibfgrp, nbfgrp, ibface;
  int egrp, elem, face;
  int f, i;
  xf_BFace BFace;
  xf_VectorGroup *Value_UG;
  xf_DGStorage *DS;
  xf_Face *Face;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // sanity check
  if (OutputEval == NULL) return xf_Error(xf_INPUT_ERROR);

  // pull off output value linearization
  Value_UG = OutputEval->Value_UG;

  // allocate local storage for linearization; Value_UG takes place of residuals
  xf_Call(xf_AllocateDGStorage(All, UG, Value_UG, NULL, NULL, &DS));

  // store pointer to OutputEval structure in DS
  DS->RS->OutputEval = OutputEval;

  // Where do we put the linearization?
  // Into Value_UG, which is in DS as a residual .. so nothing else to set

  // Loop over boundary groups and faces
  nbfgrp = ((BFGs == NULL) ? Mesh->nBFaceGroup : nBFG);
  for (ibfg=0; ibfg<nbfgrp; ibfg++){

    ibfgrp = ((BFGs == NULL) ? ibfg : BFGs[ibfg]);

    // skip zero-measure boundary groups
    if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) == 0) continue;

    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){

      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];

      // element info
      egrp = BFace.ElemGroup;
      elem = BFace.Elem;
      face = BFace.Face;

      // face as seen by the element
      Face    = Mesh->ElemGroup[egrp].Face[elem] + face;

      // skip halo elements (do not want to double count)
      if (egrp >= Mesh->nElemGroup) continue;

      // prepare for face integration
      xf_Call(xf_PrepareFaceIntDG(All, Face->Group, Face->Number, DS));

      // specify where we want output linearizations placed
      DS->RS->NeedLin = (Value_UG != NULL);
      for (f=0; f<xfe_FaceIntVecLast; f++)
        for (i=0; i<xf_FaceIntVec2nState(f); i++) DS->RS->Output_VF[f][i] = DS->RS->RF[f][i];

      // J += int_face FluxWeights{k} * (Fhat + Ghat){d,k} * n{d}
      // note, xglob[FluxMoments[k]] is included when FluxMoments[k] >= 0
      // on boundaries, Fhat is calculated from ub; Ghat becomes diffusive boundary flux
      xf_Call(xf_AddResTermFaceInt(All, NULL, Face->Group, Face->Number, DS->RS,
                                   xfe_FaceIntFluxHatBC, 1.0));

      if (OutputEval->fidDump != NULL)
        xf_Call(xf_DumpBoundaryIntegralData(All, OutputEval->fidDump, ibfgrp, ibface,
                                            DS->RS->iConv[1] - DS->RS->iConv[0],
                                            DS->RS->iDiff[1] - DS->RS->iDiff[0],
                                            DS->RS->FPI.QuadData->nquad, DS->RS->FPI.QuadData->wquad,
                                            DS->RS->FPI.wn, DS->RS->FPI.xglob,
                                            DS->RS->EqnWS.val[xfe_EWS_F], DS->RS->EqnWS.val[xfe_EWS_qb]));
    } // ibface
  } // ibfg

  // destroy local storage
  xf_Call(xf_DestroyDGStorage(DS));

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemConnectivity_NewDG
int
xf_ElemConnectivity_NewDG(xf_All *All, xf_VectorGroup *UG, xf_Vector *R)
{
  int ierr, s;
  int iiface, rt, rb;
  int FGroup;
  int nIFaceRegular, nIFaceTot;
  real norm;
  xf_FaceInfo FI;
  xf_VectorGroup *RG;
  xf_DGStorage *DS=NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // Locate Residual vector group
  xf_Call(xf_FindSimilarVectorGroup(All, UG, "Residual", xfe_False, xfe_True,
                                    NULL, &RG, NULL));

  // begin halo exchange of vector group UG
  xf_Call(xf_HaloExchangeVectorGroupBegin(UG, -1));

  // determine memory requirements and allocate local storage for linearization
  xf_Call(xf_AllocateDGStorage(All, UG, RG, NULL, NULL, &DS));

  /* residuals on interior faces */

  // regular faces are ones for which L and R elems are on current proc
  // total faces > self faces if some faces are next to halos (in parallel)
  nIFaceRegular = (Mesh->ParallelInfo != NULL) ? Mesh->ParallelInfo->nIFaceRegular : -1;
  nIFaceTot = Mesh->nIFace;
  if (Mesh->ParallelInfo != NULL) nIFaceTot += Mesh->ParallelInfo->nIFaceHalo;

  for (iiface=0; iiface < nIFaceTot; iiface++){
    if (iiface == nIFaceRegular){ // halo exchange must be completed at this point
      xf_Call(xf_HaloExchangeVectorGroupEnd(UG, -1));
    }
    if (!xf_FullLinkageIFace(Mesh, iiface)) continue;

    xf_GetFaceGroup(Mesh, xfe_FaceInterior, -1, xfe_False, &FGroup);
    xf_Call(xf_CalculateResidualFaceDG(All, FGroup, iiface, DS, NULL));

    xf_Call(xf_GetFaceInfo(Mesh, FGroup, iiface, &FI));
    for (s=0; s<2; s++){
      if (s == 0){
        rt = xf_GetVectorRank(DS->U, FI.egrpL, FI.elemL);
        rb = xf_GetVectorRank(DS->U, FI.egrpR, FI.elemR);
      }
      else{
        rt = xf_GetVectorRank(DS->U, FI.egrpR, FI.elemR);
        rb = xf_GetVectorRank(DS->U, FI.egrpL, FI.elemL);
      }
      norm = xf_MatrixFrobeniusNorm(DS->RF_UF[s%2][(s+1)%2].A, rt, rb);
      R->GenArray[FGroup].rValue[iiface][s] = norm;
    } // side
  } // iiface

  // destroy local storage
  xf_Call(xf_DestroyDGStorage(DS));

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PingResidual_NewDG
int
xf_PingResidual_NewDG(xf_All *All, xf_VectorGroup *UG, real eps, real tol)
{
  int ierr, j, face, faceN;
  int egrp, elem, egrpN, elemN;
  int sr;
  int nn, nnN, r, rN;
  int ieps, rUmax;
  enum xfe_Verbosity Verbosity;
  real ep;
  real *ER, *ER0, *EU;
  real **veps = NULL;
  xf_Mesh *Mesh;
  xf_VectorGroup *RG0, *RG;
  xf_Vector *U, *R0, *R;
  xf_JacobianMatrix *R_U0, *R_U;
  xf_SolverData *SolverData;

  // pull off Mesh and state rank
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;

  // locate two residual vectors
  xf_Call(xf_FindSimilarVectorGroup(All, UG, "Residual0", xfe_False, xfe_True, NULL, &RG0, NULL));
  xf_Call(xf_FindSimilarVectorGroup(All, UG, "Residual", xfe_False, xfe_True, NULL, &RG, NULL));

  // initialize residuals to zero
  xf_Call(xf_SetZeroVectorGroup(RG0));
  xf_Call(xf_SetZeroVectorGroup(RG));

  // locate two residual Jacobian matrices, R_U
  xf_Call(xf_FindJacobianMatrix(All, All->DataSet, UG, "R_U0", xfe_True, NULL, NULL, &R_U0, NULL));
  xf_Call(xf_FindJacobianMatrix(All, All->DataSet, UG, "R_U", xfe_True, NULL, NULL, &R_U, NULL));


  // initialize Jacobians to zero
  xf_Call(xf_SetZeroJacobian(R_U0));
  xf_Call(xf_SetZeroJacobian(R_U));


  // create/allocate SolverData
  xf_Call(xf_CreateSolverData(All, &SolverData));

  // determine verbosity level
  xf_Call(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName,
                             (int ) xfe_VerbosityLast, (int *) &Verbosity));

  // calculate baseline residual and Jacobian
  xf_Call(xf_CalculateResidual_NewDG(All, UG, RG0, R_U0, SolverData));

  // get elem state vector, U
  xf_Call(xf_GetR_UState( UG, xfe_SpaceSchemeDG, &U));

  // get residuals, R and R0
  xf_Call(xf_GetR_UState( RG0, xfe_SpaceSchemeDG, &R0));
  xf_Call(xf_GetR_UState( RG , xfe_SpaceSchemeDG, &R ));

  // loop over arrays of U (element groups)
  for (egrp=0; egrp<U->nArray; egrp++){

    // loop over elements
    for (elem=0; elem<U->GenArray[egrp].n; elem++){

      if (Verbosity != xfe_VerbosityLow)
        xf_printf("\n *** egrp = %d, elem = %d *** \n\n", egrp, elem);


      // rank of U vector on current element
      nn = xf_Jacobian_n(R_U, egrp, elem);
      r = nn*sr;

      // pointer to U state on face
      EU = U->GenArray[egrp].rValue[elem];

      // max rank of self + neighbors
      for (face=0, rUmax=r; face<R_U->nFace[egrp][elem]; face++){
        nnN = xf_Jacobian_n(R_U, R_U->egrpN[egrp][elem][face], R_U->elemN[egrp][elem][face]);
        rN = nnN*sr;
        rUmax = max(rUmax, rN);
      }

      // reallocate space for storing ping diffs for checking rates
      xf_Call(xf_ReAlloc2((void ***) &veps, R_U->nFace[egrp][elem]+1, rUmax, sizeof(real)));

      // loop over dofs on current face
      for (j=0; j<r; j++){

        // loop over two epsilons
        for (ieps=0, ep=eps; ieps<2; ieps++, ep*=0.5){

          // perturb state
          EU[j] += ep;

          // recalculate residual and Jacobian
          xf_Call(xf_CalculateResidual_NewDG(All, UG, RG, R_U, SolverData));

          // loop over self and neighbor faces
          for (face=-1; face<R_U->nFace[egrp][elem]; face++){
            //for (face=-1; face<0; face++){

            // set neighbor face info
            if (face == -1){
              egrpN = egrp;
              elemN = elem;
              faceN = face;
            }
            else{
              // continue if no block here
              if (R_U->Value[egrp][elem][face] == NULL) continue;
              egrpN = R_U->egrpN[egrp][elem][face];
              elemN = R_U->elemN[egrp][elem][face];
              faceN = R_U->faceN[egrp][elem][face];
            }

            // skip boundary faces
            if (egrpN < 0) continue;

            // pointer to neighbor residual (orig and perturbed)
            ER0 = R0->GenArray[egrpN].rValue[elemN];
            ER  =  R->GenArray[egrpN].rValue[elemN];

            // rank of U vector on neighbor
            nnN = xf_Jacobian_n(R_U, egrpN, elemN);
            rN = nnN*sr;

            // ping derivatives
            xf_Call(xf_PingVector(ER, ER0, rN, R_U->Value[egrpN][elemN][1+faceN],
                                  R_U0->Value[egrpN][elemN][1+faceN], r, sr, j,
                                  Verbosity, "RU_U", ep, tol, ieps, veps[face+1]));

          } // face

          // unperturb state
          EU[j] -= ep;

        } // ieps

      } // j

    } // ibface

  } // ibfgrp

  // destroy SolverData
  xf_Call(xf_DestroySolverData(SolverData));

  xf_Release2((void **) veps);

  return xf_OK;
}




#if( UNIT_TEST==1 )
#include "xf_ResidualNewDG.test.in"
#endif
