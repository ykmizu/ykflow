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
  FILE:  xf_Residual.c

  This file contains high-level functions for evaluating residuals and
  Jacobians.  Fork between different spatial discretizations happens
  here.

*/

#include "xf_String.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Data.h"
#include "xf_MeshDistance.h"
#include "xf_ResidualDG.h"
#include "xf_ResidualHDG.h"
#include "xf_ResidualNewDG.h"
#include "xf_ResidualNewHDG.h"
#include "xf_Penalty.h"


/******************************************************************/
//   FUNCTION Definition: xf_GetSpaceScheme
enum xfe_SpaceSchemeType
xf_GetSpaceScheme(xf_KeyValue *KeyValue){
  int ierr;
  enum xfe_SpaceSchemeType SpaceScheme;

  // Determine spatial discretization scheme
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "SpaceScheme",
                                     xfe_SpaceSchemeName, (int ) xfe_SpaceSchemeLast,
                                     (int *) &SpaceScheme));

  if (ierr != xf_OK) SpaceScheme = xfe_SpaceSchemeLast;

  return SpaceScheme;
}



/******************************************************************/
//   FUNCTION Definition: xf_ProcessResTerms
int
xf_ProcessResTerms(xf_ResTerms *ResTerms, int iConv[], int iDiff[], int iSource[])
{
/*
PURPOSE:

  Arranges ResTerms so that they are in sequence: first Conv, then
  Diff, then Source.  Returns index ranges of all terms; i.e.

    ResTerms->ResTerm[iConv[0] .. iConv[1]-1] are convection terms
    ResTerms->ResTerm[iDiff[0] .. iDiff[1]-1] are diffusion terms
    ...

  The number of (e.g.) convection terms is iConv[1]-iConv[0].  If, for
  example, no convection terms are present, iConv[1] == iConv[0].
  Only Active residual terms are considered.

INPUTS:

  ResTerms: set of residual terms

OUTPUTS:

  iConv  : [start,end+1] range for convection terms (preallocated 2-int vector)
  iDiff  : [start,end+1] range for diffusion terms (preallocated 2-int vector)
  iSource: [start,end+1] range for source terms (preallocated 2-int vector)

RETURNS:

  Error code

*/
  int ierr, nResTerm, pos, i;
  xf_ResTerm *ResTermNew;

  nResTerm = ResTerms->nResTerm;

  if (nResTerm == 0){
    iConv[0] = iConv[1] = 0;
    iDiff[0] = iDiff[1] = 0;
    iSource[0] = iSource[1] = 0;
    return xf_OK;
  }

  ierr = xf_Error(xf_Alloc((void **) &ResTermNew, nResTerm, sizeof(xf_ResTerm)));
  if (ierr != xf_OK) return ierr;

  pos = 0;

  iConv[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermConv) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iConv[1] = pos;

  iDiff[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermDiff) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iDiff[1] = pos;

  iSource[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermSource) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iSource[1] = pos;

  for (i=0; i<nResTerm; i++)
    if (!ResTerms->ResTerm[i].Active){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }

  /* Replace original pointer vec with new pointer vec */
  xf_Release( (void *) ResTerms->ResTerm);
  ResTerms->ResTerm = ResTermNew;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateSupportedVector
static int
xf_CreateSupportedVector(xf_All *All, const char Title[])
{
/*
PURPOSE:

  Creates a supported vector with name Title.  An example is the
  "WallDistance" vector.

INPUTS:

  All   : All structure
  Title : Name of vector to create

OUTPUTS:

  None  : appropriate vector is stored in All->DataSet if created.

RETURN:

  Error Code
*/
  int ierr;

  if (strncmp(Title, "WallDistance", 12) == 0){
    // create a wall distance vector
    ierr = xf_Error(xf_CalculateDistFcn(All, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  else{
    xf_printf("Auxiliary vector = %s is not supported.\n", Title);
    return xf_Error(xf_NOT_SUPPORTED);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindSupportedVector
int
xf_FindSupportedVector(xf_All *All, const char Title[], xf_Vector **pV)
{
  int ierr;
  xf_Data *D;

  ierr = xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND){
    // try creating vector if supported (e.g. wall distance)
    ierr = xf_Error(xf_CreateSupportedVector(All, Title));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D));
  }
  if (ierr != xf_OK) return xf_Error(ierr);

  (*pV) = (xf_Vector *) D->Data;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemConnectivity
int
xf_ElemConnectivity(xf_All *All, xf_VectorGroup *UG, xf_Vector *R)
{
  int ierr;

  // Call appropriate function
  switch (xf_GetSpaceScheme(All->Param->KeyValue)){
    case xfe_SpaceSchemeDG:
    case xfe_SpaceSchemeNewDG:
      xf_Call(xf_ElemConnectivity_DG(All, UG, R));
      break;
    case xfe_SpaceSchemeHDG:
    case xfe_SpaceSchemeNewHDG:
      xf_Call(xf_ElemConnectivity_NewHDG(All, UG, R));
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RetrieveFcnParams
int
xf_RetrieveFcnParams(xf_All *All, xf_EqnSet *EqnSet,
		     int **pIParam, real **pRParam, int *pnAuxU,
		     xf_Vector ***pAuxU)
{
  int ierr, i, iAux;
  char *Name;
  xf_KeyValue *KeyValue;
  xf_Data *D;

  KeyValue = ((All == NULL) ? NULL : All->Param->KeyValue);

  // Re-allocate memory for IParam and RParam
  ierr = xf_Error(xf_Alloc((void **) pIParam, EqnSet->nIParam, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) pRParam, EqnSet->nRParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<EqnSet->nIParam; i++){
    ierr = xf_GetKeyValueInt(EqnSet->KeyValue, EqnSet->IParamKey[i], (*pIParam)+i);
    if (ierr == xf_STRING_ERROR) // try a boolean read
      ierr = xf_GetKeyValueBool(EqnSet->KeyValue, EqnSet->IParamKey[i],
				(enum xfe_Bool *) (*pIParam)+i);
    if ((ierr == xf_NOT_FOUND) && (KeyValue != NULL)){
      ierr = xf_GetKeyValueInt(KeyValue, EqnSet->IParamKey[i], (*pIParam)+i);
      if (ierr == xf_STRING_ERROR) // try a boolean read
	ierr = xf_Error(xf_GetKeyValueBool(KeyValue, EqnSet->IParamKey[i],
					   (enum xfe_Bool *) (*pIParam)+i));
      if (ierr != xf_OK){
	xf_printf("Error, could not find desired key = %s.\n", EqnSet->IParamKey[i]);
	return xf_Error(ierr);
      }
    }
    else
      if (ierr != xf_OK) return ierr;
  }

  for (i=0; i<EqnSet->nRParam; i++){
    ierr = xf_GetKeyValueReal(EqnSet->KeyValue, EqnSet->RParamKey[i], (*pRParam)+i);
    if ((ierr == xf_NOT_FOUND) && (KeyValue != NULL)){
      ierr = xf_Error(xf_GetKeyValueReal(KeyValue, EqnSet->RParamKey[i], (*pRParam)+i));
      if (ierr != xf_OK){
	xf_printf("Error, could not find desired key = %s.\n", EqnSet->RParamKey[i]);
	return xf_Error(ierr);
      }
    }
    else
      if (ierr != xf_OK) return ierr;
  }

  // Auxiliary vectors to be passed into eqnset functions
  if (pAuxU != NULL) (*pAuxU) = NULL;
  if ((pnAuxU != NULL) &&  (((*pnAuxU) = EqnSet->nAuxU) != 0)){

    if ((pAuxU == NULL) || (All == NULL)) return xf_Error(xf_INPUT_ERROR);

    ierr = xf_Error(xf_Alloc((void **) pAuxU, (*pnAuxU), sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;

    for (iAux=0; iAux<(*pnAuxU); iAux++){
      (*pAuxU)[iAux] = NULL;
      Name = EqnSet->AuxUNames[iAux];

      ierr = xf_FindDataByTitle(All->DataSet, Name, xfe_Vector, &D);
      if (ierr == xf_NOT_FOUND){
	// try creating vector if supported (e.g. wall distance)
	ierr = xf_Error(xf_CreateSupportedVector(All, Name));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_FindDataByTitle(All->DataSet, Name, xfe_Vector, &D));
      }
      if (ierr != xf_OK) return xf_Error(ierr);

      (*pAuxU)[iAux] = (xf_Vector *) D->Data;

    } // iAux
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SortEqnSetBCs
int
xf_SortEqnSetBCs(const xf_Mesh *Mesh, xf_BCs *BCs)
{
  int ierr, nBC, nbfgrp, iBC, ibfgrp;
  enum xfe_Bool found;
  xf_BC *BCnew;

  nBC = BCs->nBC;
  nbfgrp = Mesh->nBFaceGroup;

  if (nBC != nbfgrp) return xf_Error(xf_BOUNDARY_CONDITION_ERROR);

  ierr = xf_Error(xf_Alloc((void **) &BCnew, nBC, sizeof(xf_BC)));
  if (ierr != xf_OK) return ierr;

  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    found = xfe_False;
    for (iBC=0; iBC<nBC; iBC++)
      if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, BCs->BC[iBC].BFGTitle) == 0){
        BCnew[ibfgrp] = BCs->BC[iBC];
        found = xfe_True;
        break;
      }

    if (!found){
      xf_printf("Error. Was not able to find Mesh bfg %s in the EqnSet BC list.\n",
		Mesh->BFaceGroup[ibfgrp].Title);
      return xf_Error(xf_BOUNDARY_CONDITION_ERROR);
    }
  }

  xf_Release( (void *) BCs->BC);
  BCs->BC = BCnew;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_BoundaryFluxIntegral
int
xf_BoundaryFluxIntegral(xf_All *All, xf_VectorGroup *UG,
                        xf_OutputEvalData *OutputEval,
                        int nBFG, char **BFGTitles)
{
  int ierr;
  int i, k, nset;
  int *BFGs = NULL;
  enum xfe_SpaceSchemeType SpaceScheme;

  // determine which boundary face groups to integrate over
  ierr = xf_Error(xf_Alloc((void **) &BFGs, nBFG, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  nset = 0;
  for (i=0; i<nBFG; i++){
    for (k=0; k<All->Mesh->nBFaceGroup; k++)
      if (strcmp(BFGTitles[i], All->Mesh->BFaceGroup[k].Title) == 0){
        BFGs[i] = k;
        nset++;
      }
  }
  if ((nset != nBFG) || (nBFG > All->Mesh->nBFaceGroup))
    return xf_Error(xf_OUT_OF_BOUNDS);

  // Determine spatial discretization scheme
  xf_Call(xf_GetKeyValueEnum(All->Param->KeyValue, "SpaceScheme", xfe_SpaceSchemeName,
                             (int) xfe_SpaceSchemeLast, (int *) &SpaceScheme));

  // Call appropriate function
  switch (SpaceScheme) {
  case xfe_SpaceSchemeDG:
    if (UG->nVector != 1) return xf_Error(xf_NOT_SUPPORTED);
    xf_Call(xf_BoundaryFluxIntegral_DG(All, UG->Vector[0], OutputEval, nBFG, BFGs));
    break;
  case xfe_SpaceSchemeNewDG:
    if (UG->nVector != 1) return xf_Error(xf_NOT_SUPPORTED);
    xf_Call(xf_BoundaryFluxIntegral_NewDG(All, UG, OutputEval, nBFG, BFGs));
    break;
  case xfe_SpaceSchemeHDG:
    xf_Call(xf_BoundaryFluxIntegral_HDG(All, UG, OutputEval, nBFG, BFGs));
    break;
  case xfe_SpaceSchemeNewHDG:
    xf_Call(xf_BoundaryFluxIntegral_NewHDG(All, UG, OutputEval, nBFG, BFGs));
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  xf_Release( (void *) BFGs);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidual
int
xf_CalculateResidual(xf_All *All, xf_VectorGroup *UG, xf_VectorGroup *RG,
		     xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  int ierr;
  enum xfe_SpaceSchemeType SpaceScheme;
  clock_t clock_Start, clock_End;

  if (xf_TIMING) clock_Start = clock();

  // Determine spatial discretization scheme
  xf_Call(xf_GetKeyValueEnum(All->Param->KeyValue, "SpaceScheme", xfe_SpaceSchemeName,
                             (int) xfe_SpaceSchemeLast, (int *) &SpaceScheme));

  // Call appropriate function
  switch (SpaceScheme) {
  case xfe_SpaceSchemeDG:
    xf_Call(xf_CalculateResidual_DG(All, UG, RG, R_U, SolverData));
    break;
  case xfe_SpaceSchemeNewDG:
    xf_Call(yk_CalculateResidual_NewDG(All, UG, RG, R_U, SolverData));
    break;
  case xfe_SpaceSchemeHDG:
    xf_Call(xf_CalculateResidual_HDG(All, UG, RG, R_U, SolverData));
    break;
  case xfe_SpaceSchemeNewHDG:
    xf_Call(xf_CalculateResidual_NewHDG(All, UG, RG, R_U, SolverData));
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  if (xf_TIMING){
    clock_End = clock();
    xf_printf("Residual time = %.10E\n", ((real)(clock_End - clock_Start)) / CLOCKS_PER_SEC);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PingResidual
int
xf_PingResidual(xf_All *All, xf_VectorGroup *UG, real eps, real tol)
{
  int ierr;
  enum xfe_SpaceSchemeType SpaceScheme;

  // Determine spatial discretization scheme
  xf_Call(xf_GetKeyValueEnum(All->Param->KeyValue, "SpaceScheme", xfe_SpaceSchemeName,
                             (int) xfe_SpaceSchemeLast, (int *) &SpaceScheme));

  // Call appropriate function
  switch (SpaceScheme) {
  case xfe_SpaceSchemeDG:
    xf_Call(xf_PingResidual_DG(All, UG, eps, tol));
    break;
  case xfe_SpaceSchemeNewDG:
    xf_Call(xf_PingResidual_NewDG(All, UG, eps, tol));
    break;
  case xfe_SpaceSchemeHDG:
    xf_Call(xf_PingResidual_HDG(All, UG, eps, tol));
    break;
  case xfe_SpaceSchemeNewHDG:
    xf_Call(xf_PingResidual_NewHDG(All, UG, eps, tol));
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_Residual.test.in"
#endif
