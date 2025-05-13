#include "Projection.hh"

// moved from $PETSC_DIR/src/dm/impls/plex/plex.c
static inline void add(PetscScalar *x, PetscScalar y)
{
  *x += y;
}

// moved from $PETSC_DIR/src/dm/impls/plex/plex.c
static inline void insert(PetscScalar *x, PetscScalar y)
{
  *x = y;
}

// moved from $PETSC_DIR/src/dm/impls/plex/plex.c
static inline PetscErrorCode CheckPoint_Private(DMLabel label, PetscInt labelId, PetscSection section, PetscInt point, PetscInt f, PetscInt *offset, PetscBool *contains)
{
  PetscFunctionBegin;
  *contains = PETSC_TRUE;
  if (label) {
    PetscInt fdof;

    PetscCall(DMLabelStratumHasPoint(label, labelId, point, contains));
    if (!*contains) {
      PetscCall(PetscSectionGetFieldDof(section, point, f, &fdof));
      *offset += fdof;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// moved from $PETSC_DIR/src/dm/impls/plex/plex.c
static inline PetscErrorCode updatePointFields_private(PetscSection section, PetscInt point, const PetscInt *perm, const PetscScalar *flip, PetscInt f, void (*fuse)(PetscScalar *, PetscScalar), PetscBool setBC, const PetscInt clperm[], const PetscScalar values[], PetscInt *offset, PetscScalar array[])
{
  PetscScalar    *a;
  PetscInt        fdof, foff, fcdof, foffset = *offset;
  const PetscInt *fcdofs; /* The indices of the constrained dofs for field f on this point */
  PetscInt        cind = 0, b;

  PetscFunctionBegin;
  PetscCall(PetscSectionGetFieldDof(section, point, f, &fdof));
  PetscCall(PetscSectionGetFieldConstraintDof(section, point, f, &fcdof));
  PetscCall(PetscSectionGetFieldOffset(section, point, f, &foff));
  a = &array[foff];
  if (!fcdof || setBC) {
    if (clperm) {
      if (perm) {
        for (b = 0; b < fdof; b++) fuse(&a[b], values[clperm[foffset + perm[b]]] * (flip ? flip[perm[b]] : 1.));
      } else {
        for (b = 0; b < fdof; b++) fuse(&a[b], values[clperm[foffset + b]] * (flip ? flip[b] : 1.));
      }
    } else {
      if (perm) {
        for (b = 0; b < fdof; b++) fuse(&a[b], values[foffset + perm[b]] * (flip ? flip[perm[b]] : 1.));
      } else {
        for (b = 0; b < fdof; b++) fuse(&a[b], values[foffset + b] * (flip ? flip[b] : 1.));
      }
    }
  } else {
    PetscCall(PetscSectionGetFieldConstraintIndices(section, point, f, &fcdofs));
    if (clperm) {
      if (perm) {
        for (b = 0; b < fdof; b++) {
          if ((cind < fcdof) && (b == fcdofs[cind])) {
            ++cind;
            continue;
          }
          fuse(&a[b], values[clperm[foffset + perm[b]]] * (flip ? flip[perm[b]] : 1.));
        }
      } else {
        for (b = 0; b < fdof; b++) {
          if ((cind < fcdof) && (b == fcdofs[cind])) {
            ++cind;
            continue;
          }
          fuse(&a[b], values[clperm[foffset + b]] * (flip ? flip[b] : 1.));
        }
      }
    } else {
      if (perm) {
        for (b = 0; b < fdof; b++) {
          if ((cind < fcdof) && (b == fcdofs[cind])) {
            ++cind;
            continue;
          }
          fuse(&a[b], values[foffset + perm[b]] * (flip ? flip[perm[b]] : 1.));
        }
      } else {
        for (b = 0; b < fdof; b++) {
          if ((cind < fcdof) && (b == fcdofs[cind])) {
            ++cind;
            continue;
          }
          fuse(&a[b], values[foffset + b] * (flip ? flip[b] : 1.));
        }
      }
    }
  }
  *offset += fdof;
  PetscFunctionReturn(PETSC_SUCCESS);
}

// moved from $PETSC_DIR/src/dm/impls/plex/plex.c
static inline PetscErrorCode updatePointFieldsBC_private(PetscSection section, PetscInt point, const PetscInt perm[], const PetscScalar flip[], PetscInt f, PetscInt Ncc, const PetscInt comps[], void (*fuse)(PetscScalar *, PetscScalar), const PetscInt clperm[], const PetscScalar values[], PetscInt *offset, PetscScalar array[])
{
  PetscScalar    *a;
  PetscInt        fdof, foff, fcdof, foffset = *offset;
  const PetscInt *fcdofs; /* The indices of the constrained dofs for field f on this point */
  PetscInt        Nc, cind = 0, ncind = 0, b;
  PetscBool       ncSet, fcSet;

  PetscFunctionBegin;
  PetscCall(PetscSectionGetFieldComponents(section, f, &Nc));
  PetscCall(PetscSectionGetFieldDof(section, point, f, &fdof));
  PetscCall(PetscSectionGetFieldConstraintDof(section, point, f, &fcdof));
  PetscCall(PetscSectionGetFieldOffset(section, point, f, &foff));
  a = &array[foff];
  if (fcdof) {
    /* We just override fcdof and fcdofs with Ncc and comps */
    PetscCall(PetscSectionGetFieldConstraintIndices(section, point, f, &fcdofs));
    if (clperm) {
      if (perm) {
        if (comps) {
          for (b = 0; b < fdof; b++) {
            ncSet = fcSet = PETSC_FALSE;
            if (b % Nc == comps[ncind]) {
              ncind = (ncind + 1) % Ncc;
              ncSet = PETSC_TRUE;
            }
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              ++cind;
              fcSet = PETSC_TRUE;
            }
            if (ncSet && fcSet) fuse(&a[b], values[clperm[foffset + perm[b]]] * (flip ? flip[perm[b]] : 1.));
          }
        } else {
          for (b = 0; b < fdof; b++) {
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              fuse(&a[b], values[clperm[foffset + perm[b]]] * (flip ? flip[perm[b]] : 1.));
              ++cind;
            }
          }
        }
      } else {
        if (comps) {
          for (b = 0; b < fdof; b++) {
            ncSet = fcSet = PETSC_FALSE;
            if (b % Nc == comps[ncind]) {
              ncind = (ncind + 1) % Ncc;
              ncSet = PETSC_TRUE;
            }
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              ++cind;
              fcSet = PETSC_TRUE;
            }
            if (ncSet && fcSet) fuse(&a[b], values[clperm[foffset + b]] * (flip ? flip[b] : 1.));
          }
        } else {
          for (b = 0; b < fdof; b++) {
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              fuse(&a[b], values[clperm[foffset + b]] * (flip ? flip[b] : 1.));
              ++cind;
            }
          }
        }
      }
    } else {
      if (perm) {
        if (comps) {
          for (b = 0; b < fdof; b++) {
            ncSet = fcSet = PETSC_FALSE;
            if (b % Nc == comps[ncind]) {
              ncind = (ncind + 1) % Ncc;
              ncSet = PETSC_TRUE;
            }
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              ++cind;
              fcSet = PETSC_TRUE;
            }
            if (ncSet && fcSet) fuse(&a[b], values[foffset + perm[b]] * (flip ? flip[perm[b]] : 1.));
          }
        } else {
          for (b = 0; b < fdof; b++) {
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              fuse(&a[b], values[foffset + perm[b]] * (flip ? flip[perm[b]] : 1.));
              ++cind;
            }
          }
        }
      } else {
        if (comps) {
          for (b = 0; b < fdof; b++) {
            ncSet = fcSet = PETSC_FALSE;
            if (b % Nc == comps[ncind]) {
              ncind = (ncind + 1) % Ncc;
              ncSet = PETSC_TRUE;
            }
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              ++cind;
              fcSet = PETSC_TRUE;
            }
            if (ncSet && fcSet) fuse(&a[b], values[foffset + b] * (flip ? flip[b] : 1.));
          }
        } else {
          for (b = 0; b < fdof; b++) {
            if ((cind < fcdof) && (b == fcdofs[cind])) {
              fuse(&a[b], values[foffset + b] * (flip ? flip[b] : 1.));
              ++cind;
            }
          }
        }
      }
    }
  }
  *offset += fdof;
  PetscFunctionReturn(PETSC_SUCCESS);
}

// moved from $PETSC_DIR/src/dm/impls/plex/plexfem.c
static PetscErrorCode DMPlexBasisTransformField_Internal(DM dm, DM tdm, Vec tv, PetscInt p, PetscInt f, PetscBool l2g, PetscScalar *a)
{
  PetscSection       ts;
  const PetscScalar *ta, *tva;
  PetscInt           dof;

  PetscFunctionBeginHot;
  PetscCall(DMGetLocalSection(tdm, &ts));
  PetscCall(PetscSectionGetFieldDof(ts, p, f, &dof));
  PetscCall(VecGetArrayRead(tv, &ta));
  PetscCall(DMPlexPointLocalFieldRead(tdm, p, f, ta, &tva));
  if (l2g) {
    switch (dof) {
    case 4:
      DMPlex_Mult2D_Internal(tva, 1, a, a);
      break;
    case 9:
      DMPlex_Mult3D_Internal(tva, 1, a, a);
      break;
    }
  } else {
    switch (dof) {
    case 4:
      DMPlex_MultTranspose2D_Internal(tva, 1, a, a);
      break;
    case 9:
      DMPlex_MultTranspose3D_Internal(tva, 1, a, a);
      break;
    }
  }
  PetscCall(VecRestoreArrayRead(tv, &ta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMProjectFieldLabelLocal' in $PETSC_DIR/src/dm/interface/dm.c
PetscErrorCode DMProjectFieldLabelLocal_User(DM dm, PetscReal time, DMLabel label, PetscInt numIds, const PetscInt ids[], PetscInt Nc, const PetscInt comps[], Vec localU, void (**funcs)(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar f[]), void **ctxs, InsertMode mode, Vec localX)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(localU, VEC_CLASSID, 8);
  PetscValidHeaderSpecific(localX, VEC_CLASSID, 11);
  // dm->ops->projectfieldlabellocal 대신 DMProjectFieldLabelLocal_Plex 씀
  // 변수에 ctxs 추가함
  PetscCall(DMProjectFieldLabelLocal_Plex_User(dm, time, label, numIds, ids, Nc, comps, localU, funcs, ctxs, mode, localX));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMProjectFieldLabelLocal_Plex' in $PETSC_DIR/src/dm/impls/plex/plexproject.c
PetscErrorCode DMProjectFieldLabelLocal_Plex_User(DM dm, PetscReal time, DMLabel label, PetscInt numIds, const PetscInt ids[], PetscInt Ncc, const PetscInt comps[], Vec localU, void (**funcs)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], PetscInt, const PetscScalar[], void *, PetscScalar[]), void **ctxs, InsertMode mode, Vec localX)
{
  PetscFunctionBegin;
  // 아래에 NULL 대산 ctxs 씀
  PetscCall(DMProjectLocal_Generic_Plex_User(dm, time, localU, Ncc, comps, label, numIds, ids, DM_BC_ESSENTIAL_FIELD, (void (**)(void))funcs, ctxs, mode, localX));
  PetscFunctionReturn(PETSC_SUCCESS);    
}

// modified from 'DMProjectLocal_Generic_Plex' in $PETSC_DIR/src/dm/impls/plex/plexproject.c
static PetscErrorCode DMProjectLocal_Generic_Plex_User(DM dm, PetscReal time, Vec localU, PetscInt Ncc, const PetscInt comps[], DMLabel label, PetscInt numIds, const PetscInt ids[], DMBoundaryConditionType type, void (**funcs)(void), void **ctxs, InsertMode mode, Vec localX)
{
  DM               plex, dmIn, plexIn, dmAux = NULL, plexAux = NULL, tdm;
  DMEnclosureType  encIn, encAux;
  PetscDS          ds = NULL, dsIn = NULL, dsAux = NULL;
  Vec              localA = NULL, tv;
  IS               fieldIS;
  PetscSection     section;
  PetscDualSpace  *sp, *cellsp, *spIn, *cellspIn;
  PetscTabulation *T = NULL, *TAux = NULL;
  PetscInt        *Nc;
  PetscInt         dim, dimEmbed, depth, htInc = 0, htIncIn = 0, htIncAux = 0, minHeight, maxHeight, h, regionNum, Nf, NfIn, NfAux = 0, NfTot, f;
  PetscBool       *isFE, hasFE = PETSC_FALSE, hasFV = PETSC_FALSE, isCohesive = PETSC_FALSE, isCohesiveIn = PETSC_FALSE, transform;
  DMField          coordField;
  DMLabel          depthLabel;
  PetscQuadrature  allPoints = NULL;

  PetscFunctionBegin;
  if (localU) PetscCall(VecGetDM(localU, &dmIn));
  else dmIn = dm;
  PetscCall(DMGetAuxiliaryVec(dm, label, numIds ? ids[0] : 0, 0, &localA));
  if (localA) PetscCall(VecGetDM(localA, &dmAux));
  else dmAux = NULL;
  PetscCall(DMConvert(dm, DMPLEX, &plex));
  PetscCall(DMConvert(dmIn, DMPLEX, &plexIn));
  PetscCall(DMGetEnclosureRelation(dmIn, dm, &encIn));
  PetscCall(DMGetEnclosureRelation(dmAux, dm, &encAux));
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMPlexGetVTKCellHeight(plex, &minHeight));
  PetscCall(DMGetBasisTransformDM_Internal(dm, &tdm));
  PetscCall(DMGetBasisTransformVec_Internal(dm, &tv));
  PetscCall(DMHasBasisTransform(dm, &transform));
  /* Auxiliary information can only be used with interpolation of field functions */
  if (dmAux) {
    PetscCall(DMConvert(dmAux, DMPLEX, &plexAux));
    if (type == DM_BC_ESSENTIAL_FIELD || type == DM_BC_ESSENTIAL_BD_FIELD || type == DM_BC_NATURAL_FIELD) PetscCheck(localA, PETSC_COMM_SELF, PETSC_ERR_USER, "Missing localA vector");
  }
  if (localU && localU != localX) PetscCall(DMPlexInsertBoundaryValues(plex, PETSC_TRUE, localU, time, NULL, NULL, NULL));
  PetscCall(DMGetCoordinateField(dm, &coordField));
  /**** No collective calls below this point ****/
  /* Determine height for iteration of all meshes */
  {
    DMPolytopeType ct, ctIn, ctAux;
    PetscInt       minHeightIn, minHeightAux, lStart, pStart, pEnd, p, pStartIn, pStartAux, pEndAux;
    PetscInt       dim = -1, dimIn = -1, dimAux = -1;

    PetscCall(DMPlexGetSimplexOrBoxCells(plex, minHeight, &pStart, &pEnd));
    if (pEnd > pStart) {
      PetscCall(DMGetFirstLabeledPoint(dm, dm, label, numIds, ids, minHeight, &lStart, NULL));
      p = lStart < 0 ? pStart : lStart;
      PetscCall(DMPlexGetCellType(plex, p, &ct));
      dim = DMPolytopeTypeGetDim(ct);
      PetscCall(DMPlexGetVTKCellHeight(plexIn, &minHeightIn));
      PetscCall(DMPlexGetSimplexOrBoxCells(plexIn, minHeightIn, &pStartIn, NULL));
      PetscCall(DMPlexGetCellType(plexIn, pStartIn, &ctIn));
      dimIn = DMPolytopeTypeGetDim(ctIn);
      if (dmAux) {
        PetscCall(DMPlexGetVTKCellHeight(plexAux, &minHeightAux));
        PetscCall(DMPlexGetSimplexOrBoxCells(plexAux, minHeightAux, &pStartAux, &pEndAux));
        if (pStartAux < pEndAux) {
          PetscCall(DMPlexGetCellType(plexAux, pStartAux, &ctAux));
          dimAux = DMPolytopeTypeGetDim(ctAux);
        }
      } else dimAux = dim;
    } else {
      PetscCall(DMDestroy(&plex));
      PetscCall(DMDestroy(&plexIn));
      if (dmAux) PetscCall(DMDestroy(&plexAux));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
    if (dim < 0) {
      DMLabel spmap = NULL, spmapIn = NULL, spmapAux = NULL;

      /* Fall back to determination based on being a submesh */
      PetscCall(DMPlexGetSubpointMap(plex, &spmap));
      PetscCall(DMPlexGetSubpointMap(plexIn, &spmapIn));
      if (plexAux) PetscCall(DMPlexGetSubpointMap(plexAux, &spmapAux));
      dim    = spmap ? 1 : 0;
      dimIn  = spmapIn ? 1 : 0;
      dimAux = spmapAux ? 1 : 0;
    }
    {
      PetscInt dimProj   = PetscMin(PetscMin(dim, dimIn), (dimAux < 0 ? PETSC_MAX_INT : dimAux));
      PetscInt dimAuxEff = dimAux < 0 ? dimProj : dimAux;

      PetscCheck(PetscAbsInt(dimProj - dim) <= 1 && PetscAbsInt(dimProj - dimIn) <= 1 && PetscAbsInt(dimProj - dimAuxEff) <= 1, PETSC_COMM_SELF, PETSC_ERR_SUP, "Do not currently support differences of more than 1 in dimension");
      if (dimProj < dim) minHeight = 1;
      htInc    = dim - dimProj;
      htIncIn  = dimIn - dimProj;
      htIncAux = dimAuxEff - dimProj;
    }
  }
  PetscCall(DMPlexGetDepth(plex, &depth));
  PetscCall(DMPlexGetDepthLabel(plex, &depthLabel));
  PetscCall(DMPlexGetMaxProjectionHeight(plex, &maxHeight));
  maxHeight = PetscMax(maxHeight, minHeight);
  PetscCheck(maxHeight >= 0 && maxHeight <= dim, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Maximum projection height %" PetscInt_FMT " not in [0, %" PetscInt_FMT ")", maxHeight, dim);
  PetscCall(DMGetFirstLabeledPoint(dm, dm, label, numIds, ids, 0, NULL, &ds));
  if (!ds) PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMGetFirstLabeledPoint(dmIn, dm, label, numIds, ids, 0, NULL, &dsIn));
  if (!dsIn) PetscCall(DMGetDS(dmIn, &dsIn));
  PetscCall(PetscDSGetNumFields(ds, &Nf));
  PetscCall(PetscDSGetNumFields(dsIn, &NfIn));
  PetscCall(PetscDSIsCohesive(dsIn, &isCohesiveIn));
  if (isCohesiveIn) --htIncIn; // Should be rearranged
  PetscCall(DMGetNumFields(dm, &NfTot));
  PetscCall(DMFindRegionNum(dm, ds, &regionNum));
  PetscCall(DMGetRegionNumDS(dm, regionNum, NULL, &fieldIS, NULL, NULL));
  PetscCall(PetscDSIsCohesive(ds, &isCohesive));
  PetscCall(DMGetCoordinateDim(dm, &dimEmbed));
  PetscCall(DMGetLocalSection(dm, &section));
  if (dmAux) {
    PetscCall(DMGetDS(dmAux, &dsAux));
    PetscCall(PetscDSGetNumFields(dsAux, &NfAux));
  }
  PetscCall(PetscDSGetComponents(ds, &Nc));
  PetscCall(PetscMalloc3(Nf, &isFE, Nf, &sp, NfIn, &spIn));
  if (maxHeight > 0) PetscCall(PetscMalloc2(Nf, &cellsp, NfIn, &cellspIn));
  else {
    cellsp   = sp;
    cellspIn = spIn;
  }
  /* Get cell dual spaces */
  for (f = 0; f < Nf; ++f) {
    PetscDiscType disctype;
    //원래 PetscCall(PetscDSGetDiscType_Internal(ds, f, &disctype));
    PetscCall(PetscDSGetDiscType_Internal_User(ds, f, &disctype));
    if (disctype == PETSC_DISC_FE) {
      PetscFE fe;

      isFE[f] = PETSC_TRUE;
      hasFE   = PETSC_TRUE;
      PetscCall(PetscDSGetDiscretization(ds, f, (PetscObject *)&fe));
      PetscCall(PetscFEGetDualSpace(fe, &cellsp[f]));
    } else if (disctype == PETSC_DISC_FV) {
      PetscFV fv;

      isFE[f] = PETSC_FALSE;
      hasFV   = PETSC_TRUE;
      PetscCall(PetscDSGetDiscretization(ds, f, (PetscObject *)&fv));
      PetscCall(PetscFVGetDualSpace(fv, &cellsp[f]));
    } else {
      isFE[f]   = PETSC_FALSE;
      cellsp[f] = NULL;
    }
  }
  for (f = 0; f < NfIn; ++f) {
    PetscDiscType disctype;
    //원래 PetscCall(PetscDSGetDiscType_Internal(ds, f, &disctype));
    PetscCall(PetscDSGetDiscType_Internal_User(dsIn, f, &disctype));
    if (disctype == PETSC_DISC_FE) {
      PetscFE fe;

      PetscCall(PetscDSGetDiscretization(dsIn, f, (PetscObject *)&fe));
      PetscCall(PetscFEGetDualSpace(fe, &cellspIn[f]));
    } else if (disctype == PETSC_DISC_FV) {
      PetscFV fv;

      PetscCall(PetscDSGetDiscretization(dsIn, f, (PetscObject *)&fv));
      PetscCall(PetscFVGetDualSpace(fv, &cellspIn[f]));
    } else {
      cellspIn[f] = NULL;
    }
  }
  for (f = 0; f < Nf; ++f) {
    if (!htInc) {
      sp[f] = cellsp[f];
    } else PetscCall(PetscDualSpaceGetHeightSubspace(cellsp[f], htInc, &sp[f]));
  }
  if (type == DM_BC_ESSENTIAL_FIELD || type == DM_BC_ESSENTIAL_BD_FIELD || type == DM_BC_NATURAL_FIELD) {
    PetscFE          fem, subfem;
    PetscDiscType    disctype;
    const PetscReal *points;
    PetscInt         numPoints;

    PetscCheck(maxHeight <= minHeight, PetscObjectComm((PetscObject)dm), PETSC_ERR_SUP, "Field projection not supported for face interpolation");
    //_User로 수정
    // PetscCall(PetscDualSpaceGetAllPointsUnion(Nf, sp, dim - htInc, funcs, &allPoints));
    PetscCall(PetscDualSpaceGetAllPointsUnion_User(Nf, sp, dim - htInc, funcs, &allPoints));
    PetscCall(PetscQuadratureGetData(allPoints, NULL, NULL, &numPoints, &points, NULL));
    PetscCall(PetscMalloc2(NfIn, &T, NfAux, &TAux));
    for (f = 0; f < NfIn; ++f) {
      if (!htIncIn) {
        spIn[f] = cellspIn[f];
      } else PetscCall(PetscDualSpaceGetHeightSubspace(cellspIn[f], htIncIn, &spIn[f]));
    //원래 PetscCall(PetscDSGetDiscType_Internal(ds, f, &disctype));
      PetscCall(PetscDSGetDiscType_Internal_User(dsIn, f, &disctype));
      if (disctype != PETSC_DISC_FE) continue;
      PetscCall(PetscDSGetDiscretization(dsIn, f, (PetscObject *)&fem));
      if (!htIncIn) {
        subfem = fem;
      } else PetscCall(PetscFEGetHeightSubspace(fem, htIncIn, &subfem));
      PetscCall(PetscFECreateTabulation(subfem, 1, numPoints, points, 1, &T[f]));
    }
    for (f = 0; f < NfAux; ++f) {
     //원래 PetscCall(PetscDSGetDiscType_Internal(ds, f, &disctype));       
      PetscCall(PetscDSGetDiscType_Internal_User(dsAux, f, &disctype));
      if (disctype != PETSC_DISC_FE) continue;
      PetscCall(PetscDSGetDiscretization(dsAux, f, (PetscObject *)&fem));
      if (!htIncAux) {
        subfem = fem;
      } else PetscCall(PetscFEGetHeightSubspace(fem, htIncAux, &subfem));
      PetscCall(PetscFECreateTabulation(subfem, 1, numPoints, points, 1, &TAux[f]));
    }
  }
  /* Note: We make no attempt to optimize for height. Higher height things just overwrite the lower height results. */
  for (h = minHeight; h <= maxHeight; h++) {
    PetscInt     hEff     = h - minHeight + htInc;
    PetscInt     hEffIn   = h - minHeight + htIncIn;
    PetscInt     hEffAux  = h - minHeight + htIncAux;
    PetscDS      dsEff    = ds;
    PetscDS      dsEffIn  = dsIn;
    PetscDS      dsEffAux = dsAux;
    PetscScalar *values;
    PetscBool   *fieldActive;
    PetscInt     maxDegree;
    PetscInt     pStart, pEnd, p, lStart, spDim, totDim, numValues;
    IS           heightIS;

    if (h > minHeight) {
      for (f = 0; f < Nf; ++f) PetscCall(PetscDualSpaceGetHeightSubspace(cellsp[f], hEff, &sp[f]));
    }
    PetscCall(DMPlexGetSimplexOrBoxCells(plex, h, &pStart, &pEnd));
    PetscCall(DMGetFirstLabeledPoint(dm, dm, label, numIds, ids, h, &lStart, NULL));
    PetscCall(DMLabelGetStratumIS(depthLabel, depth - h, &heightIS));
    if (pEnd <= pStart) {
      PetscCall(ISDestroy(&heightIS));
      continue;
    }
    /* Compute totDim, the number of dofs in the closure of a point at this height */
    totDim = 0;
    for (f = 0; f < Nf; ++f) {
      PetscBool cohesive;

      if (!sp[f]) continue;
      PetscCall(PetscDSGetCohesive(ds, f, &cohesive));
      PetscCall(PetscDualSpaceGetDimension(sp[f], &spDim));
      totDim += spDim;
      if (isCohesive && !cohesive) totDim += spDim;
    }
    p = lStart < 0 ? pStart : lStart;
    PetscCall(DMPlexVecGetClosure(plex, section, localX, p, &numValues, NULL));
    PetscCheck(numValues == totDim, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "The output section point (%" PetscInt_FMT ") closure size %" PetscInt_FMT " != dual space dimension %" PetscInt_FMT " at height %" PetscInt_FMT " in [%" PetscInt_FMT ", %" PetscInt_FMT "]", p, numValues, totDim, h, minHeight, maxHeight);
    if (!totDim) {
      PetscCall(ISDestroy(&heightIS));
      continue;
    }
    if (htInc) PetscCall(PetscDSGetHeightSubspace(ds, hEff, &dsEff));
    /* Compute totDimIn, the number of dofs in the closure of a point at this height */
    if (localU) {
      PetscInt totDimIn, pIn, numValuesIn;

      totDimIn = 0;
      for (f = 0; f < NfIn; ++f) {
        PetscBool cohesive;

        if (!spIn[f]) continue;
        PetscCall(PetscDSGetCohesive(dsIn, f, &cohesive));
        PetscCall(PetscDualSpaceGetDimension(spIn[f], &spDim));
        totDimIn += spDim;
        if (isCohesiveIn && !cohesive) totDimIn += spDim;
      }
      PetscCall(DMGetEnclosurePoint(dmIn, dm, encIn, lStart < 0 ? pStart : lStart, &pIn));
      PetscCall(DMPlexVecGetClosure(plexIn, NULL, localU, pIn, &numValuesIn, NULL));
      // TODO We could check that pIn is a cohesive cell for this check
      PetscCheck(isCohesiveIn || (numValuesIn == totDimIn), PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "The input section point (%" PetscInt_FMT ") closure size %" PetscInt_FMT " != dual space dimension %" PetscInt_FMT " at height %" PetscInt_FMT, pIn, numValuesIn, totDimIn, htIncIn);
      if (htIncIn) PetscCall(PetscDSGetHeightSubspace(dsIn, hEffIn, &dsEffIn));
    }
    if (htIncAux) PetscCall(PetscDSGetHeightSubspace(dsAux, hEffAux, &dsEffAux));
    /* Loop over points at this height */
    PetscCall(DMGetWorkArray(dm, numValues, MPIU_SCALAR, &values));
    PetscCall(DMGetWorkArray(dm, NfTot, MPI_INT, &fieldActive));
    {
      const PetscInt *fields;

      PetscCall(ISGetIndices(fieldIS, &fields));
      for (f = 0; f < NfTot; ++f) fieldActive[f] = PETSC_FALSE;
      for (f = 0; f < Nf; ++f) fieldActive[fields[f]] = (funcs[f] && sp[f]) ? PETSC_TRUE : PETSC_FALSE;
      PetscCall(ISRestoreIndices(fieldIS, &fields));
    }
    if (label) {
      PetscInt i;

      for (i = 0; i < numIds; ++i) {
        IS              pointIS, isectIS;
        const PetscInt *points;
        PetscInt        n;
        PetscFEGeom    *fegeom = NULL, *chunkgeom = NULL;
        PetscQuadrature quad = NULL;

        PetscCall(DMLabelGetStratumIS(label, ids[i], &pointIS));
        if (!pointIS) continue; /* No points with that id on this process */
        PetscCall(ISIntersect(pointIS, heightIS, &isectIS));
        PetscCall(ISDestroy(&pointIS));
        if (!isectIS) continue;
        PetscCall(ISGetLocalSize(isectIS, &n));
        PetscCall(ISGetIndices(isectIS, &points));
        PetscCall(DMFieldGetDegree(coordField, isectIS, NULL, &maxDegree));
        if (maxDegree <= 1) PetscCall(DMFieldCreateDefaultQuadrature(coordField, isectIS, &quad));
        if (!quad) {
          if (!h && allPoints) {
            quad      = allPoints;
            allPoints = NULL;
          } else {
            // PetscCall(PetscDualSpaceGetAllPointsUnion(Nf, sp, isCohesive ? dim - htInc - 1 : dim - htInc, funcs, &quad));
            // _User로 수정
            PetscCall(PetscDualSpaceGetAllPointsUnion_User(Nf, sp, isCohesive ? dim - htInc - 1 : dim - htInc, funcs, &quad));
          }
        }
        PetscCall(DMFieldCreateFEGeom(coordField, isectIS, quad, (htInc && h == minHeight) ? PETSC_TRUE : PETSC_FALSE, &fegeom));
        for (p = 0; p < n; ++p) {
          const PetscInt point = points[p];

          PetscCall(PetscArrayzero(values, numValues));
          PetscCall(PetscFEGeomGetChunk(fegeom, p, p + 1, &chunkgeom));
          PetscCall(DMPlexSetActivePoint(dm, point));
          PetscCall(DMProjectPoint_Private_User(dm, dsEff, plexIn, encIn, dsEffIn, plexAux, encAux, dsEffAux, chunkgeom, htInc, time, localU, localA, hasFE, hasFV, isFE, sp, point, T, TAux, type, funcs, ctxs, fieldActive, values));
        //   if (transform) PetscCall(DMPlexBasisTransformPoint_Internal(plex, tdm, tv, point, fieldActive, PETSC_FALSE, values));
        // 수정
          if (transform) PetscCall(DMPlexBasisTransformPoint_Internal_User(plex, tdm, tv, point, fieldActive, PETSC_FALSE, values));
          // 수정
        //   PetscCall(DMPlexVecSetFieldClosure_Internal(plex, section, localX, fieldActive, point, Ncc, comps, label, ids[i], values, mode));
          PetscCall(DMPlexVecSetFieldClosure_Internal_User(plex, section, localX, fieldActive, point, Ncc, comps, label, ids[i], values, mode));
        }
        PetscCall(PetscFEGeomRestoreChunk(fegeom, p, p + 1, &chunkgeom));
        PetscCall(PetscFEGeomDestroy(&fegeom));
        PetscCall(PetscQuadratureDestroy(&quad));
        PetscCall(ISRestoreIndices(isectIS, &points));
        PetscCall(ISDestroy(&isectIS));
      }
    } else {
      PetscFEGeom    *fegeom = NULL, *chunkgeom = NULL;
      PetscQuadrature quad = NULL;
      IS              pointIS;

      PetscCall(ISCreateStride(PETSC_COMM_SELF, pEnd - pStart, pStart, 1, &pointIS));
      PetscCall(DMFieldGetDegree(coordField, pointIS, NULL, &maxDegree));
      if (maxDegree <= 1) PetscCall(DMFieldCreateDefaultQuadrature(coordField, pointIS, &quad));
      if (!quad) {
        if (!h && allPoints) {
          quad      = allPoints;
          allPoints = NULL;
        } else {
        //   PetscCall(PetscDualSpaceGetAllPointsUnion(Nf, sp, dim - htInc, funcs, &quad));
        // _User로 수정
          PetscCall(PetscDualSpaceGetAllPointsUnion_User(Nf, sp, dim - htInc, funcs, &quad));
        }
      }
      PetscCall(DMFieldCreateFEGeom(coordField, pointIS, quad, (htInc && h == minHeight) ? PETSC_TRUE : PETSC_FALSE, &fegeom));
      for (p = pStart; p < pEnd; ++p) {
        PetscCall(PetscArrayzero(values, numValues));
        PetscCall(PetscFEGeomGetChunk(fegeom, p - pStart, p - pStart + 1, &chunkgeom));
        PetscCall(DMPlexSetActivePoint(dm, p));
        PetscCall(DMProjectPoint_Private_User(dm, dsEff, plexIn, encIn, dsEffIn, plexAux, encAux, dsEffAux, chunkgeom, htInc, time, localU, localA, hasFE, hasFV, isFE, sp, p, T, TAux, type, funcs, ctxs, fieldActive, values));
        // if (transform) PetscCall(DMPlexBasisTransformPoint_Internal(plex, tdm, tv, p, fieldActive, PETSC_FALSE, values));
        // 수정
        if (transform) PetscCall(DMPlexBasisTransformPoint_Internal_User(plex, tdm, tv, p, fieldActive, PETSC_FALSE, values));
        //원래
        // PetscCall(DMPlexVecSetFieldClosure_Internal(plex, section, localX, fieldActive, p, Ncc, comps, NULL, -1, values, mode));
        PetscCall(DMPlexVecSetFieldClosure_Internal_User(plex, section, localX, fieldActive, p, Ncc, comps, NULL, -1, values, mode));
      }
      PetscCall(PetscFEGeomRestoreChunk(fegeom, p - pStart, pStart - p + 1, &chunkgeom));
      PetscCall(PetscFEGeomDestroy(&fegeom));
      PetscCall(PetscQuadratureDestroy(&quad));
      PetscCall(ISDestroy(&pointIS));
    }
    PetscCall(ISDestroy(&heightIS));
    PetscCall(DMRestoreWorkArray(dm, numValues, MPIU_SCALAR, &values));
    PetscCall(DMRestoreWorkArray(dm, Nf, MPI_INT, &fieldActive));
  }
  /* Cleanup */
  if (type == DM_BC_ESSENTIAL_FIELD || type == DM_BC_ESSENTIAL_BD_FIELD || type == DM_BC_NATURAL_FIELD) {
    for (f = 0; f < NfIn; ++f) PetscCall(PetscTabulationDestroy(&T[f]));
    for (f = 0; f < NfAux; ++f) PetscCall(PetscTabulationDestroy(&TAux[f]));
    PetscCall(PetscFree2(T, TAux));
  }
  PetscCall(PetscQuadratureDestroy(&allPoints));
  PetscCall(PetscFree3(isFE, sp, spIn));
  if (maxHeight > 0) PetscCall(PetscFree2(cellsp, cellspIn));
  PetscCall(DMDestroy(&plex));
  PetscCall(DMDestroy(&plexIn));
  if (dmAux) PetscCall(DMDestroy(&plexAux));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMProjectPoint_Private' in $PETSC_DIR/src/dm/impls/plex/plexproject.c
static PetscErrorCode DMProjectPoint_Private_User(DM dm, PetscDS ds, DM dmIn, DMEnclosureType encIn, PetscDS dsIn, DM dmAux, DMEnclosureType encAux, PetscDS dsAux, PetscFEGeom *fegeom, PetscInt effectiveHeight, PetscReal time, Vec localU, Vec localA, PetscBool hasFE, PetscBool hasFV, PetscBool isFE[], PetscDualSpace sp[], PetscInt p, PetscTabulation *T, PetscTabulation *TAux, DMBoundaryConditionType type, void (**funcs)(void), void **ctxs, PetscBool fieldActive[], PetscScalar values[])
{
  PetscFVCellGeom fvgeom;
  PetscInt        dim, dimEmbed;

  PetscFunctionBeginHot;
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMGetCoordinateDim(dm, &dimEmbed));
  if (hasFV) PetscCall(DMPlexComputeCellGeometryFVM(dm, p, &fvgeom.volume, fvgeom.centroid, NULL));
  switch (type) {
  case DM_BC_ESSENTIAL:
  case DM_BC_NATURAL:
    // PetscCall(DMProjectPoint_Func_Private(dm, ds, dmIn, dsIn, time, fegeom, &fvgeom, isFE, sp, (PetscErrorCode(**)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar *, void *))funcs, ctxs, values));
    break;
  case DM_BC_ESSENTIAL_FIELD:
  case DM_BC_NATURAL_FIELD:
  // 아래 함수에 void * 추가
    PetscCall(DMProjectPoint_Field_Private_User(dm, ds, dmIn, encIn, dsIn, dmAux, encAux, dsAux, time, localU, localA, fegeom, sp, p, T, TAux, (void (**)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], PetscInt, const PetscScalar[], void *, PetscScalar[]))funcs, ctxs, values));
    break;
  case DM_BC_ESSENTIAL_BD_FIELD:
    // PetscCall(DMProjectPoint_BdField_Private(dm, ds, dmIn, dsIn, dmAux, encAux, dsAux, time, localU, localA, fegeom, sp, p, T, TAux, (void (**)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]))funcs, ctxs, values));
    break;
  default:
    SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "Unknown boundary condition type: %d", (int)type);
  }
  PetscFunctionReturn(PETSC_SUCCESS);    
}

// modified from 'DMProjectPoint_Field_Private' in $PETSC_DIR/src/dm/impls/plex/plexproject.c:
static PetscErrorCode DMProjectPoint_Field_Private_User(DM dm, PetscDS ds, DM dmIn, DMEnclosureType encIn, PetscDS dsIn, DM dmAux, DMEnclosureType encAux, PetscDS dsAux, PetscReal time, Vec localU, Vec localA, PetscFEGeom *cgeom, PetscDualSpace sp[], PetscInt p, PetscTabulation *T, PetscTabulation *TAux, void (**funcs)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], PetscInt, const PetscScalar[], void *, PetscScalar[]), void **ctxs, PetscScalar values[])
{
  PetscSection       section, sectionAux = NULL;
  PetscScalar       *u, *u_t = NULL, *u_x, *a = NULL, *a_t = NULL, *a_x = NULL, *bc;
  PetscScalar       *coefficients = NULL, *coefficientsAux = NULL;
  PetscScalar       *coefficients_t = NULL, *coefficientsAux_t = NULL;
  const PetscScalar *constants;
  PetscReal         *x;
  PetscInt          *uOff, *uOff_x, *aOff = NULL, *aOff_x = NULL, *Nc, face[2];
  PetscFEGeom        fegeom;
  const PetscInt     dE = cgeom->dimEmbed, *cone, *ornt;
  PetscInt           numConstants, Nf, NfIn, NfAux = 0, f, spDim, d, v, inp, tp = 0;
  PetscBool          isAffine, isCohesive, isCohesiveIn, transform;
  DMPolytopeType     qct;

  PetscFunctionBeginHot;
  PetscCall(PetscDSGetNumFields(ds, &Nf));
  PetscCall(PetscDSGetComponents(ds, &Nc));
  PetscCall(PetscDSIsCohesive(ds, &isCohesive));
  PetscCall(PetscDSGetNumFields(dsIn, &NfIn));
  PetscCall(PetscDSIsCohesive(dsIn, &isCohesiveIn));
  PetscCall(PetscDSGetComponentOffsets(dsIn, &uOff));
  PetscCall(PetscDSGetComponentDerivativeOffsets(dsIn, &uOff_x));
  PetscCall(PetscDSGetEvaluationArrays(dsIn, &u, &bc /*&u_t*/, &u_x));
  PetscCall(PetscDSGetWorkspace(dsIn, &x, NULL, NULL, NULL, NULL));
  PetscCall(PetscDSGetConstants(dsIn, &numConstants, &constants));
  PetscCall(DMHasBasisTransform(dmIn, &transform));
  PetscCall(DMGetLocalSection(dmIn, &section));
  PetscCall(DMGetEnclosurePoint(dmIn, dm, encIn, p, &inp));
  // Get cohesive cell hanging off face
  if (isCohesiveIn) {
    PetscCall(DMPlexGetCellType(dmIn, inp, &qct));
    if ((qct != DM_POLYTOPE_POINT_PRISM_TENSOR) && (qct != DM_POLYTOPE_SEG_PRISM_TENSOR) && (qct != DM_POLYTOPE_TRI_PRISM_TENSOR) && (qct != DM_POLYTOPE_QUAD_PRISM_TENSOR)) {
      DMPolytopeType  ct;
      const PetscInt *support;
      PetscInt        Ns, s;

      PetscCall(DMPlexGetSupport(dmIn, inp, &support));
      PetscCall(DMPlexGetSupportSize(dmIn, inp, &Ns));
      for (s = 0; s < Ns; ++s) {
        PetscCall(DMPlexGetCellType(dmIn, support[s], &ct));
        if ((ct == DM_POLYTOPE_POINT_PRISM_TENSOR) || (ct == DM_POLYTOPE_SEG_PRISM_TENSOR) || (ct == DM_POLYTOPE_TRI_PRISM_TENSOR) || (ct == DM_POLYTOPE_QUAD_PRISM_TENSOR)) {
          inp = support[s];
          break;
        }
      }
      PetscCheck(s < Ns, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cohesive cell not found from face %" PetscInt_FMT, inp);
      PetscCall(PetscDSGetComponentOffsetsCohesive(dsIn, 2, &uOff));
      PetscCall(DMPlexGetOrientedCone(dmIn, inp, &cone, &ornt));
      face[0] = 0;
      face[1] = 0;
    }
  }
  if (localU) PetscCall(DMPlexVecGetClosure(dmIn, section, localU, inp, NULL, &coefficients));
  if (dmAux) {
    PetscInt subp;

    PetscCall(DMGetEnclosurePoint(dmAux, dm, encAux, p, &subp));
    PetscCall(PetscDSGetNumFields(dsAux, &NfAux));
    PetscCall(DMGetLocalSection(dmAux, &sectionAux));
    PetscCall(PetscDSGetComponentOffsets(dsAux, &aOff));
    PetscCall(PetscDSGetComponentDerivativeOffsets(dsAux, &aOff_x));
    PetscCall(PetscDSGetEvaluationArrays(dsAux, &a, NULL /*&a_t*/, &a_x));
    PetscCall(DMPlexVecGetClosure(dmAux, sectionAux, localA, subp, NULL, &coefficientsAux));
  }
  /* Get values for closure */
  isAffine        = cgeom->isAffine;
  fegeom.dim      = cgeom->dim;
  fegeom.dimEmbed = cgeom->dimEmbed;
  if (isAffine) {
    fegeom.v    = x;
    fegeom.xi   = cgeom->xi;
    fegeom.J    = cgeom->J;
    fegeom.invJ = cgeom->invJ;
    fegeom.detJ = cgeom->detJ;
  }
  for (f = 0, v = 0; f < Nf; ++f) {
    PetscQuadrature  allPoints;
    PetscInt         q, dim, numPoints;
    const PetscReal *points;
    PetscScalar     *pointEval;
    PetscBool        cohesive;
    DM               dm;
    // 추가
    void *const ctx = ctxs ? ctxs[f] : NULL;
    // 추가

    if (!sp[f]) continue;
    PetscCall(PetscDSGetCohesive(ds, f, &cohesive));
    PetscCall(PetscDualSpaceGetDimension(sp[f], &spDim));
    if (!funcs[f]) {
      for (d = 0; d < spDim; d++, v++) values[v] = 0.;
      if (isCohesive && !cohesive) {
        for (d = 0; d < spDim; d++, v++) values[v] = 0.;
      }
      continue;
    }
    PetscCall(PetscDualSpaceGetDM(sp[f], &dm));
    PetscCall(PetscDualSpaceGetAllData(sp[f], &allPoints, NULL));
    PetscCall(PetscQuadratureGetData(allPoints, &dim, NULL, &numPoints, &points, NULL));
    PetscCall(DMGetWorkArray(dm, numPoints * Nc[f], MPIU_SCALAR, &pointEval));
    for (q = 0; q < numPoints; ++q, ++tp) {
      PetscInt qpt[2];

      if (isCohesiveIn) {
        PetscCall(PetscDSPermuteQuadPoint(dsIn, ornt[0], f, q, &qpt[0]));
        PetscCall(PetscDSPermuteQuadPoint(dsIn, DMPolytopeTypeComposeOrientationInv(qct, ornt[1], 0), f, q, &qpt[1]));
      }
      if (isAffine) {
        CoordinatesRefToReal(dE, cgeom->dim, fegeom.xi, cgeom->v, fegeom.J, &points[q * dim], x);
      } else {
        fegeom.v    = &cgeom->v[tp * dE];
        fegeom.J    = &cgeom->J[tp * dE * dE];
        fegeom.invJ = &cgeom->invJ[tp * dE * dE];
        fegeom.detJ = &cgeom->detJ[tp];
      }
      if (coefficients) {
        // 수정
        // if (isCohesiveIn) PetscCall(PetscFEEvaluateFieldJets_Hybrid_Internal(dsIn, NfIn, 0, tp, T, face, qpt, T, &fegeom, coefficients, coefficients_t, u, u_x, u_t));
        // else PetscCall(PetscFEEvaluateFieldJets_Internal(dsIn, NfIn, 0, tp, T, &fegeom, coefficients, coefficients_t, u, u_x, u_t));
        if (isCohesiveIn) PetscCall(PetscFEEvaluateFieldJets_Hybrid_Internal_User(dsIn, NfIn, 0, tp, T, face, qpt, T, &fegeom, coefficients, coefficients_t, u, u_x, u_t));
        else PetscCall(PetscFEEvaluateFieldJets_Internal_User(dsIn, NfIn, 0, tp, T, &fegeom, coefficients, coefficients_t, u, u_x, u_t));
      }
    //   if (dsAux) PetscCall(PetscFEEvaluateFieldJets_Internal(dsAux, NfAux, 0, tp, TAux, &fegeom, coefficientsAux, coefficientsAux_t, a, a_x, a_t));
      if (dsAux) PetscCall(PetscFEEvaluateFieldJets_Internal_User(dsAux, NfAux, 0, tp, TAux, &fegeom, coefficientsAux, coefficientsAux_t, a, a_x, a_t));
    //   if (transform) PetscCall(DMPlexBasisTransformApplyReal_Internal(dmIn, fegeom.v, PETSC_TRUE, dE, fegeom.v, fegeom.v, dm->transformCtx));
      if (transform) PetscCall(DMPlexBasisTransformApplyReal_Internal_User(dmIn, fegeom.v, PETSC_TRUE, dE, fegeom.v, fegeom.v, dm->transformCtx));
      // 아래 함수의 변수로 ctx 추가
      (*funcs[f])(dE, NfIn, NfAux, uOff, uOff_x, u, u_t, u_x, aOff, aOff_x, a, a_t, a_x, time, fegeom.v, numConstants, constants, ctx, &pointEval[Nc[f] * q]);
    }
    PetscCall(PetscDualSpaceApplyAll(sp[f], pointEval, &values[v]));
    PetscCall(DMRestoreWorkArray(dm, numPoints * Nc[f], MPIU_SCALAR, &pointEval));
    v += spDim;
    /* TODO: For now, set both sides equal, but this should use info from other support cell */
    if (isCohesive && !cohesive) {
      for (d = 0; d < spDim; d++, v++) values[v] = values[v - spDim];
    }
  }
  if (localU) PetscCall(DMPlexVecRestoreClosure(dmIn, section, localU, inp, NULL, &coefficients));
  if (dmAux) PetscCall(DMPlexVecRestoreClosure(dmAux, sectionAux, localA, p, NULL, &coefficientsAux));
  if (isCohesiveIn) PetscCall(DMPlexRestoreOrientedCone(dmIn, inp, &cone, &ornt));
  PetscFunctionReturn(PETSC_SUCCESS);    
}

// modified from 'PetscDualSpaceGetAllPointsUnion' in $PETSC_DIR/src/dm/impls/plex/plexproject.c
static PetscErrorCode PetscDualSpaceGetAllPointsUnion_User(PetscInt Nf, PetscDualSpace *sp, PetscInt dim, void (**funcs)(void), PetscQuadrature *allPoints)
{
  PetscReal *points;
  PetscInt   f, numPoints;

  PetscFunctionBegin;
  if (!dim) {
    PetscCall(PetscQuadratureCreate(PETSC_COMM_SELF, allPoints));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  numPoints = 0;
  for (f = 0; f < Nf; ++f) {
    if (funcs[f]) {
      PetscQuadrature fAllPoints;
      PetscInt        fNumPoints;

      PetscCall(PetscDualSpaceGetAllData(sp[f], &fAllPoints, NULL));
      PetscCall(PetscQuadratureGetData(fAllPoints, NULL, NULL, &fNumPoints, NULL, NULL));
      numPoints += fNumPoints;
    }
  }
  PetscCall(PetscMalloc1(dim * numPoints, &points));
  numPoints = 0;
  for (f = 0; f < Nf; ++f) {
    if (funcs[f]) {
      PetscQuadrature  fAllPoints;
      PetscInt         qdim, fNumPoints, q;
      const PetscReal *fPoints;

      PetscCall(PetscDualSpaceGetAllData(sp[f], &fAllPoints, NULL));
      PetscCall(PetscQuadratureGetData(fAllPoints, &qdim, NULL, &fNumPoints, &fPoints, NULL));
      PetscCheck(qdim == dim, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Spatial dimension %" PetscInt_FMT " for dual basis does not match input dimension %" PetscInt_FMT, qdim, dim);
      for (q = 0; q < fNumPoints * dim; ++q) points[numPoints * dim + q] = fPoints[q];
      numPoints += fNumPoints;
    }
  }
  PetscCall(PetscQuadratureCreate(PETSC_COMM_SELF, allPoints));
  PetscCall(PetscQuadratureSetData(*allPoints, dim, 0, numPoints, points, NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'PetscDSGetDiscType_Internal' in $PETSC_DIR/src/dm/dt/interface/dtds.c
PetscErrorCode PetscDSGetDiscType_Internal_User(PetscDS ds, PetscInt f, PetscDiscType *disctype)
{
  PetscObject  obj;
  PetscClassId id;
  PetscInt     Nf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds, PETSCDS_CLASSID, 1);
  PetscAssertPointer(disctype, 3);
  *disctype = PETSC_DISC_NONE;
  PetscCall(PetscDSGetNumFields(ds, &Nf));
  PetscCheck(f < Nf, PetscObjectComm((PetscObject)ds), PETSC_ERR_ARG_SIZ, "Field %" PetscInt_FMT " must be in [0, %" PetscInt_FMT ")", f, Nf);
  PetscCall(PetscDSGetDiscretization(ds, f, &obj));
  if (obj) {
    PetscCall(PetscObjectGetClassId(obj, &id));
    if (id == PETSCFE_CLASSID) *disctype = PETSC_DISC_FE;
    else *disctype = PETSC_DISC_FV;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMPlexVecSetFieldClosure_Internal' in $PETSC_DIR/src/dm/impls/plex/plex.c
PetscErrorCode DMPlexVecSetFieldClosure_Internal_User(DM dm, PetscSection section, Vec v, PetscBool fieldActive[], PetscInt point, PetscInt Ncc, const PetscInt comps[], DMLabel label, PetscInt labelId, const PetscScalar values[], InsertMode mode)
{
  PetscSection    clSection;
  IS              clPoints;
  PetscScalar    *array;
  PetscInt       *points = NULL;
  const PetscInt *clp;
  PetscInt        numFields, numPoints, p;
  PetscInt        offset = 0, f;

  PetscFunctionBeginHot;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  if (!section) PetscCall(DMGetLocalSection(dm, &section));
  PetscValidHeaderSpecific(section, PETSC_SECTION_CLASSID, 2);
  PetscValidHeaderSpecific(v, VEC_CLASSID, 3);
  PetscCall(PetscSectionGetNumFields(section, &numFields));
  /* Get points */
  PetscCall(DMPlexGetCompressedClosure(dm, section, point, 0, &numPoints, &points, &clSection, &clPoints, &clp));
  /* Get array */
  PetscCall(VecGetArray(v, &array));
  /* Get values */
  for (f = 0; f < numFields; ++f) {
    const PetscInt    **perms = NULL;
    const PetscScalar **flips = NULL;
    PetscBool           contains;

    if (!fieldActive[f]) {
      for (p = 0; p < numPoints * 2; p += 2) {
        PetscInt fdof;
        PetscCall(PetscSectionGetFieldDof(section, points[p], f, &fdof));
        offset += fdof;
      }
      continue;
    }
    PetscCall(PetscSectionGetFieldPointSyms(section, f, numPoints, points, &perms, &flips));
    switch (mode) {
    case INSERT_VALUES:
      for (p = 0; p < numPoints; p++) {
        const PetscInt     point = points[2 * p];
        const PetscInt    *perm  = perms ? perms[p] : NULL;
        const PetscScalar *flip  = flips ? flips[p] : NULL;
        PetscCall(CheckPoint_Private(label, labelId, section, point, f, &offset, &contains));
        if (!contains) continue;
        PetscCall(updatePointFields_private(section, point, perm, flip, f, insert, PETSC_FALSE, NULL, values, &offset, array));
      }
      break;
    case INSERT_ALL_VALUES:
      for (p = 0; p < numPoints; p++) {
        const PetscInt     point = points[2 * p];
        const PetscInt    *perm  = perms ? perms[p] : NULL;
        const PetscScalar *flip  = flips ? flips[p] : NULL;
        PetscCall(CheckPoint_Private(label, labelId, section, point, f, &offset, &contains));
        if (!contains) continue;
        PetscCall(updatePointFields_private(section, point, perm, flip, f, insert, PETSC_TRUE, NULL, values, &offset, array));
      }
      break;
    case INSERT_BC_VALUES:
      for (p = 0; p < numPoints; p++) {
        const PetscInt     point = points[2 * p];
        const PetscInt    *perm  = perms ? perms[p] : NULL;
        const PetscScalar *flip  = flips ? flips[p] : NULL;
        PetscCall(CheckPoint_Private(label, labelId, section, point, f, &offset, &contains));
        if (!contains) continue;
        PetscCall(updatePointFieldsBC_private(section, point, perm, flip, f, Ncc, comps, insert, NULL, values, &offset, array));
      }
      break;
    case ADD_VALUES:
      for (p = 0; p < numPoints; p++) {
        const PetscInt     point = points[2 * p];
        const PetscInt    *perm  = perms ? perms[p] : NULL;
        const PetscScalar *flip  = flips ? flips[p] : NULL;
        PetscCall(CheckPoint_Private(label, labelId, section, point, f, &offset, &contains));
        if (!contains) continue;
        PetscCall(updatePointFields_private(section, point, perm, flip, f, add, PETSC_FALSE, NULL, values, &offset, array));
      }
      break;
    case ADD_ALL_VALUES:
      for (p = 0; p < numPoints; p++) {
        const PetscInt     point = points[2 * p];
        const PetscInt    *perm  = perms ? perms[p] : NULL;
        const PetscScalar *flip  = flips ? flips[p] : NULL;
        PetscCall(CheckPoint_Private(label, labelId, section, point, f, &offset, &contains));
        if (!contains) continue;
        PetscCall(updatePointFields_private(section, point, perm, flip, f, add, PETSC_TRUE, NULL, values, &offset, array));
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_OUTOFRANGE, "Invalid insert mode %d", mode);
    }
    PetscCall(PetscSectionRestoreFieldPointSyms(section, f, numPoints, points, &perms, &flips));
  }
  /* Cleanup points */
  PetscCall(DMPlexRestoreCompressedClosure(dm, section, point, &numPoints, &points, &clSection, &clPoints, &clp));
  /* Cleanup array */
  PetscCall(VecRestoreArray(v, &array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMPlexBasisTransformPoint_Internal' in $PETSC_DIR/src/dm/impls/plex/plexfem.c
PetscErrorCode DMPlexBasisTransformPoint_Internal_User(DM dm, DM tdm, Vec tv, PetscInt p, PetscBool fieldActive[], PetscBool l2g, PetscScalar *a)
{
  PetscSection    s;
  PetscSection    clSection;
  IS              clPoints;
  const PetscInt *clp;
  PetscInt       *points = NULL;
  PetscInt        Nf, f, Np, cp, dof, d = 0;

  PetscFunctionBegin;
  PetscCall(DMGetLocalSection(dm, &s));
  PetscCall(PetscSectionGetNumFields(s, &Nf));
  PetscCall(DMPlexGetCompressedClosure(dm, s, p, 0, &Np, &points, &clSection, &clPoints, &clp));
  for (f = 0; f < Nf; ++f) {
    for (cp = 0; cp < Np * 2; cp += 2) {
      PetscCall(PetscSectionGetFieldDof(s, points[cp], f, &dof));
      if (!dof) continue;
      if (fieldActive[f]) PetscCall(DMPlexBasisTransformField_Internal(dm, tdm, tv, points[cp], f, l2g, &a[d]));
      d += dof;
    }
  }
  PetscCall(DMPlexRestoreCompressedClosure(dm, s, p, &Np, &points, &clSection, &clPoints, &clp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'PetscFEEvaluateFieldJets_Internal' in $PETSC_DIR/src/dm/dt/fe/interface/fe.c
PetscErrorCode PetscFEEvaluateFieldJets_Internal_User(PetscDS ds, PetscInt Nf, PetscInt r, PetscInt q, PetscTabulation T[], PetscFEGeom *fegeom, const PetscScalar coefficients[], const PetscScalar coefficients_t[], PetscScalar u[], PetscScalar u_x[], PetscScalar u_t[])
{
  PetscInt dOffset = 0, fOffset = 0, f, g;

  for (f = 0; f < Nf; ++f) {
    PetscCheck(r < T[f]->Nr, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Replica number %" PetscInt_FMT " should be in [0, %" PetscInt_FMT ")", r, T[f]->Nr);
    PetscCheck(q < T[f]->Np, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point number %" PetscInt_FMT " should be in [0, %" PetscInt_FMT ")", q, T[f]->Np);
    PetscFE          fe;
    const PetscInt   k       = ds->jetDegree[f];
    const PetscInt   cdim    = T[f]->cdim;
    const PetscInt   dE      = fegeom->dimEmbed;
    const PetscInt   Nq      = T[f]->Np;
    const PetscInt   Nbf     = T[f]->Nb;
    const PetscInt   Ncf     = T[f]->Nc;
    const PetscReal *Bq      = &T[f]->T[0][(r * Nq + q) * Nbf * Ncf];
    const PetscReal *Dq      = &T[f]->T[1][(r * Nq + q) * Nbf * Ncf * cdim];
    const PetscReal *Hq      = k > 1 ? &T[f]->T[2][(r * Nq + q) * Nbf * Ncf * cdim * cdim] : NULL;
    PetscInt         hOffset = 0, b, c, d;

    PetscCall(PetscDSGetDiscretization(ds, f, (PetscObject *)&fe));
    for (c = 0; c < Ncf; ++c) u[fOffset + c] = 0.0;
    for (d = 0; d < dE * Ncf; ++d) u_x[fOffset * dE + d] = 0.0;
    for (b = 0; b < Nbf; ++b) {
      for (c = 0; c < Ncf; ++c) {
        const PetscInt cidx = b * Ncf + c;

        u[fOffset + c] += Bq[cidx] * coefficients[dOffset + b];
        for (d = 0; d < cdim; ++d) u_x[(fOffset + c) * dE + d] += Dq[cidx * cdim + d] * coefficients[dOffset + b];
      }
    }
    if (k > 1) {
      for (g = 0; g < Nf; ++g) hOffset += T[g]->Nc * dE;
      for (d = 0; d < dE * dE * Ncf; ++d) u_x[hOffset + fOffset * dE * dE + d] = 0.0;
      for (b = 0; b < Nbf; ++b) {
        for (c = 0; c < Ncf; ++c) {
          const PetscInt cidx = b * Ncf + c;

          for (d = 0; d < cdim * cdim; ++d) u_x[hOffset + (fOffset + c) * dE * dE + d] += Hq[cidx * cdim * cdim + d] * coefficients[dOffset + b];
        }
      }
      PetscCall(PetscFEPushforwardHessian(fe, fegeom, 1, &u_x[hOffset + fOffset * dE * dE]));
    }
    PetscCall(PetscFEPushforward(fe, fegeom, 1, &u[fOffset]));
    PetscCall(PetscFEPushforwardGradient(fe, fegeom, 1, &u_x[fOffset * dE]));
    if (u_t) {
      for (c = 0; c < Ncf; ++c) u_t[fOffset + c] = 0.0;
      for (b = 0; b < Nbf; ++b) {
        for (c = 0; c < Ncf; ++c) {
          const PetscInt cidx = b * Ncf + c;

          u_t[fOffset + c] += Bq[cidx] * coefficients_t[dOffset + b];
        }
      }
      PetscCall(PetscFEPushforward(fe, fegeom, 1, &u_t[fOffset]));
    }
    fOffset += Ncf;
    dOffset += Nbf;
  }
  return PETSC_SUCCESS;
}

// modified from 'PetscFEEvaluateFieldJets_Hybrid_Internal' in $PETSC_DIR/src/dm/dt/fe/interface/fe.c
PetscErrorCode PetscFEEvaluateFieldJets_Hybrid_Internal_User(PetscDS ds, PetscInt Nf, PetscInt rc, PetscInt qc, PetscTabulation Tab[], const PetscInt rf[], const PetscInt qf[], PetscTabulation Tabf[], PetscFEGeom *fegeom, const PetscScalar coefficients[], const PetscScalar coefficients_t[], PetscScalar u[], PetscScalar u_x[], PetscScalar u_t[])
{
  PetscInt dOffset = 0, fOffset = 0, f, g;

  /* f is the field number in the DS, g is the field number in u[] */
  for (f = 0, g = 0; f < Nf; ++f) {
    PetscBool isCohesive;
    PetscInt  Ns, s;

    if (!Tab[f]) continue;
    PetscCall(PetscDSGetCohesive(ds, f, &isCohesive));
    Ns = isCohesive ? 1 : 2;
    {
      PetscTabulation T   = isCohesive ? Tab[f] : Tabf[f];
      PetscFE         fe  = (PetscFE)ds->disc[f];
      const PetscInt  dEt = T->cdim;
      const PetscInt  dE  = fegeom->dimEmbed;
      const PetscInt  Nq  = T->Np;
      const PetscInt  Nbf = T->Nb;
      const PetscInt  Ncf = T->Nc;

      for (s = 0; s < Ns; ++s, ++g) {
        const PetscInt   r  = isCohesive ? rc : rf[s];
        const PetscInt   q  = isCohesive ? qc : qf[s];
        const PetscReal *Bq = &T->T[0][(r * Nq + q) * Nbf * Ncf];
        const PetscReal *Dq = &T->T[1][(r * Nq + q) * Nbf * Ncf * dEt];
        PetscInt         b, c, d;

        PetscCheck(r < T->Nr, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Field %" PetscInt_FMT " Side %" PetscInt_FMT " Replica number %" PetscInt_FMT " should be in [0, %" PetscInt_FMT ")", f, s, r, T->Nr);
        PetscCheck(q < T->Np, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Field %" PetscInt_FMT " Side %" PetscInt_FMT " Point number %" PetscInt_FMT " should be in [0, %" PetscInt_FMT ")", f, s, q, T->Np);
        for (c = 0; c < Ncf; ++c) u[fOffset + c] = 0.0;
        for (d = 0; d < dE * Ncf; ++d) u_x[fOffset * dE + d] = 0.0;
        for (b = 0; b < Nbf; ++b) {
          for (c = 0; c < Ncf; ++c) {
            const PetscInt cidx = b * Ncf + c;

            u[fOffset + c] += Bq[cidx] * coefficients[dOffset + b];
            for (d = 0; d < dEt; ++d) u_x[(fOffset + c) * dE + d] += Dq[cidx * dEt + d] * coefficients[dOffset + b];
          }
        }
        PetscCall(PetscFEPushforward(fe, fegeom, 1, &u[fOffset]));
        PetscCall(PetscFEPushforwardGradient(fe, fegeom, 1, &u_x[fOffset * dE]));
        if (u_t) {
          for (c = 0; c < Ncf; ++c) u_t[fOffset + c] = 0.0;
          for (b = 0; b < Nbf; ++b) {
            for (c = 0; c < Ncf; ++c) {
              const PetscInt cidx = b * Ncf + c;

              u_t[fOffset + c] += Bq[cidx] * coefficients_t[dOffset + b];
            }
          }
          PetscCall(PetscFEPushforward(fe, fegeom, 1, &u_t[fOffset]));
        }
        fOffset += Ncf;
        dOffset += Nbf;
      }
    }
  }
  return PETSC_SUCCESS;
}

// modified from 'DMPlexBasisTransformApplyReal_Internal' in $PETSC_DIR/src/dm/impls/plex/plexfem.c
PetscErrorCode DMPlexBasisTransformApplyReal_Internal_User(DM dm, const PetscReal x[], PetscBool l2g, PetscInt dim, const PetscReal *y, PetscReal *z, void *ctx)
{
  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  switch (dim) {
  case 2: {
    PetscScalar yt[2] = {y[0], y[1]}, zt[2] = {0.0, 0.0};

    PetscCall(DMPlexBasisTransformApply_Internal_User(dm, x, l2g, dim, yt, zt, ctx));
    z[0] = PetscRealPart(zt[0]);
    z[1] = PetscRealPart(zt[1]);
  } break;
  case 3: {
    PetscScalar yt[3] = {y[0], y[1], y[2]}, zt[3] = {0.0, 0.0, 0.0};

    PetscCall(DMPlexBasisTransformApply_Internal_User(dm, x, l2g, dim, yt, zt, ctx));
    z[0] = PetscRealPart(zt[0]);
    z[1] = PetscRealPart(zt[1]);
    z[2] = PetscRealPart(zt[2]);
  } break;
  }
#else
  PetscCall(DMPlexBasisTransformApply_Internal_User(dm, x, l2g, dim, y, z, ctx));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

// modified from 'DMPlexBasisTransformApply_Internal' in $PETSC_DIR/src/dm/impls/plex/plexfem.c
PetscErrorCode DMPlexBasisTransformApply_Internal_User(DM dm, const PetscReal x[], PetscBool l2g, PetscInt dim, const PetscScalar *y, PetscScalar *z, void *ctx)
{
  const PetscScalar *A;

  PetscFunctionBeginHot;
  PetscCall((*dm->transformGetMatrix)(dm, x, l2g, &A, ctx));
  switch (dim) {
  case 2:
    DMPlex_Mult2D_Internal(A, 1, y, z);
    break;
  case 3:
    DMPlex_Mult3D_Internal(A, 1, y, z);
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
