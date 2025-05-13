#ifndef PROJECTION_HH
#define PROJECTION_HH

#include "petscheaders.hh"
#include <petsc/private/dmpleximpl.h>
#include <petsc/private/petscfeimpl.h>
#include <petscdmplex.h>
#include <petscdmfield.h>
#include <petsc/private/petscdsimpl.h>

// Functions modified from PETSc source codes located in
// $PETSC_DIR/src/dm/interface/dm.c
// $PETSC_DIR/src/dm/impls/plex/plexproject.c
// $PETSC_DIR/src/dm/dt/interface/dtds.c
// $PETSC_DIR/src/dm/impls/plex/plex.c
// $PETSC_DIR/src/dm/impls/plex/plexfem.c
// $PETSC_DIR/src/dm/dt/fe/interface/fe.c

PetscErrorCode DMProjectFieldLabelLocal_User(DM dm, PetscReal time, DMLabel label, PetscInt numIds, const PetscInt ids[], PetscInt Nc, const PetscInt comps[], Vec localU, void (**funcs)(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar f[]), void **ctxs, InsertMode mode, Vec localX);
PetscErrorCode DMProjectFieldLabelLocal_Plex_User(DM dm, PetscReal time, DMLabel label, PetscInt numIds, const PetscInt ids[], PetscInt Ncc, const PetscInt comps[], Vec localU, void (**funcs)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], PetscInt, const PetscScalar[], void *, PetscScalar[]), void **ctxs, InsertMode mode, Vec localX);
static PetscErrorCode DMProjectLocal_Generic_Plex_User(DM dm, PetscReal time, Vec localU, PetscInt Ncc, const PetscInt comps[], DMLabel label, PetscInt numIds, const PetscInt ids[], DMBoundaryConditionType type, void (**funcs)(void), void **ctxs, InsertMode mode, Vec localX);
static PetscErrorCode DMProjectPoint_Private_User(DM dm, PetscDS ds, DM dmIn, DMEnclosureType encIn, PetscDS dsIn, DM dmAux, DMEnclosureType encAux, PetscDS dsAux, PetscFEGeom *fegeom, PetscInt effectiveHeight, PetscReal time, Vec localU, Vec localA, PetscBool hasFE, PetscBool hasFV, PetscBool isFE[], PetscDualSpace sp[], PetscInt p, PetscTabulation *T, PetscTabulation *TAux, DMBoundaryConditionType type, void (**funcs)(void), void **ctxs, PetscBool fieldActive[], PetscScalar values[]);
static PetscErrorCode DMProjectPoint_Field_Private_User(DM dm, PetscDS ds, DM dmIn, DMEnclosureType encIn, PetscDS dsIn, DM dmAux, DMEnclosureType encAux, PetscDS dsAux, PetscReal time, Vec localU, Vec localA, PetscFEGeom *cgeom, PetscDualSpace sp[], PetscInt p, PetscTabulation *T, PetscTabulation *TAux, void (**funcs)(PetscInt, PetscInt, PetscInt, const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[], PetscReal, const PetscReal[], PetscInt, const PetscScalar[], void *, PetscScalar[]), void **ctxs, PetscScalar values[]);
static PetscErrorCode PetscDualSpaceGetAllPointsUnion_User(PetscInt Nf, PetscDualSpace *sp, PetscInt dim, void (**funcs)(void), PetscQuadrature *allPoints);
PetscErrorCode PetscDSGetDiscType_Internal_User(PetscDS ds, PetscInt f, PetscDiscType *disctype);
PetscErrorCode DMPlexVecSetFieldClosure_Internal_User(DM dm, PetscSection section, Vec v, PetscBool fieldActive[], PetscInt point, PetscInt Ncc, const PetscInt comps[], DMLabel label, PetscInt labelId, const PetscScalar values[], InsertMode mode);
PetscErrorCode DMPlexBasisTransformPoint_Internal_User(DM dm, DM tdm, Vec tv, PetscInt p, PetscBool fieldActive[], PetscBool l2g, PetscScalar *a);
PetscErrorCode PetscFEEvaluateFieldJets_Internal_User(PetscDS ds, PetscInt Nf, PetscInt r, PetscInt q, PetscTabulation T[], PetscFEGeom *fegeom, const PetscScalar coefficients[], const PetscScalar coefficients_t[], PetscScalar u[], PetscScalar u_x[], PetscScalar u_t[]);
PetscErrorCode PetscFEEvaluateFieldJets_Hybrid_Internal_User(PetscDS ds, PetscInt Nf, PetscInt rc, PetscInt qc, PetscTabulation Tab[], const PetscInt rf[], const PetscInt qf[], PetscTabulation Tabf[], PetscFEGeom *fegeom, const PetscScalar coefficients[], const PetscScalar coefficients_t[], PetscScalar u[], PetscScalar u_x[], PetscScalar u_t[]);
PetscErrorCode DMPlexBasisTransformApplyReal_Internal_User(DM dm, const PetscReal x[], PetscBool l2g, PetscInt dim, const PetscReal *y, PetscReal *z, void *ctx);
PetscErrorCode DMPlexBasisTransformApply_Internal_User(DM dm, const PetscReal x[], PetscBool l2g, PetscInt dim, const PetscScalar *y, PetscScalar *z, void *ctx);

#endif