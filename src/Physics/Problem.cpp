#include "Problem.hh"

PetscErrorCode SetInitialConditions(TS ts, Vec u)
{
  DM             dm;
  PetscReal      t;

  PetscFunctionBegin;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(TSGetTime(ts, &t));
  PetscCall(DMComputeExactSolution(dm, t, u, NULL));
  PetscFunctionReturn(0);
}

PetscErrorCode CreateSolveDerivedField(DM dm, PetscFE *fe, PetscPointFunc func, Vec *sol, const char name[], PetscReal time, PetscInt step, PetscInt numComp, PetscInt dim) {
  
  PetscBool simplex;
  DM dmAux;
  DM coordDM;
  char prefix[PETSC_MAX_PATH_LEN];
  PetscFE fe_derived;
  Vec derivedsol;
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;

  DMPlexIsSimplex(dm, &simplex);
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject)dm), dim, numComp, simplex, name ? prefix : NULL, -1, &fe_derived));
  PetscCall(PetscObjectSetName((PetscObject)fe_derived, name));
  PetscCall(PetscFECopyQuadrature(*fe, fe_derived));
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  PetscCall(DMClone(dm, &dmAux));
  PetscCall(DMSetCoordinateDM(dmAux, coordDM));
  PetscCall(DMSetField(dmAux, 0, NULL, (PetscObject)fe_derived));
  PetscCall(DMCreateDS(dmAux));
  
  PetscCall(DMCreateGlobalVector(dmAux, &derivedsol));
  PetscCall(VecSet(derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) derivedsol, name));

  PetscCall(DMCopyAuxiliaryVec(dm, dmAux));
  PetscCall(DMProjectField(dmAux, time, *sol, &func, INSERT_VALUES, derivedsol));


  PetscCall(DMSetOutputSequenceNumber(dmAux, step, time));

  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(derivedsol, NULL, view_call));  

  PetscFunctionReturn(0);  
}