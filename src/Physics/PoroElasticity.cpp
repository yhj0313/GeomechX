#include "PoroElasticity.hh"


PetscErrorCode PoroElasticity::CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace)
{
  PetscFunctionBegin;
  PetscCall(DMPlexCreateRigidBody(dm, origField, nullspace));
  PetscFunctionReturn(0);
}

PetscErrorCode PoroElasticity::CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim) {
  PetscBool      simplex;
  DM             coordDM;
  char           prefix[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;

  DMPlexIsSimplex(dm, &simplex);
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, numComp, simplex, name ? prefix : NULL, -1, fe_derived));
  PetscCall(PetscObjectSetName((PetscObject) *fe_derived, name));
  PetscCall(PetscFECopyQuadrature(*fe, *fe_derived));
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  PetscCall(DMClone(dm, dmAux));
  PetscCall(DMSetCoordinateDM(*dmAux, coordDM));
  PetscCall(DMSetField(*dmAux, 0, NULL, (PetscObject) *fe_derived));
  PetscCall(DMCreateDS(*dmAux));

  PetscFunctionReturn(0);  
}

PetscErrorCode PoroElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]) {
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscInt       n;

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, 0.0, *sol, &func, INSERT_VALUES, *derivedsol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  

  PetscFunctionReturn(0);  
}

PetscErrorCode PoroElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) {
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, time, *sol, &func, INSERT_VALUES, *derivedsol));
  PetscCall(DMSetOutputSequenceNumber(*dmAux, step, time));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  


  PetscFunctionReturn(0);  
}

PetscErrorCode PoroElasticity::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) 
{
//   PetscInt       n;
  
  PetscFunctionBeginUser;
//   SolveDerivedField(dm, &dmAux, Pw_Functions::scale_darcyvelocity, sol, &derivedsols[0], "darcyvelocity");
//   SolveDerivedField(dm, &dmAux, Pw_Functions::scale_pressure, sol, &derivedsols[1], "pressure");

  // switch (this->derivedsoltype) {
  //   case PLANESTRAIN_2D:
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_planestrain, sol, &derivedsols[0], "strain");
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_planestrain, sol, &derivedsols[1], "stress");
  //   break;
  //   case CAUCHY_3D:
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_3D, sol, &derivedsols[0], "strain");
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_3D, sol, &derivedsols[1], "stress");
  //   break;
  //   case AXISYMMETRIC_2D_PLANESTRAIN:
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_strain_planestrain, sol, &derivedsols[0], "strain");
  //   SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_stress_planestrain, sol, &derivedsols[1], "stress");
  //   break;
  //   default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type");

  // }
  
  // // add
  // PetscCall(DMGetNumAuxiliaryVec(dm, &n));
  // std::cout << n << std::endl;

  PetscFunctionReturn(0);

}

PoroElasticity::PoroElasticity() {
  // this->useNearNullspace = PETSC_TRUE;
  // this->radial_tangential_2D = PETSC_FALSE;
  this->material_view = PETSC_FALSE;
}