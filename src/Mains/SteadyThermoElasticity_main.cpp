static char help[] = "Linear elasticity + Steady state heat conduction and convection in 2d and 3d with finite elements.\n\n\n";

#include "SteadyThermoElasticity_main.hh"


PetscErrorCode main_steadythermoelasticity(MPI_Comm comm)
{
  DM             dm, dmAux;
  SNES           snes; /* Nonlinear solver */
  Vec            solutions; /* temperature, displcement */
  Vec            derived_solutions[3]; /* heat flux, stress, strain */
  SteadyThermoElasticity        user;
  PetscFE        fe[2], fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  const char *name[2] = {"displacement", "temperature"};
  
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.SetMaterial(PETSC_COMM_WORLD, "granite", 0, "domain rock"));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(user.SetupFE(dm, ds, name, fe)); 
  PetscCall(user.SolveFE(dm, &solutions, "solution", &snes));
  PetscCall(DMGetDimension(dm, &dim));
  if (dim == 3) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe[0], &fe_derived, 6, dim));}
  else if (dim == 2) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe[0], &fe_derived, 4, dim));}

  user.SolveDerivedFields(dm, dmAux, &solutions, &derived_solutions[0]) ;
  PetscCall(VecDestroy(&solutions));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm));
  PetscCall(DMDestroy(&dmAux));  
  PetscCall(PetscFEDestroy(fe));
  PetscCall(PetscFEDestroy(&fe_derived));
  return 0;
}
