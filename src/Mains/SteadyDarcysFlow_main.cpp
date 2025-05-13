static char help[] = "Darcy's flow in a steady state in 2d and 3d with finite elements.\n\n\n";

#include "SteadyDarcysFlow_main.hh"

PetscErrorCode main_steadydarcysflow(MPI_Comm comm) {
    
  DM             dm, dmAux;
  SNES           snes; /* Nonlinear solver */
  Vec            solutions; /* fluid flux[dim], gravity[1] */
  Vec            derived_solutions[2]; /* pressure, darcy velocity */
  SteadyDarcysFlow        user;
  PetscFE        fe[2], fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  const char *name[2] = {"darcyvelocity", "pressure"};

  
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.SetMaterial(PETSC_COMM_WORLD, "granite", 0, "domain rock"));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 
  PetscCall(user.SolveFE(dm, &solutions, "solution", &snes));
  
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(user.CreateDerivedField(dm, &dmAux, "darcyvelocity", &fe[0], &fe_derived, dim, dim)); /* darcy velocity(# of component = dim) */
  PetscCall(user.SolveDerivedField(dm, &dmAux, Pw_Functions::scale_darcyvelocity, &solutions, &derived_solutions[0], "scaled_darcyvelocity"));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));
  PetscCall(user.CreateDerivedField(dm, &dmAux, "pressure", &fe[1], &fe_derived, 1, dim)); /* pressure(# of component = 1) */
  PetscCall(user.SolveDerivedField(dm, &dmAux, Pw_Functions::scale_pressure, &solutions, &derived_solutions[1], "scaled_pressure"));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));

  PetscCall(VecDestroy(&solutions));
  PetscCall(VecDestroy(&derived_solutions[0]));
  PetscCall(VecDestroy(&derived_solutions[1]));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));
  PetscCall(PetscFEDestroy(&fe[1]));

  return 0;
}