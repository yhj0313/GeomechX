static char help[] = "Steady state heat conduction and convection in 2d and 3d with finite elements.\n\n\n";

#include "SteadyThermal_main.hh"


PetscErrorCode main_steadythermal(MPI_Comm comm)
{
  DM             dm, dmAux; //, dmbd;
  SNES           snes; /* Nonlinear solver */
  Vec            temperature;
  Vec            derived_solutions[1]; /* heat flux */
  SteadyThermal        user;
  PetscFE        fe, fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.SetMaterials(PETSC_COMM_WORLD));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(user.SetupFE(dm, ds, "steady_heat", &fe)); 
  PetscCall(user.SolveFE(dm, &temperature, "temperature", &snes));
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(user.CreateDerivedField(dm, &dmAux, "thermal_derived", &fe, &fe_derived, dim, dim)); /* heat flux (# of component = dim) */
  user.SolveDerivedFields(dm, dmAux, &temperature, &derived_solutions[0]) ;
  PetscCall(VecDestroy(&temperature));
  PetscCall(VecDestroy(&derived_solutions[0]));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe));
  PetscCall(PetscFEDestroy(&fe_derived));

  return 0;
}
