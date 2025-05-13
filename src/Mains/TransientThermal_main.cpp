static char help[] = "Transient state heat conduction and convection in 2d and 3d with finite elements.\n\n\n";

#include "TransientThermal_main.hh"

PetscErrorCode main_transientthermal(MPI_Comm comm)
{
  DM             dm, dmAux;
  TS             ts;
  Vec            solutions;
  TransientThermal        user;
  PetscFE        fe[1]; 
  PetscInt       dim;
  PetscDS        ds;
  const char *name[1] = {"temperature"};

  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.SetInitialTemperature(PETSC_COMM_WORLD));
  PetscCall(user.SetMaterials(PETSC_COMM_WORLD));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));
  
  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 

  PetscCall(user.SetupTS(dm, &solutions, "solution", &ts));
  PetscCall(user.SolveFE(dm, &solutions, "solution", &ts));
  PetscCall(DMGetDimension(dm, &dim));

  if (user.expTS) PetscCall(DMTSDestroyRHSMassMatrix(dm));
  PetscCall(VecDestroy(&solutions));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));

  return 0;
}
