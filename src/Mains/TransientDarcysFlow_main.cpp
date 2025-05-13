static char help[] = "Darcy's flow in a transient state in 2d and 3d with finite elements.\n\n\n";

#include "TransientDarcysFlow_main.hh"

PetscErrorCode main_transientdarcysflow(MPI_Comm comm) {
    
  DM             dm, dmAux; 
  TS             ts;
  Vec            solutions; /* fluid flux[dim], gravity[1] */
  Vec            derived_solutions[2]; /* pressure, darcy velocity */
  TransientDarcysFlow        user;
  PetscFE        fe[1], fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  const char *name[1] = {"pressure"};

  
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.SetInitialPressure(PETSC_COMM_WORLD));
  PetscCall(user.SetMaterial(PETSC_COMM_WORLD, "granite", 0, "domain rock"));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));

  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));

  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 
  PetscCall(user.SolveFE(dm, &solutions, "solution", &ts));

  if (user.expTS) PetscCall(DMTSDestroyRHSMassMatrix(dm));
  PetscCall(VecDestroy(&solutions));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));

  return 0;
}