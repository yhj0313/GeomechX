static char help[] = "Transient poroelasticity - Linear elasticity + darcy's law considering fluid storage in 2d and 3d with finite elements.\n\n\n";

#include "TransientPoroElasticity_main.hh"

PetscErrorCode main_transientporoelasticity(MPI_Comm comm) {
    
  DM             dm, dmAux; 
  TS             ts;
  Vec            solutions; /* displacemnet[dim], darcyvelocity[dim], pressure[1] */
  TransientPoroElasticity        user;
  PetscFE        fe[3], fe_derived;
  PetscInt       dim;
  PetscDS        ds;

  const char *name[3] = {"displacement", "pressure", "volumetric_strain"};

  
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.transientdarcysflow.SetInitialPressure(PETSC_COMM_WORLD));
  PetscCall(user.SetMaterial(PETSC_COMM_WORLD, "granite", 0, "domain rock"));
  if (user.isotropiclinearelasticity.cryer == PETSC_TRUE) {PetscCall(Mesh::Create_cryer(PETSC_COMM_WORLD, &dm));}
  else {PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));}
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));

  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 
  PetscCall(user.SetupTS(dm, &solutions, "solution", &ts));
  PetscCall(TSMonitorSet(ts, TransientPoroElasticity::SolutionMonitor_HM, (void (*))&user, NULL));
  PetscCall(user.SolveFE(dm, &solutions, "solution", &ts));

  if (user.transientdarcysflow.expTS) PetscCall(DMTSDestroyRHSMassMatrix(dm));
  PetscCall(VecDestroy(&solutions));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));
  PetscCall(PetscFEDestroy(&fe[1]));
  PetscCall(PetscFEDestroy(&fe[2]));

  return 0;
}