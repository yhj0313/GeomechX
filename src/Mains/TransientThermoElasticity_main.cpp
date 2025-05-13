static char help[] = "Transient thermoelasticity - Linear elasticity + heat conduction and convection in 2d and 3d with finite elements.\n\n\n";

#include "TransientThermoElasticity_main.hh"


PetscErrorCode main_transientthermoelasticity(MPI_Comm comm)
{
  DM             dm, dmAux;
  TS             ts;
  Vec            solutions; /* displcement, temperature */
  TransientThermoElasticity        user;
  PetscFE        fe[2], fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  const char *name[2] = {"displacement", "temperature"};
  
  PetscFunctionBegin;
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.transientthermal.SetInitialTemperature(PETSC_COMM_WORLD));
  PetscCall(user.SetMaterials(PETSC_COMM_WORLD));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));
  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 
  PetscCall(user.SetupTS(dm, &solutions, "solution", &ts));
  PetscCall(TSMonitorSet(ts, TransientThermoElasticity::SolutionMonitor_TM, (void (*))&user, NULL));
  PetscCall(user.SolveFE(dm, &solutions, "solution", &ts));
  if (user.transientthermal.expTS) PetscCall(DMTSDestroyRHSMassMatrix(dm));
  PetscCall(VecDestroy(&solutions));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));
  PetscCall(PetscFEDestroy(&fe[1]));
   PetscFunctionReturn(0); 
}
