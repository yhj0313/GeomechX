static char help[] = "Transient thermo-poro-elasticity - Linear elasticity + heat conduction and convection + darcy's law considering fluid storage in 2d and 3d with finite elements.\n\n\n";

#include "TransientThermoPoroElasticity_main.hh"


PetscErrorCode main_transientthermoporoelasticity(MPI_Comm comm)
{
  DM             dm, dmAux;
  TS             ts;
  Vec            solutions; /* displcement[dim], temperature[1], darcyvelocity[dim], pressure[1] */
  TransientThermoPoroElasticity   user;
  PetscFE        fe[4], fe_derived;
  PetscInt       dim;
  PetscDS        ds;
  const char *name[4] = {"displacement", "temperature", "pressure", "volumetric_strain"};
  
  PetscFunctionBegin;
  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));
  PetscCall(user.transientthermal.SetInitialTemperature(PETSC_COMM_WORLD));
  PetscCall(user.transientdarcysflow.SetInitialPressure(PETSC_COMM_WORLD));
  
  PetscCall(user.SetMaterials(PETSC_COMM_WORLD));
  PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));
  PetscCall(user.SetupFE(dm, ds, name, &fe[0])); 
  PetscCall(user.SetupTS(dm, &solutions, "solution", &ts));
  PetscCall(TSMonitorSet(ts, TransientThermoPoroElasticity::SolutionMonitor_THM, (void (*))&user, NULL));
  PetscCall(user.SolveFE(dm, &solutions, "solution", &ts));
  if (user.expTS) PetscCall(DMTSDestroyRHSMassMatrix(dm));
  PetscCall(VecDestroy(&solutions));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFEDestroy(&fe[0]));
  PetscCall(PetscFEDestroy(&fe[1]));
  PetscCall(PetscFEDestroy(&fe[2]));
  PetscCall(PetscFEDestroy(&fe[3]));

  PetscFunctionReturn(0); 
}
