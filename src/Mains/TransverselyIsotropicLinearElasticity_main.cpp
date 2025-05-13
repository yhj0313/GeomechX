static char help[] = "Linear Transeversely isotropic elasticity in 2d and 3d with finite elements.\n\n\n";

#include "TransverselyIsotropicLinearElasticity_main.hh"

PetscErrorCode main_transverselyisotropiclinearelasticity(MPI_Comm comm)
{
  DM             dm, dmAux;   /* Problem specification */
  SNES           snes; /* Nonlinear solver */
  Vec            displacement;    /* Solutions */
  Vec            derived_solutions[2]; /* Stress and strain */
  TransverselyIsotropicLinearElasticity           user; /* User-defined work context */
  PetscFE        fe, fe_derived;
  PetscInt       dim; 
  PetscDS        ds;
  PetscInt       n;

  PetscCall(user.ProcessOptions(PETSC_COMM_WORLD));
  PetscCall(user.InputParameterOptions(PETSC_COMM_WORLD, "granite"));

  user.SetMaterial(PETSC_COMM_WORLD, "granite", 0, "domain rock");
  if (user.TIStripLineLoad == PETSC_TRUE) {PetscCall(Mesh::Create_striplineload(PETSC_COMM_WORLD, &dm));}
  else {PetscCall(Mesh::Create(PETSC_COMM_WORLD, &dm));}
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(user.SetupFE(dm, ds, "elasticity", &fe)); 
  PetscCall(user.SolveFE(dm, &displacement, "displacement", &snes));
  PetscCall(DMGetDimension(dm, &dim));
  if (dim == 3)
  {
    if (user.GetDerivedSolType() == 3 ) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe, &fe_derived, 4, dim));}  // if (user.GetDerivedSolType() == AXISYMMETRIC_3D )
    else {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe, &fe_derived, 6, dim));}
  }
  else if (dim == 2) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe, &fe_derived, 4, dim));}
  user.SolveDerivedFields(dm, dmAux, &displacement, &derived_solutions[0]) ;

  // add
  PetscCall(DMGetNumAuxiliaryVec(dmAux, &n));

  /* Cleanup */
  PetscCall(VecDestroy(&displacement));
  PetscCall(VecDestroy(&derived_solutions[0]));
  PetscCall(VecDestroy(&derived_solutions[1]));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm));
  PetscCall(DMDestroy(&dmAux));
  PetscCall(PetscFEDestroy(&fe));
  PetscCall(PetscFEDestroy(&fe_derived));

  return 0;
}
