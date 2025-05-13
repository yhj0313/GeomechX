static char help[] = "Linear isotropic elasticity in 2d and 3d with finite elements.\n\n\n";

#include "IsotropicLinearElasticity_main.hh"


PetscErrorCode main_isotropiclinearelasticity(MPI_Comm comm)
{
  DM             dm, dmAux;   /* Problem specification */
  SNES           snes; /* Nonlinear solver */
  Vec            displacement;    /* Solutions */
  Vec            derived_solutions[2]; /* Stress and strain */
  IsotropicLinearElasticity           user; /* User-defined work context */
  PetscFE        fe, fe_derived;
  PetscInt       dim; 
  PetscDS        ds;

  // add 
  PetscInt       n;

  PetscFunctionBegin;  

  PetscCall(user.ProcessOptions(comm));
  PetscCall(user.InputParameterOptions(comm, "granite"));

  user.SetMaterial(comm, "granite", 0, "domain rock");
  PetscCall(Mesh::Create(comm, &dm));
  PetscCall(SNESCreate(comm, &snes));
  PetscCall(SNESSetDM(snes, dm));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(user.SetupFE(dm, ds, "elasticity", &fe)); 
  PetscCall(user.SolveFE(dm, &displacement, "displacement", &snes));
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, NULL));
  if (dim == 3) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe, &fe_derived, 6, dim));}
  else if (dim == 2) {PetscCall(user.CreateDerivedField(dm, &dmAux, "elasticity_derived", &fe, &fe_derived, 4, dim));}

  user.SolveDerivedFields(dm, dmAux, &displacement, &derived_solutions[0]) ;
  
  /* Cleanup */
  PetscCall(VecDestroy(&displacement));
  PetscCall(VecDestroy(&derived_solutions[0]));
  PetscCall(VecDestroy(&derived_solutions[1]));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm));
  PetscCall(DMDestroy(&dmAux));
  PetscCall(PetscFEDestroy(&fe));
  PetscCall(PetscFEDestroy(&fe_derived));

  PetscFunctionReturn(0);
}
