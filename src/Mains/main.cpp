static char help[] = "GeomechX v0.0.\n\n\n";


#include "IsotropicLinearElasticity_main.hh"
#include "TransverselyIsotropicLinearElasticity_main.hh"
#include "SteadyThermal_main.hh"
#include "TransientThermal_main.hh"
#include "SteadyThermoElasticity_main.hh"
#include "TransientThermoElasticity_main.hh"
#include "SteadyDarcysFlow_main.hh"
#include "TransientDarcysFlow_main.hh"
#include "TransientPoroElasticity_main.hh"
#include "TransientThermoPoroElasticity_main.hh"


int main(int argc, char **argv) {

    PetscFunctionBegin;
    PetscCall(PetscOptionsSetValue(NULL, "-skip_petscrc", "true"));
    // PetscCall(PetscInitialize(&argc, &argv, NULL, help)); 
    PetscCall(PetscInitialize(&argc, &argv, ".geomechxrc", help)); // default option file name : .geomechxrc

    PetscMPIInt rank;
    // MPI_Comm comm;
    PetscInt physics = 0;
    PetscBool flg;
    enum PhysicsType
    {
        ISOTROPICLINEARELASTICITY,
        TRANSVERSELYISOTROPICLINEARELASTICITY,
        STEADYTHERMAL,
        TRANSIENTTHERMAL,
        STEADYTHERMOELASTICITY,
        TRANSIENTTHERMOELASTICITY,
        STEADYDARCYSFLOW,
        TRANSIENTDARCYSFLOW,
        TRANSIENTPOROELASTICITY,
        TRANSIENTTHERMOPOROELASTICITY,
        NUM_INPUT_PHYSICS_TYPES
    }; /* Type of two input elastic parameters */

    PhysicsType physicstype = NUM_INPUT_PHYSICS_TYPES;
    const char *PhysicsTypes[NUM_INPUT_PHYSICS_TYPES + 1] = {"IsotropicLinearElasticity", "TransverselyIsotropicLinearElasticity", "SteadyThermal", "TransientThermal", "SteadyThermoElasticity", "TransientThermoElasticity", "SteadyDarcysFlow", "TransientDarcysFlow", "TransientPoroElasticity", "TransientThermoPoroElasticity", "unknown"};

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Type of physics", "");
    physics = physicstype;
    PetscCall(PetscOptionsEList("-physics_type", "Type of physics", NULL, PhysicsTypes, NUM_INPUT_PHYSICS_TYPES, PhysicsTypes[physicstype], &physics, &flg));
    if (flg)
    {
      physicstype = (PhysicsType)physics;

    }
    else
    {
      CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
      if (rank == 0)
      {
          std::cout << "ERROR: Physics is not given." << std::endl;
          break;
      }
    }
    
    PetscOptionsEnd();

    if (flg)
    {
      switch (physicstype)
      {
      case ISOTROPICLINEARELASTICITY:
          PetscCall(main_isotropiclinearelasticity(PETSC_COMM_WORLD));
          break;
      case TRANSVERSELYISOTROPICLINEARELASTICITY:
          PetscCall(main_transverselyisotropiclinearelasticity(PETSC_COMM_WORLD));
          break;
      case STEADYTHERMAL:
          PetscCall(main_steadythermal(PETSC_COMM_WORLD));
          break;
      case TRANSIENTTHERMAL:
          PetscCall(main_transientthermal(PETSC_COMM_WORLD));
          break;
      case STEADYTHERMOELASTICITY:
          PetscCall(main_steadythermoelasticity(PETSC_COMM_WORLD));
          break;
      case TRANSIENTTHERMOELASTICITY:
          PetscCall(main_transientthermoelasticity(PETSC_COMM_WORLD));
          break;
      case STEADYDARCYSFLOW:
          PetscCall(main_steadydarcysflow(PETSC_COMM_WORLD));
          break;
      case TRANSIENTDARCYSFLOW:
          PetscCall(main_transientdarcysflow(PETSC_COMM_WORLD));
          break;    
      case TRANSIENTPOROELASTICITY:
          PetscCall(main_transientporoelasticity(PETSC_COMM_WORLD));
          break;
      case TRANSIENTTHERMOPOROELASTICITY:
          PetscCall(main_transientthermoporoelasticity(PETSC_COMM_WORLD));
          break;
      }
    }

    PetscCall(PetscFinalize());
    PetscFunctionReturn(0);
}