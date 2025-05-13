#ifndef STEADYTHERMAL_HH
#define STEADYTHERMAL_HH

#include "Thermal.hh"

class SteadyThermal : public Thermal {

public:

    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], SNES* snes);
    PetscErrorCode SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode SetMaterial(MPI_Comm comm);
    PetscErrorCode SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description);
    PetscErrorCode SetMaterials(MPI_Comm comm);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup_mat_const(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup_mat_fields(DM dm, PetscDS ds, PetscWeakForm wf);    
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);


};

#endif