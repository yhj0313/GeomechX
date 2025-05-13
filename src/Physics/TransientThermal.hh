#ifndef TRANSIENTTHERMAL_HH
#define TRANSIENTTHERMAL_HH

#include "Thermal.hh"

class TransientThermal : public Thermal {

private:

public:

    PetscScalar init_temperature;
    std::string init_temperature_equation;
    PetscBool flg_iniv_equation;
    const char *inivalue_call_T[2] = {"-inivalue_temperature", "-inivalue_temperature_equation"};

    PetscBool expTS;      /* Flag for explicit timestepping */
    PetscBool lumped;     /* Lump the mass matrix */
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], TS* ts);
    virtual PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode SetupTS(DM dm, Vec *sol, const char name[], TS* ts);
    PetscErrorCode SetInitialTemperature(MPI_Comm comm);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode SetMaterial(MPI_Comm comm);
    PetscErrorCode SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode SetMaterials(MPI_Comm comm);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup_mat_const(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup_mat_fields(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup_TM_mat_const(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
    PetscErrorCode BCSetup_TM_mat_fields(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
    
    TransientThermal();
    
};

#endif