#ifndef STEADYTHERMOELASTICITY_HH
#define STEADYTHERMOELASTICITY_HH

#include "ThermoElasticity.hh"
#include "SteadyThermal.hh"
#include "IsotropicLinearElasticity.hh"


class SteadyThermoElasticity : public ThermoElasticity {

protected:
    SteadyThermal steadythermal;
    IsotropicLinearElasticity isotropiclinearelasticity; 


public:
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], SNES* snes);
    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);        
    PetscErrorCode DomainSetup(PetscDS ds) ;
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name) ;
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols);

};

#endif