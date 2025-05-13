#ifndef TRANSIENTPOROELASTICITY_HH
#define TRANSIENTPOROELASTICITY_HH

#include "PoroElasticity.hh"
#include "TransientDarcysFlow.hh"
#include "IsotropicLinearElasticity.hh"


class TransientPoroElasticity : public PoroElasticity {

protected:


public:
    TransientDarcysFlow transientdarcysflow;
    IsotropicLinearElasticity isotropiclinearelasticity; 

    static PetscErrorCode SolutionMonitor_HM(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx);

    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], TS* ts);
    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);        
    PetscErrorCode SetMaterial(MPI_Comm comm);
    PetscErrorCode DomainSetup(PetscDS ds) ;
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupTS(DM dm, Vec *sol, const char name[], TS* ts);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name) ;
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);

    TransientPoroElasticity();
    // virtual ~Elasticity();
};

#endif