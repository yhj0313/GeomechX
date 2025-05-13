#ifndef STEADYDARCYSFLOW_HH
#define STEADYDARCYSFLOW_HH

#include "PorousFluidFlow.hh"

class SteadyDarcysFlow : public PorousFluidFlow {

public:
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], SNES* snes);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
};


#endif