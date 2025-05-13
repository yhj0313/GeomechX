#ifndef TRANSIENTDARCYSFLOW_HH
#define TRANSIENTDARCYSFLOW_HH

#include "PorousFluidFlow.hh"


static PetscErrorCode SolutionMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx);
    

class TransientDarcysFlow : public PorousFluidFlow {

private:
    
public:
    PetscScalar init_pressure;
    std::string init_pressure_equation;
    PetscBool flg_iniv_equation;
    const char *inivalue_call_H[2] = {"-inivalue_pressure", "-inivalue_pressure_equation"};

    PetscBool expTS;      /* Flag for explicit timestepping */
    PetscBool lumped;     /* Lump the mass matrix */

    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode SetInitialPressure(MPI_Comm comm);
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], TS* ts);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
    PetscErrorCode BCSetup_HM(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
    TransientDarcysFlow();
};


#endif