#ifndef TRANSIENTTHERMOPOROELASTICITY_HH
#define TRANSIENTTHERMOPOROELASTICITY_HH

#include "ThermoPoroElasticity.hh"
#include "TransientThermal.hh"
#include "IsotropicLinearElasticity.hh"
#include "TransientDarcysFlow.hh"

class TransientThermoPoroElasticity : public ThermoPoroElasticity {

protected:
    

public:
    TransientThermal transientthermal;
    TransientDarcysFlow transientdarcysflow;
    IsotropicLinearElasticity isotropiclinearelasticity; 
    PetscBool expTS;      /* Flag for explicit timestepping */
    PetscBool lumped;     /* Lump the mass matrix */

    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step, PetscFE feAux[]) ;
    static PetscErrorCode SolutionMonitor_THM(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx);
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], TS* ts);
    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode SetupTS(DM dm, Vec *sol, const char name[], TS* ts) ;
    PetscErrorCode SetMaterial(MPI_Comm comm);
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);        
    PetscErrorCode SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description);
    PetscErrorCode SetMaterials(MPI_Comm comm);  
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode DomainSetup(PetscDS ds) ;
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name);
    TransientThermoPoroElasticity();
};

#endif