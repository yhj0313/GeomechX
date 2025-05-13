#ifndef PROBLEM_HH
#define PROBLEM_HH

#include "../Base/petscheaders.hh"

PetscErrorCode SetInitialConditions(TS ts, Vec u);
PetscErrorCode CreateSolveDerivedField(DM dm, PetscFE *fe, PetscPointFunc func, Vec *sol, const char name[], PetscReal time, PetscInt step, PetscInt numComp, PetscInt dim) ;
    
class Problem{
protected:
    PetscBool    material_view;
    PetscBool    material_fields = PETSC_FALSE;
    char         dmType[256]; /* DM type for the solve */
    
public:
    virtual PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim) = 0;
    virtual PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]) = 0;
    virtual PetscErrorCode ProcessOptions(MPI_Comm comm) = 0;
    virtual PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) = 0;
    virtual PetscErrorCode DomainSetup(PetscDS ds) = 0;
    virtual PetscErrorCode ParameterSetup(PetscDS ds) = 0;
    virtual PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds) = 0;
    virtual PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe) = 0;
    virtual PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name) = 0;
    virtual PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) = 0;
    virtual PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf) = 0;
/*     DM  dm;
    PetscBool    isTimedependent;
    PetscInt    dim;
    // @todo: add range of time? */
};

#endif