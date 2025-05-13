#ifndef POROELASTICITY_HH
#define POROELASTICITY_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../PointwiseFunctions/PF_PoroElasticity.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "Problem.hh"

class PoroElasticity : public Problem {

protected:
    Material material;
    
public:
   
    static PetscErrorCode CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace);
    PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim); 
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) ;
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) ;
    
    PoroElasticity();
    // virtual ~Elasticity();
};

#endif
