#ifndef THERMOPOROELASTICITY_HH
#define THERMOPOROELASTICITY_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../PointwiseFunctions/PF_ThermoPoroElasticity.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Base/Projection.hh"   // 추가
#include "Problem.hh"

class ThermoPoroElasticity : public Problem {

protected:
    Material material;
    std::vector <Material> materials;
    
public:
   
    static PetscErrorCode CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace);
    void AddMaterial(Material& mtrl);
    PetscInt GetNumofMaterials();
    PetscErrorCode SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step);  
    PetscErrorCode SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step, Vec *sol);
    PetscErrorCode SetupMaterial(DM dm, DM dmAux, PetscScalar time, PetscInt step, Vec *sol) ;
    PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols);
    ThermoPoroElasticity();
    // virtual ~Elasticity();
};

#endif