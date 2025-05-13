#ifndef THERMAL_HH
#define THERMAL_HH

#include "../Base/petscheaders.hh"
#include "../PointwiseFunctions/PF_Thermal.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Base/Projection.hh"   // 
#include "../BoundaryCondition/BCThermal.hh"
#include "Problem.hh"

class Thermal : public Problem {


protected: 
//  change to private ??
    Material material;
    std::vector <Material> materials;

    const char *BCTypes[4] = {"insulation", "temperature", "heat_flux", "unknown"};
    std::vector<BCThermal> BoundaryConditions;    

public:
    enum DerivedSolutionType {HEATFLUX, NUM_DERIVED_SOLUTION_TYPES} ; /* soluion type */
    DerivedSolutionType derivedsoltype{HEATFLUX};
    const char *DerivedSolutionTypes[NUM_DERIVED_SOLUTION_TYPES+1] = {"heatflux", "unknown"};
    
    void AddMaterial(Material& mtrl);
    PetscInt GetNumofMaterials();
    PetscErrorCode SetupMaterial(DM dm, DM dmAux);
    PetscErrorCode SetupMaterial(DM dm, DM dmAux, PetscScalar time, PetscInt step, Vec *sol);
    PetscErrorCode SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[]);
    PetscErrorCode SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step);
    PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]);
    
    virtual PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name);
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols);

    Thermal();



};

#endif