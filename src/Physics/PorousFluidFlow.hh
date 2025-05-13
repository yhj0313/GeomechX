#ifndef POROUSFLUIDFLOW_HH
#define POROUSFLUIDFLOW_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../PointwiseFunctions/PF_PorousFluidFlow.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../BoundaryCondition/BCPorousFluidFlow.hh"
#include "Problem.hh"



class PorousFluidFlow: public Problem {

protected:
    Material material;
    
    const char *BCTypes[6] = {"noflow", "pressure", "massflux", "massflowrate", "darcyvelocity", "unknown"};
    std::vector<BCPorousFluidFlow> BoundaryConditions;

public:
    PetscBool    useGravity; 
    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name);
    PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) ;
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) ;
    
    PorousFluidFlow();
    // virtual ~Elasticity();
};

#endif

