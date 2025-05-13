#ifndef ELASTICITY_HH
#define ELASTICITY_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../PointwiseFunctions/PF_Elasticity.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../BoundaryCondition/BCElasticity.hh"
#include "Problem.hh"



class Elasticity : public Problem {

protected:
    Material material;
    
    const char *BCTypes[7] = {"free", "roller", "fixed", "displacement", "normal_stress", "traction", "unknown"};
    std::vector<BCElasticity> BoundaryConditions;

public:

    PetscBool    useNearNullspace; /* Use the rigid body modes as a near nullspace for AMG */ 
    enum DerivedSolutionType {PLANESTRAIN_2D, CAUCHY_3D, AXISYMMETRIC_2D_PLANESTRAIN, AXISYMMETRIC_3D, NUM_DERIVED_SOLUTION_TYPES} ; /* soluion type */
    DerivedSolutionType derivedsoltype{PLANESTRAIN_2D};
    PetscBool TIStripLineLoad;
    PetscBool cryer;
    
    const char *DerivedSolutionTypes[NUM_DERIVED_SOLUTION_TYPES+1] = {"planestrain", "cauchy_3d", "axisymmetric_2d_planestrain", "axisymmetric_3d", "unknown"};
    
   
    static PetscErrorCode CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace);
    PetscErrorCode SolveFE(DM dm, Vec *sol, const char name[], SNES* snes);
    PetscErrorCode CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim);
    PetscErrorCode SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]);
    PetscErrorCode ProcessOptions(MPI_Comm comm);
    PetscErrorCode ParameterSetup(PetscDS ds);
    PetscErrorCode SetupPrimalProblem(DM dm, PetscDS ds);
    PetscErrorCode CreateAuxiliaryVecv22(DM dm, DM *auxdm, Vec *la, std::vector<BCElasticity> bc);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf);
    PetscErrorCode BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field);
    DerivedSolutionType GetDerivedSolType();
    Elasticity();
    // virtual ~Elasticity();
};

#endif