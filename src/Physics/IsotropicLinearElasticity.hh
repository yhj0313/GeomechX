#ifndef ISOTROPICLINEARELASTICITY_HH
#define ISOTROPICLINEARELASTICITY_HH

#include "Elasticity.hh"
#include <vector>

class IsotropicLinearElasticity : public Elasticity{

    public:    
    enum InputParameterType {YOUNGSMODULUS_POISSONSRATIO, YOUNGSMODULUS_SHEARMODULUS, BULKMODULUS_SHAERMODULUS, LAMESPARAMETERS, NUM_INPUT_PARA_TYPES} ; /* Type of two input elastic parameters */
    InputParameterType inputpara = LAMESPARAMETERS;
    const char *inputParameterTypes[NUM_INPUT_PARA_TYPES+1] = {"young_poisson", "young_shear", "bulk_shear", "lames", "unknown"};
    
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name);
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols);


};

#endif 