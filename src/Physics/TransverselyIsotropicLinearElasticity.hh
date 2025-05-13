#ifndef TRANSVERSELYISOTROPICLINEARELASTICITY_HH
#define TRANSVERSELYISOTROPICLINEARELASTICITY_HH

#include "Elasticity.hh"
#include <vector>

class TransverselyIsotropicLinearElasticity : public Elasticity{

    public:    
    enum InputParameterType {FIVEPARAMETERS, NUM_INPUT_PARA_TYPES} ; /* Type of two input elastic parameters, Five parameters are E, E', nu, nu', G'. */
    InputParameterType inputpara = FIVEPARAMETERS;
    const char *inputParameterTypes[NUM_INPUT_PARA_TYPES+1] = {"fiveparameters", "unknown"};
    
    PetscErrorCode SetMaterial (MPI_Comm comm, const std::string name, PetscInt label_number, std::string description);
    PetscErrorCode DomainSetup(PetscDS ds);
    PetscErrorCode SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe);
    PetscErrorCode InputParameterOptions(MPI_Comm comm, const std::string name);
    PetscErrorCode SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols);

};

#endif 