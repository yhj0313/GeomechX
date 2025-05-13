#ifndef BCTHERMAL_HH
#define BCTHERMAL_HH

#include "BoundaryCondition.hh"

class BCThermal : public BoundaryCondition {
public:
    enum BoundaryConditionType {INSULATION, TEMPERATURE, HEAT_FLUX, NUM_BC_TYPES};
    // const char *BCTypes[4] = {"insulation", "temperature", "heat_flux", "unknown"};
    BCThermal(PetscInt bctype_index);
    virtual PetscErrorCode SetBCType(PetscInt bctype_index); 
    BCThermal::BoundaryConditionType GetBCType();

protected:    
    BoundaryConditionType bctype;
};

#endif 

