#ifndef BCPOROUSFLUIDFLOW_HH
#define BCPOROUSFLUIDFLOW_HH

#include "BoundaryCondition.hh"

class BCPorousFluidFlow : public BoundaryCondition {
public:
    enum BoundaryConditionType {NOFLOW, PRESSURE, MASSFLUX, MASSFLOWRATE, DARCYVELOCITY, NUM_BC_TYPES};
    BCPorousFluidFlow(PetscInt bctype_index);
    virtual PetscErrorCode SetBCType(PetscInt bctype_index); 
    BCPorousFluidFlow::BoundaryConditionType GetBCType();

protected:    
    BoundaryConditionType bctype;
};

#endif 