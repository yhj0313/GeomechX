#include "BCThermal.hh"

BCThermal::BCThermal(PetscInt bctype_index) {
    this->bctype = (BCThermal::BoundaryConditionType)bctype_index;

}

PetscErrorCode BCThermal::SetBCType(PetscInt bctype_index){

    PetscFunctionBegin;
    this->bctype = (BCThermal::BoundaryConditionType)bctype_index;
    PetscFunctionReturn(0);  
}

BCThermal::BoundaryConditionType BCThermal::GetBCType(){
    return this->bctype;
}