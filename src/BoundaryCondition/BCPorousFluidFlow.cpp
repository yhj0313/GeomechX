#include "BCPorousFluidFlow.hh"

BCPorousFluidFlow::BCPorousFluidFlow(PetscInt bctype_index) {
    this->bctype = (BCPorousFluidFlow::BoundaryConditionType)bctype_index;

}

PetscErrorCode BCPorousFluidFlow::SetBCType(PetscInt bctype_index){

    PetscFunctionBegin;
    this->bctype = (BCPorousFluidFlow::BoundaryConditionType)bctype_index;
    PetscFunctionReturn(0);  
}

BCPorousFluidFlow::BoundaryConditionType BCPorousFluidFlow::GetBCType(){
    return this->bctype;
}