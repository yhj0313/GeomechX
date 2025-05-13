#include "BCElasticity.hh"

BCElasticity::BCElasticity(PetscInt bctype_index) {
    this->bctype = (BCElasticity::BoundaryConditionType)bctype_index;

}

PetscErrorCode BCElasticity::SetBCType(PetscInt bctype_index){

    PetscFunctionBegin;
    this->bctype = (BCElasticity::BoundaryConditionType)bctype_index;
    PetscFunctionReturn(0);  
}

BCElasticity::BoundaryConditionType BCElasticity::GetBCType(){
    return this->bctype;
}


