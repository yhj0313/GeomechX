#ifndef BOUNDARYCONDITION_HH
#define BOUNDARYCONDITION_HH

#include "../Base/petscheaders.hh"

class BoundaryCondition {
    public:
    PetscScalar bdvalue ;
    protected:
    PetscInt bdindex;
    PetscScalar bdvalues[3] = {0.0, 0.0, 0.0};
    PetscBool bdbool = PETSC_FALSE; 
    PetscScalar area = 0.0; // the total surface area of the boundary label, area in 3D and length in 2D.  
    // 
    public:
    virtual PetscErrorCode SetBCType(PetscInt bctype_index) = 0; 
    PetscErrorCode SetBdIndex(PetscInt bd_index);
    PetscErrorCode SetBdValue(PetscScalar value);
    PetscErrorCode SetBdValues(PetscScalar *values, PetscInt dimension);
    PetscInt GetBdIndex();
    PetscScalar GetBdValue();
    PetscScalar* GetBdValues();
    PetscErrorCode SetBdArea(DM dm, DMLabel label, PetscInt bd_index);
    PetscErrorCode SetBdValuePerArea(DM dm, DMLabel label, PetscInt bd_index);
    PetscErrorCode SetBdValueTimesArea(DM dm, DMLabel label, PetscInt bd_index);

};

#endif