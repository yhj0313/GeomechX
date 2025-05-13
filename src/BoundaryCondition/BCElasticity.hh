#ifndef BCELASTICITY_HH
#define BCELASTICITY_HH

#include "BoundaryCondition.hh"
// #include "../Physics/Elasticity.hh"

class BCElasticity : public BoundaryCondition {
public:
    enum BoundaryConditionType {FREE, ROLLER, FIXED, DISPLACEMENT, NORMAL_STRESS, TRACTION, NUM_BC_TYPES};
    // const char *BCTypes[NUM_BC_TYPES + 1] = {"free", "roller", "fixed", "displacement", "normal_stress", "traction", "unknown"};
    BCElasticity(PetscInt bctype_index);
    virtual PetscErrorCode SetBCType(PetscInt bctype_index); 
    BCElasticity::BoundaryConditionType GetBCType();
    static void f0_normal_stress_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]);
    BCElasticity();
    friend PetscErrorCode PetscWeakFormSetIndexBdResidual(PetscWeakForm, DMLabel, PetscInt, PetscInt, PetscInt, PetscInt, void (*)(PetscInt, PetscInt, PetscInt, const PetscInt *, const PetscInt *, const PetscScalar *, const PetscScalar *, const PetscScalar *, const PetscInt *, const PetscInt *, const PetscScalar *, const PetscScalar *, const PetscScalar *, PetscReal, const PetscReal *, const PetscReal *, PetscInt, const PetscScalar *, PetscScalar *), PetscInt, void (*)(PetscInt, PetscInt, PetscInt, const PetscInt *, const PetscInt *, const PetscScalar *, const PetscScalar *, const PetscScalar *, const PetscInt *, const PetscInt *, const PetscScalar *, const PetscScalar *, const PetscScalar *, PetscReal, const PetscReal *, const PetscReal *, PetscInt, const PetscScalar *, PetscScalar *));

protected:    
    BoundaryConditionType bctype;
};

#endif