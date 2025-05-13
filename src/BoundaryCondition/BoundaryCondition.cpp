#include "BoundaryCondition.hh"

// PetscScalar BoundaryCondition::bdvalue = 0.0;
 
PetscErrorCode BoundaryCondition::SetBdIndex(PetscInt bd_index)
{
  PetscFunctionBegin;
  this->bdindex = bd_index;
  PetscFunctionReturn(0);  
}

PetscErrorCode BoundaryCondition::SetBdValue(PetscScalar value)
{
  PetscFunctionBegin;
  this->bdvalue = value;
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryCondition::SetBdValues(PetscScalar *values, PetscInt dimension)
{
  PetscFunctionBegin;
  for (int i = 0; i < dimension; i++) {this->bdvalues[i] = values[i];};
  PetscFunctionReturn(0);
}

PetscInt BoundaryCondition::GetBdIndex()
{
  return this->bdindex; 
}

PetscScalar BoundaryCondition::GetBdValue()
{
  return this->bdvalue;
}

// void BoundaryCondition::GetBdValue(PetscScalar &boundaryvalue)
// {
//   boundaryvalue = bdvalue;
// }

PetscScalar* BoundaryCondition::GetBdValues()
{
  static PetscScalar pvalues[3];
  // PetscScalar *pvalues = new PetscScalar[3];
  for (int i=0; i<3; i++) {
    pvalues[i] = this->bdvalues[i]; 
  } 
  return pvalues;
}

PetscErrorCode BoundaryCondition::SetBdArea(DM dm, DMLabel label, PetscInt bd_index)
{
  IS  is;
  const PetscInt  *faces;
  PetscInt Nf, f, *value;
  PetscScalar larea = 0.0, garea = 0.0, tarea = 0.0; 

  PetscFunctionBeginUser;
  PetscCall(DMLabelGetStratumIS(label, 1,  &is));

  if (is) {
    PetscCall(ISGetLocalSize(is, &Nf));
    PetscCall(ISGetIndices(is, &faces));
    for (f = 0; f < Nf; ++f) {
      PetscCall(DMLabelGetValue(label, faces[f], value));
      if (bd_index == *value) {
        PetscCall(DMPlexComputeCellGeometryFVM(dm, faces[f], &tarea, NULL, NULL));
        larea += tarea;
      }
    }
    PetscCall(ISRestoreIndices(is, &faces));
  }
  PetscCall(ISDestroy(&is));
  PetscCall(MPIU_Allreduce(&larea, &garea, 1, MPIU_SCALAR, MPIU_SUM, PetscObjectComm((PetscObject)dm)));
  this->area = garea;
  PetscFunctionReturn(0);

}

PetscErrorCode BoundaryCondition::SetBdValuePerArea(DM dm, DMLabel label, PetscInt bd_index)
{
  PetscFunctionBeginUser;

  PetscCall(SetBdArea(dm, label, bd_index));
  PetscCheck(this->area >= 0.0, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "The area or length of the boundary should be larger than 0.0");
  this->bdvalue /= area;

  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryCondition::SetBdValueTimesArea(DM dm, DMLabel label, PetscInt bd_index)
{
  PetscFunctionBeginUser;

  PetscCall(SetBdArea(dm, label, bd_index));
  PetscCheck(this->area >= 0.0, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "The area or length of the boundary should be larger than 0.0");
  this->bdvalue *= area;

  PetscFunctionReturn(0);
}

