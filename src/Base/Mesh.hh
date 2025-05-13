#ifndef MESH_HH
#define MESH_HH

#include "petscheaders.hh"
#include <vector>

namespace Mesh
{
    PetscErrorCode Create(MPI_Comm comm, DM *dm);
    PetscErrorCode BoundarySetup(DM *dm);
    PetscErrorCode Create_striplineload(MPI_Comm comm, DM *dm);
    PetscErrorCode BoundarySetup_striplineload(DM *dm);
    PetscErrorCode Create_cryer(MPI_Comm comm, DM *dm);
    PetscErrorCode BoundarySetup_cryer(DM *dm);
    
    
}

#endif