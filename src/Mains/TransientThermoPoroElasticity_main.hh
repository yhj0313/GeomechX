#ifndef TRANSIENTTHERMOPOROELASTICITY_MAIN_HH
#define TRANSIENTTHERMOPOROELASTICITY_MAIN_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Physics/TransientThermoPoroElasticity.hh"


PetscErrorCode main_transientthermoporoelasticity(MPI_Comm comm);

#endif