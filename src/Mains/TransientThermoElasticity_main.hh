#ifndef TRANSIENTTHERMOELASTICITY_MAIN_HH
#define TRANSIENTTHERMOELASTICITY_MAIN_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Physics/TransientThermoElasticity.hh"

PetscErrorCode main_transientthermoelasticity(MPI_Comm comm);

#endif