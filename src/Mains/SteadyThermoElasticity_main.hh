#ifndef STEADYTHERMOELASTICITY_MAIN_HH
#define STEADYTHERMOELASTICITY_MAIN_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Physics/SteadyThermoElasticity.hh"
#include "../PointwiseFunctions/PF_Thermal.hh"


PetscErrorCode main_steadythermoelasticity(MPI_Comm comm);

#endif