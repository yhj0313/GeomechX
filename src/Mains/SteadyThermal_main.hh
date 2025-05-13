#ifndef STEADYTHERMAL_MAIN_HH
#define STEADYTHERMAL_MAIN_HH

#include "../Base/petscheaders.hh"
#include "../Base/Property.hh"
#include "../Base/Mesh.hh"
#include "../Base/Material.hh"
#include "../Physics/SteadyThermal.hh"
#include "../PointwiseFunctions/PF_Thermal.hh"


PetscErrorCode main_steadythermal(MPI_Comm comm);

#endif