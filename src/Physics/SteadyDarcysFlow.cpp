#include "SteadyDarcysFlow.hh"

PetscErrorCode SteadyDarcysFlow::SolveFE(DM dm, Vec *sol, const char name[], SNES* snes)
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 
  PetscPointFunc scale[2] = {Pw_Functions::scale_darcyvelocity, Pw_Functions::scale_pressure};
  Vec scaledsol;
  
  PetscFunctionBegin;

  PetscCall(DMCreateGlobalVector(dm, sol));
  PetscCall(PetscObjectSetName((PetscObject) *sol, name));
  PetscCall(VecSet(*sol, 0.0));
  PetscCall(DMPlexSetSNESLocalFEM(dm, NULL, NULL, NULL));
  PetscCall(SNESSetFromOptions(*snes));
  PetscCall(SNESSolve(*snes, NULL, *sol));
  PetscCall(SNESGetSolution(*snes, sol));
  PetscCall(DMCreateGlobalVector(dm, &scaledsol));
  PetscCall(VecSet(scaledsol, 0.0));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*sol, NULL, view_call));

  PetscCall(VecDestroy(&scaledsol));

  PetscFunctionReturn(0); 
}

PetscErrorCode SteadyDarcysFlow::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar permeability_darcy(1e-15) /* m^2 */, fluid_density_darcy(1000.0)/* kg/m^3 */, fluid_viscosity_darcy(0.001)/* Pa s */;
    PetscScalar gravity_darcy(0.0) ;/* m/s^2 */
    char option_call[PETSC_MAX_PATH_LEN];

    PetscFunctionBegin;   

    // std::cout << name << std::endl;
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_permeability", const_cast<char*>(name.c_str())));
    // std::cout << option_call << std::endl;
    PetscOptionsGetScalar(NULL, NULL, option_call, &permeability_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_density", const_cast<char*>(name.c_str())));
    // std::cout << option_call << std::endl;
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_viscosity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_viscosity_darcy, NULL);
    if (useGravity == PETSC_TRUE) {gravity_darcy = 9.8;}
    
    Property permeability(permeability_darcy, "m2", "Permeability" ); 
    Property fluid_density(fluid_density_darcy, "kg/m3", "Fluid density" ); 
    Property fluid_viscosity(fluid_viscosity_darcy, "Pa s", "Fluid viscosity" ); 
    Property gravity(gravity_darcy, "m/s^2", "Gravity" ); 
       
    material.SetName(name);
    material.SetLabel_number(label_number);
    material.SetDescription(description);
    material.AddProperty(permeability);
    material.AddProperty(fluid_density);
    material.AddProperty(fluid_viscosity);
    material.AddProperty(gravity);

    PetscFunctionReturn(0); 

}

PetscErrorCode SteadyDarcysFlow::ParameterSetup(PetscDS ds)
{
  const PetscInt num = material.GetNumofProperties();
  PetscScalar constants[num];

  PetscFunctionBeginUser;

  for (int i = 0; i < num; i++)
  {
      constants[i] = material.GetProperty(i);
  }

  constants[0] = material.GetProperty(0) * 1.0e18; /* scale permeability, 1.0e-18 m2 */
  constants[1] = material.GetProperty(1) * 1.0e-6; /* scale fluid density, 1.0e6 kg/m3 */
  constants[2] = material.GetProperty(2) * 1.0e3; /* scale fluid viscosity, 1.0e-3 Pa s */
  constants[3] = material.GetProperty(3); /* gravity */

  PetscCall(PetscDSSetConstants(ds, num, constants));

  PetscFunctionReturn(0);
}    

PetscErrorCode SteadyDarcysFlow::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  PetscCall(PetscDSSetResidual(ds, 0, Pw_Functions::f0_fluidflux, Pw_Functions::f1_fluidflux));
  PetscCall(PetscDSSetJacobian(ds, 0, 0, Pw_Functions::g0_fluidflux_flux, NULL, NULL, NULL));
  PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_fluidflux_pressure, NULL));
  PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_pressure, NULL));
  PetscCall(PetscDSSetJacobian(ds, 1, 0, NULL, Pw_Functions::g1_pressure_flux, NULL, NULL));
  PetscFunctionReturn(0);
}

PetscErrorCode SteadyDarcysFlow::SetupPrimalProblem(DM dm, PetscDS ds)
{
  // PetscDS          ds;
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));
  // PetscCall(ParameterSetup(ds));

  PetscFunctionReturn(0);
}


PetscErrorCode SteadyDarcysFlow::SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe)
{
  char           prefix[PETSC_MAX_PATH_LEN];
  PetscBool      simplex;
  PetscInt       dim;  
  PetscMPIInt    rank;


  PetscFunctionBegin;
  /* Create finite element */
  PetscCall(DMGetDimension(dm, &dim));
  DMPlexIsSimplex(dm, &simplex);

/* fluid flux */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[0]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, dim, simplex, name[0] ? prefix : NULL, -1, &fe[0]));
  PetscCall(PetscObjectSetName((PetscObject) fe[0], name[0]));
 
  /* pressure */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[1]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[1] ? prefix : NULL, -1, &fe[1]));
  PetscCall(PetscFECopyQuadrature(fe[0], fe[1]));
  PetscCall(PetscObjectSetName((PetscObject) fe[1], name[1]));
  
  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe[0]));
  PetscCall(DMSetField(dm, 1, NULL, (PetscObject) fe[1]));

  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  
  PetscCall(SetupPrimalProblem(dm, ds));
  // PetscCall(SetupPrimalProblem(dm, ds, fe, febd)); //added 230508
  PetscCall(ParameterSetup(ds));

  CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  if ((this->material_view) & (rank == 0)) PetscCall(material.ViewMaterial());
  
  PetscFunctionReturn(0);
}

PetscErrorCode SteadyDarcysFlow::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


PetscErrorCode SteadyDarcysFlow::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf) 
{
    PetscInt bd;
    const PetscInt id_2d[] = {1, 2, 3, 4};
    const PetscInt id_3d[] = {1, 2, 5, 6, 3, 4};
    const PetscInt cmp_2d[] = {0, 0, 1, 1};
    const PetscInt cmp_3d[] = {0, 0, 2, 2, 1, 1};
    // Pw_Functions       pf;
    const PetscInt components_2d[] = {0, 1};
    const PetscInt components_3d[] = {0, 1, 2};
    DMLabel label;
    PetscInt dim;
    PetscInt n_bd;
    PetscInt bc = 0;
    PetscScalar bdv = 0.0;
    PetscScalar *bdvs;
    PetscBool flg_bctype;
    PetscBool flg_bdvalue, flg_bdvalues;
    const char *boundary_name[6] = {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
    // const char *boundary_call[6] = {"-bc_type_left", "-bc_type_right", "-bc_type_bottom", "-bc_type_top", "-bc_type_front", "-bc_type_back"};
    
    const char *bctype_call[6] = {"-bctype_darcy_left", "-bctype_darcy_right", "-bctype_darcy_bottom", "-bctype_darcy_top", "-bctype_darcy_front", "-bctype_darcy_back"};
    const char *bdvalue_call[6] = {"-bdvalue_darcy_left", "-bdvalue_darcy_right", "-bdvalue_darcy_bottom", "-bdvalue_darcy_top", "-bdvalue_darcy_front", "-bdvalue_darcy_back"};
    const char *bdvalues_call[6] = {"-bdvalues_darcy_left", "-bdvalues_darcy_right", "-bdvalues_darcy_bottom", "-bdvalues_darcy_top", "-bdvalues_darcy_front", "-bdvalues_darcy_back"};
    
    void (*pressure[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_pressureboundary_fluidflux_bd_0, Pw_Functions::f0_pressureboundary_fluidflux_bd_1, Pw_Functions::f0_pressureboundary_fluidflux_bd_2, Pw_Functions::f0_pressureboundary_fluidflux_bd_3, Pw_Functions::f0_pressureboundary_fluidflux_bd_4, Pw_Functions::f0_pressureboundary_fluidflux_bd_5};

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCPorousFluidFlow(0)); /* Default boundary condition is the first option. (no flow) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");
    PetscCall(PetscMalloc1(dim, &bdvs));
    
    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}
    
for (int i = 0; i < n_bd; i++)
    {
        PetscCall(BoundaryConditions[i].SetBdIndex(i));
        PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCPorousFluidFlow::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
        if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        
        if (bc == 1) {  // pressure boundary
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 2) { //massflux 
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv /= material.GetProperty(1) * 1.0e-6; // divided by scaled fluid density
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 3) { // massflowrate
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv /= material.GetProperty(1) * 1.0e-6; // divided by scaled fluid density
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            PetscCall(BoundaryConditions[i].SetBdValuePerArea(dm, label, i));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 4) { // darcyvelocity 
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv *= 1.0e9;
            PetscPrintf(PETSC_COMM_WORLD, "bdv =  '%f' . \n", bdv);
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        }
        flg_bctype = PETSC_FALSE;
        bc = 0;
    }

    PetscCall(PetscFree(bdvs));
    PetscOptionsEnd();
    
 
    if (dim == 3)
    {
        for (int i = 0; i < 6; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCPorousFluidFlow::BoundaryConditionType::NOFLOW:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::PRESSURE:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, pressure[i], 0, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLOWRATE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::DARCYVELOCITY:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        default://NOFLOW
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCPorousFluidFlow::BoundaryConditionType::NOFLOW:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::PRESSURE:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, pressure[i], 0, NULL));  
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLOWRATE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::DARCYVELOCITY:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        default://NOFLOW
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        }
        }
    }

    PetscFunctionReturn(0);
}

PetscErrorCode SteadyDarcysFlow::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field) 
{
    PetscInt bd;
    const PetscInt id_2d[] = {1, 2, 3, 4};
    const PetscInt id_3d[] = {1, 2, 5, 6, 3, 4};
    const PetscInt cmp_2d[] = {0, 0, 1, 1};
    const PetscInt cmp_3d[] = {0, 0, 2, 2, 1, 1};
    // Pw_Functions       pf;
    const PetscInt components_2d[] = {0, 1};
    const PetscInt components_3d[] = {0, 1, 2};
    DMLabel label;
    PetscInt dim;
    PetscInt n_bd;
    PetscInt bc = 0;
    PetscScalar bdv = 0.0;
    PetscScalar *bdvs;
    PetscBool flg_bctype;
    PetscBool flg_bdvalue, flg_bdvalues;
    const char *boundary_name[6] = {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
    
    const char *bctype_call[6] = {"-bctype_darcy_left", "-bctype_darcy_right", "-bctype_darcy_bottom", "-bctype_darcy_top", "-bctype_darcy_front", "-bctype_darcy_back"};
    const char *bdvalue_call[6] = {"-bdvalue_darcy_left", "-bdvalue_darcy_right", "-bdvalue_darcy_bottom", "-bdvalue_darcy_top", "-bdvalue_darcy_front", "-bdvalue_darcy_back"};
    const char *bdvalues_call[6] = {"-bdvalues_darcy_left", "-bdvalues_darcy_right", "-bdvalues_darcy_bottom", "-bdvalues_darcy_top", "-bdvalues_darcy_front", "-bdvalues_darcy_back"};
    
    void (*pressure[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_pressureboundary_fluidflux_bd_0, Pw_Functions::f0_pressureboundary_fluidflux_bd_1, Pw_Functions::f0_pressureboundary_fluidflux_bd_2, Pw_Functions::f0_pressureboundary_fluidflux_bd_3, Pw_Functions::f0_pressureboundary_fluidflux_bd_4, Pw_Functions::f0_pressureboundary_fluidflux_bd_5};

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCPorousFluidFlow(0)); /* Default boundary condition is the first option. (no flow) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");

    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}
    
    for (int i = 0; i < n_bd; i++)
    {
        PetscCall(BoundaryConditions[i].SetBdIndex(i));
        PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCPorousFluidFlow::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
        if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        
        if (bc == 1) {  // pressure boundary
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 2) { //massflux 
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv /= material.GetProperty(1) * 1.0e-6; // divided by scaled fluid density
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 3) { // massflowrate
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv /= material.GetProperty(1) * 1.0e-6; // divided by scaled fluid density
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            PetscCall(BoundaryConditions[i].SetBdValuePerArea(dm, label, i));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 4) { // darcyvelocity 
            PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
            if (flg_bdvalue) {
            bdv *= 1.0e9;
            PetscPrintf(PETSC_COMM_WORLD, "bdv =  '%f' . \n", bdv);
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            }
            else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        }
        flg_bctype = PETSC_FALSE;
        bc = 0;
    }

    PetscCall(PetscFree(bdvs));
    PetscOptionsEnd();
    
 
    if (dim == 3)
    {
        for (int i = 0; i < 6; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCPorousFluidFlow::BoundaryConditionType::NOFLOW:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::PRESSURE:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, pressure[i], 0, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLOWRATE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::DARCYVELOCITY:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        default://NOFLOW
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCPorousFluidFlow::BoundaryConditionType::NOFLOW:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::PRESSURE:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, pressure[i], 0, NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::MASSFLOWRATE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCPorousFluidFlow::BoundaryConditionType::DARCYVELOCITY:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_value, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        default://NOFLOW
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero, NULL, NULL, NULL));
        break;
        
        }
        }
    }

    PetscFunctionReturn(0);
}
