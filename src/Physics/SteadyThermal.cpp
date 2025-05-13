#include "SteadyThermal.hh"


PetscErrorCode SteadyThermal::SolveFE(DM dm, Vec *sol, const char name[], SNES* snes)
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;

  PetscCall(DMCreateGlobalVector(dm, sol));
  PetscCall(PetscObjectSetName((PetscObject) *sol, name));
  PetscCall(VecSet(*sol, 0.0));
  PetscCall(DMPlexSetSNESLocalFEM(dm, NULL, NULL, NULL));
  PetscCall(SNESSetFromOptions(*snes));
  PetscCall(SNESSolve(*snes, NULL, *sol));
  PetscCall(SNESGetSolution(*snes, sol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*sol, NULL, view_call));

  PetscFunctionReturn(0); 
}



PetscErrorCode SteadyThermal::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar ther_cond_temp(3.5), fluid_density_temp(1000.0), heat_capa_temp(4182.0), flow_vel_temp(0.0);
    char option_call[PETSC_MAX_PATH_LEN];

    PetscFunctionBegin;   

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_thermal_conductivity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &ther_cond_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_density", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_heat_capacity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_flow_velocity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);

    
    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" ); 
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property flow_vel(flow_vel_temp, "m/s", "Flow velocity" ); 
       
    material.SetName(name);
    material.SetLabel_number(label_number);
    material.SetDescription(description);
    material.AddProperty(ther_cond);
    material.AddProperty(fluid_density);
    material.AddProperty(heat_capa);
    material.AddProperty(flow_vel);

    PetscFunctionReturn(0); 

}

PetscErrorCode SteadyThermal::SetMaterial(MPI_Comm comm) {

    PetscScalar ther_cond_temp(3.5), fluid_density_temp(1000.0), heat_capa_temp(4182.0), flow_vel_temp(0.0);
    char option_call[PETSC_MAX_PATH_LEN];
    char material_name[128];
    char material_description[128];
    PetscBool flg_material_name, flg_material_description;
    std::string name(""), description("");
    PetscInt label_number(0);

    PetscFunctionBegin;   

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-thermal_conductivity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &ther_cond_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_density"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-heat_capacity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-flow_velocity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-name"));
    PetscOptionsGetString(NULL, NULL, option_call, material_name, sizeof(material_name), &flg_material_name);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-description"));
    PetscOptionsGetString(NULL, NULL, option_call, material_description, sizeof(material_description), &flg_material_description);

    if (flg_material_name) name = material_name; 
    if (flg_material_description) description = material_description;

    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" ); 
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property flow_vel(flow_vel_temp, "m/s", "Flow velocity" ); 
       
    material.SetName(name);
    material.SetLabel_number(label_number);
    material.SetDescription(description);
    material.AddProperty(ther_cond);
    material.AddProperty(fluid_density);
    material.AddProperty(heat_capa);
    material.AddProperty(flow_vel);

    PetscFunctionReturn(0); 

}


PetscErrorCode SteadyThermal::SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description) {
    //changed the last varialbe of PetscSNPrintf from name to field number
    PetscScalar ther_cond_temp(3.5), fluid_density_temp(1000.0), heat_capa_temp(4182.0), flow_vel_temp(0.0);
    char option_call[PETSC_MAX_PATH_LEN];
    char material_name[128];
    char material_description[128];
    PetscBool flg_material_name, flg_material_description;
    std::string name("");
    Material mtrl;
    

    PetscFunctionBegin;   
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-thermal_conductivity_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &ther_cond_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_density_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-heat_capacity_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-flow_velocity_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-name_%d", label_number));
    PetscOptionsGetString(NULL, NULL, option_call, material_name, sizeof(material_name), &flg_material_name);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-description_%d", label_number));
    PetscOptionsGetString(NULL, NULL, option_call, material_description, sizeof(material_description), &flg_material_description);

    if (!flg_material_name) {name = std::to_string(label_number);}
    else {
      name = material_name; 
      name.append("_");
      name.append(std::to_string(label_number));
      }

    if (flg_material_description) description = material_description;

    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" ); 
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property flow_vel(flow_vel_temp, "m/s", "Flow velocity" ); 
       
    mtrl.SetName(name);
    mtrl.SetLabel_number(label_number);
    mtrl.SetDescription(description);
    mtrl.AddProperty(ther_cond);
    mtrl.AddProperty(fluid_density);
    mtrl.AddProperty(heat_capa);
    mtrl.AddProperty(flow_vel);

    AddMaterial(mtrl);

    PetscFunctionReturn(0); 

}

PetscErrorCode SteadyThermal::SetMaterials(MPI_Comm comm) {
  
  PetscInt num_materials(1);

  PetscFunctionBegin;

  if (this->material_fields) {
    PetscOptionsGetInt(NULL, NULL, "-num_materials", &num_materials, NULL);
    for (int i = 0; i < num_materials; i++) {
      PetscCall(SetMaterial(comm, i, ""));
    }
  }
  else {PetscCall(SetMaterial(comm)); }

  
  PetscFunctionReturn(0); 

}


PetscErrorCode SteadyThermal::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  if (!material_fields) {
    PetscCall(PetscDSSetResidual(ds, 0, Pw_Functions::f0_temp, Pw_Functions::f1_temp));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, Pw_Functions::g1_temp, NULL, Pw_Functions::g3_temp));
  }
  else {
    PetscCall(PetscDSSetResidual(ds, 0, Pw_Functions::f0_temp_mat_fields, Pw_Functions::f1_temp_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, Pw_Functions::g1_temp_mat_fields, NULL, Pw_Functions::g3_temp_mat_fields));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermal::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));
  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermal::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  DM             cdm  = dm;
  char           prefix[PETSC_MAX_PATH_LEN];
  PetscBool      simplex;
  PetscInt       dim;  
  PetscMPIInt    rank;
  PetscInt properties_num = 0;

  PetscFunctionBegin;
  /* Create finite element */
  PetscCall(DMGetDimension(dm, &dim));
  DMPlexIsSimplex(dm, &simplex);
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name ? prefix : NULL, -1, fe));
  PetscCall(PetscObjectSetName((PetscObject) *fe, name));

  /* Create finite element for material properties */
  if (this->material_fields) {
    properties_num = this->materials[0].GetNumofProperties();
  }
  PetscFE  feAux[properties_num];
//option bool을 넣고, true 이면 아래 수정해서 property 수만큼 수행. q
  if (this->material_fields) {
    // properties_num = this->materials[0].GetNumofProperties();
    // PetscFE  feAux[properties_num];
    for (int i = 0; i < properties_num ; i++) {
      PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, "Material_", PETSC_DEFAULT, &feAux[i]));
      PetscCall(PetscFECopyQuadrature(fe[0], feAux[i]));
    } 
  }

  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) *fe));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  
  PetscCall(SetupPrimalProblem(dm, ds));

  if (!this->material_fields) {
    PetscCall(ParameterSetup(ds));
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if ((this->material_view) & (rank == 0)) PetscCall(material.ViewMaterial());
  }
  else {
    while (cdm) {
      PetscCall(SetupAuxDM(cdm, properties_num, feAux));
      PetscCall(DMCopyDisc(dm, cdm));
      PetscCall(DMGetCoarseDM(cdm, &cdm));
    }
  }

  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFEDestroy(&feAux[i]));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermal::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
{
  PetscFunctionBeginUser;

  if (material_fields) PetscCall(BCSetup_mat_fields(dm, ds, wf));
  else PetscCall(BCSetup_mat_const(dm, ds, wf));

  PetscFunctionReturn(0);

}

PetscErrorCode SteadyThermal::BCSetup_mat_const(DM dm, PetscDS ds, PetscWeakForm wf)
{

    PetscInt bd;
    const PetscInt id_2d[] = {1, 2, 3, 4};
    const PetscInt id_3d[] = {1, 2, 5, 6, 3, 4};
    const PetscInt cmp_2d[] = {0, 0, 1, 1};
    const PetscInt cmp_3d[] = {0, 0, 2, 2, 1, 1};
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
    
    const char *bctype_call[6] = {"-bctype_heat_left", "-bctype_heat_right", "-bctype_heat_bottom", "-bctype_heat_top", "-bctype_heat_front", "-bctype_heat_back"};
    const char *bdvalue_call[6] = {"-bdvalue_heat_left", "-bdvalue_heat_right", "-bdvalue_heat_bottom", "-bdvalue_heat_top", "-bdvalue_heat_front", "-bdvalue_heat_back"};
    const char *bdvalues_call[6] = {"-bdvalues_heat_left", "-bdvalues_heat_right", "-bdvalues_heat_bottom", "-bdvalues_heat_top", "-bdvalues_heat_front", "-bdvalues_heat_back"};
    
    void (*heatflux[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_heatflux_0, Pw_Functions::f0_heatflux_1, Pw_Functions::f0_heatflux_2, Pw_Functions::f0_heatflux_3, Pw_Functions::f0_heatflux_4, Pw_Functions::f0_heatflux_5};

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCThermal(0)); /* Default boundary condition is the first option. (insulation) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");
    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}
    
    for (int i = 0; i < n_bd; i++)
    {
      PetscCall(BoundaryConditions[i].SetBdIndex(i));
      PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCThermal::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
      if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));

        
        if (bc == 1) {  // temperature boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 2) { // heat flux boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
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
        case BCThermal::BoundaryConditionType::INSULATION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_insulation, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation, 0, NULL, 0, NULL)); 
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        *BDVALUE = BoundaryConditions[i].GetBdValue();
        {PetscPrintf(PETSC_COMM_WORLD, "BDVALUE is %f \n", *BDVALUE);}
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_heatflux, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_insulation, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCThermal::BoundaryConditionType::INSULATION:
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));        
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, (void*)(&BoundaryConditions[i]), &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, heatflux[i], 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, Pw_Functions::f0_insulation, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SteadyThermal::BCSetup_mat_fields(DM dm, PetscDS ds, PetscWeakForm wf)
{
    // SteadyThermal test;
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
    const char *bctype_call[6] = {"-bctype_heat_left", "-bctype_heat_right", "-bctype_heat_bottom", "-bctype_heat_top", "-bctype_heat_front", "-bctype_heat_back"};
    const char *bdvalue_call[6] = {"-bdvalue_heat_left", "-bdvalue_heat_right", "-bdvalue_heat_bottom", "-bdvalue_heat_top", "-bdvalue_heat_front", "-bdvalue_heat_back"};
    const char *bdvalues_call[6] = {"-bdvalues_heat_left", "-bdvalues_heat_right", "-bdvalues_heat_bottom", "-bdvalues_heat_top", "-bdvalues_heat_front", "-bdvalues_heat_back"};
    
    void (*heatflux[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_heatflux_0_mat_fields, Pw_Functions::f0_heatflux_1_mat_fields, Pw_Functions::f0_heatflux_2_mat_fields, Pw_Functions::f0_heatflux_3_mat_fields, Pw_Functions::f0_heatflux_4_mat_fields, Pw_Functions::f0_heatflux_5_mat_fields};

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCThermal(0)); /* Default boundary condition is the first option. (insulation) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");
    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}
    
    for (int i = 0; i < n_bd; i++)
    {
      PetscCall(BoundaryConditions[i].SetBdIndex(i));
      PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCThermal::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
      if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        
        if (bc == 1) {  // temperature boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 2) { // heat flux boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            std::cout << BoundaryConditions[i].GetBdValue() << std::endl;
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
        case BCThermal::BoundaryConditionType::INSULATION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_insulation_mat_fields, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_mat_fields, 0, NULL, 0, NULL)); 
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        *BDVALUE = BoundaryConditions[i].GetBdValue();
        {PetscPrintf(PETSC_COMM_WORLD, "BDVALUE is %f \n", *BDVALUE);}
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_heatflux_mat_fields, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_mat_fields, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_insulation_mat_fields, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_mat_fields, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCThermal::BoundaryConditionType::INSULATION:
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));        
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, (void*)(&BoundaryConditions[i]), &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, heatflux[i], 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_mat_fields, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, Pw_Functions::f0_insulation_mat_fields, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], 0, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_mat_fields, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SteadyThermal::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field)
{/* for TM boundary conditions */
    PetscInt bd;
    const PetscInt id_2d[] = {1, 2, 3, 4};
    const PetscInt id_3d[] = {1, 2, 5, 6, 3, 4};
    const PetscInt cmp_2d[] = {0, 0, 1, 1};
    const PetscInt cmp_3d[] = {0, 0, 2, 2, 1, 1};
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
    const char *bctype_call[6] = {"-bctype_heat_left", "-bctype_heat_right", "-bctype_heat_bottom", "-bctype_heat_top", "-bctype_heat_front", "-bctype_heat_back"};
    const char *bdvalue_call[6] = {"-bdvalue_heat_left", "-bdvalue_heat_right", "-bdvalue_heat_bottom", "-bdvalue_heat_top", "-bdvalue_heat_front", "-bdvalue_heat_back"};
    const char *bdvalues_call[6] = {"-bdvalues_heat_left", "-bdvalues_heat_right", "-bdvalues_heat_bottom", "-bdvalues_heat_top", "-bdvalues_heat_front", "-bdvalues_heat_back"};
    
    void (*heatflux[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_heatflux_0_TM, Pw_Functions::f0_heatflux_1_TM, Pw_Functions::f0_heatflux_2_TM, Pw_Functions::f0_heatflux_3_TM, Pw_Functions::f0_heatflux_4_TM, Pw_Functions::f0_heatflux_5_TM};

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));


    BoundaryConditions.assign(6, BCThermal(0)); /* Default boundary condition is the first option. (insulation) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");
    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}
    
    for (int i = 0; i < n_bd; i++)
    {
      PetscCall(BoundaryConditions[i].SetBdIndex(i));
      PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCThermal::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
      if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        
        if (bc == 1) {  // temperature boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 2) { // heat flux boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
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
        case BCThermal::BoundaryConditionType::INSULATION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, Pw_Functions::f0_insulation_TM, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, heatflux[i], 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, Pw_Functions::f0_insulation_TM, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_3d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCThermal::BoundaryConditionType::INSULATION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, Pw_Functions::f0_insulation_TM, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        case BCThermal::BoundaryConditionType::TEMPERATURE:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))Pw_Functions::boundary_temperature, NULL, (void*)(&BoundaryConditions[i]), NULL));        
        break;
        case BCThermal::BoundaryConditionType::HEAT_FLUX:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, (void*)(&BoundaryConditions[i]), &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, heatflux[i], 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        default://insulation
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, Pw_Functions::f0_insulation_TM, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdJacobian(wf, label, id_2d[i], field, 0, 0, 0, NULL, 0, Pw_Functions::g1_insulation_TM, 0, NULL, 0, NULL)); 
        break;
        }
        }
    }

    PetscFunctionReturn(0);
}

