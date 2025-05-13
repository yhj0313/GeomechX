#include "TransientThermoPoroElasticity.hh"


PetscErrorCode TransientThermoPoroElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) {
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, time, *sol, &func, INSERT_VALUES, *derivedsol));
  PetscCall(DMSetOutputSequenceNumber(*dmAux, step, time));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  


  PetscFunctionReturn(0);  
}

PetscErrorCode TransientThermoPoroElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step, PetscFE feAux[]) {
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscInt properties_num_t = this->materials[0].GetNumofProperties();

  PetscCall(DMCopyAuxiliaryVec(dm, *dmAux));
  PetscCall(DMProjectField(*dmAux, time, *sol, &func, INSERT_VALUES, *derivedsol));

  PetscCall(DMSetOutputSequenceNumber(*dmAux, step, time));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  

  PetscFunctionReturn(0);  
}

PetscErrorCode TransientThermoPoroElasticity::SolutionMonitor_THM(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
  PetscBool      simplex;
  PetscInt       dim, comp;
  PetscDS        ds;
  PetscFE        fe[4];  //, fe_derived;
  DM             dm;//, dmAux;

  TransientThermoPoroElasticity* ttpe = (TransientThermoPoroElasticity*)ctx; 
  PetscInt properties_num(0);

/* displacement, temperature, heat flux, total stress, effective stress, strain  */

  PetscFunctionBeginUser;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMGetDS(dm, &ds));

  PetscCall(DMGetField(dm, 0, NULL, (PetscObject *)&fe[0])); /* displacement */
  PetscCall(DMGetField(dm, 1, NULL, (PetscObject *)&fe[1])); /* temperature */
  PetscCall(DMGetField(dm, 2, NULL, (PetscObject *)&fe[2])); /* pore pressure */
  PetscCall(DMGetField(dm, 3, NULL, (PetscObject *)&fe[3])); /* volumetric strain */

  
  if (ttpe->material_fields) {
  properties_num = ttpe->materials[0].GetNumofProperties();
  }
  PetscFE  feAux[properties_num];
  PetscCall(DMGetDimension(dm, &dim));
  DMPlexIsSimplex(dm, &simplex);

  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, "Material_", PETSC_DEFAULT, &feAux[i]));
    PetscCall(PetscFECopyQuadrature(fe[0], feAux[i]));
    PetscCall(PetscObjectSetName((PetscObject) feAux[i], const_cast<char*>(ttpe->materials[0].GetPropertyDescription(i).c_str())));
  } 
 
  PetscCall(ttpe->SetupAuxDM(dm, properties_num, feAux, time, step, &u));
  
  PetscInt remain;
  remain = step % 1;

  if (remain == 0 || 3.1536e8 - 3.1536e7 < time && time < 3.1536e8 + 3.1536e7 || 3.1536e9 - 3.1536e7 < time && time < 3.1536e9 + 3.1536e7 || 3.1536e10 - 3.1536e7 < time && time < 3.1536e10 + 3.1536e7 || 3.1536e11 - 3.1536e7 < time && time < 3.1536e11 + 3.1536e7 || time == 3.1536e12) {
    PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::displacement_THM, &u, "displacement", time, step, dim, dim));
    PetscCall(CreateSolveDerivedField(dm, &fe[1], Pw_Functions::temperature_THM, &u, "temperature", time, step, 1, dim));
    PetscCall(CreateSolveDerivedField(dm, &fe[2], Pw_Functions::porepressure_THM, &u, "pressure", time, step, 1, dim));
    if (ttpe->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::heatflux_THM_mat_fields, &u, "heat_flux", time, step, dim, dim));
    else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::heatflux_THM, &u, "heat_flux", time, step, dim, dim));
    
    if (ttpe->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::darcyvelocity_THM_mat_fields, &u, "darcyvelocity", time, step, dim, dim));
    else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::darcyvelocity_THM, &u, "darcyvelocity", time, step, dim, dim));

    if (dim == 2) {comp = 4;}
    else if (dim == 3) {comp = 6;}
/* stress / strain (# of component = 4 or 6) */
    switch (ttpe->isotropiclinearelasticity.derivedsoltype) {
      case Elasticity::DerivedSolutionType::PLANESTRAIN_2D:
      if (ttpe->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_planestrain_THM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_planestrain_THM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystrain_planestrain_THM, &u, "strain", time, step, comp, dim));
      break;
      case Elasticity::DerivedSolutionType::CAUCHY_3D:
      if (ttpe->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_3D_THM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_3D_THM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystrain_3D_THM, &u, "strain", time, step, comp, dim));
      break;
      case Elasticity::DerivedSolutionType::AXISYMMETRIC_2D_PLANESTRAIN:
      if (ttpe->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_stress_planestrain_THM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_stress_planestrain_THM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_strain_planestrain_THM, &u, "strain", time, step, comp, dim));
      break;
      default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type for elasticity");
    }   
  }
  
  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFEDestroy(&feAux[i]));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoPoroElasticity::SolveFE(DM dm, Vec *sol, const char name[], TS* ts)
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;

  PetscCall(TSSolve(*ts, *sol));
  PetscCall(DMTSCheckFromOptions(*ts, *sol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*sol, NULL, view_call));

  PetscFunctionReturn(0); 
}


PetscErrorCode TransientThermoPoroElasticity::ProcessOptions(MPI_Comm comm) {
    
  PetscInt       derived_sol_elas_type = 0, derived_sol_heat_type = 0, rank; 
  PetscBool      flg_elas, flg_heat;
  
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-near_nullspace", "Use the rigid body modes as an AMG near nullspace", NULL, this->isotropiclinearelasticity.useNearNullspace, &this->isotropiclinearelasticity.useNearNullspace, NULL)); 
  PetscCall(PetscOptionsBool("-explicit", "Use explicit timestepping", NULL, this->expTS, &this->expTS, NULL));
  PetscCall(PetscOptionsBool("-lumped", "Lump the mass matrix", NULL, this->lumped, &this->lumped, NULL));
  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));
  PetscCall(PetscOptionsBool("-material_fields", "Material properties are assigned to auxiliary fields", NULL, this->material_fields, &this->material_fields, NULL));

  PetscCall(PetscOptionsBool("-gravity", "Use gravity in the y direction in 2D and in the z direction in 3D", NULL, this->transientdarcysflow.useGravity, &this->transientdarcysflow.useGravity, NULL));

  derived_sol_elas_type = this->isotropiclinearelasticity.derivedsoltype;
  
  PetscCall(PetscOptionsEList("-derived_solution_elasicity_type", "Type of derived solutions", NULL, this->isotropiclinearelasticity.DerivedSolutionTypes, isotropiclinearelasticity.NUM_DERIVED_SOLUTION_TYPES, this->isotropiclinearelasticity.DerivedSolutionTypes[this->isotropiclinearelasticity.derivedsoltype], &derived_sol_elas_type, &flg_elas));
  if (flg_elas)
  {
    this->isotropiclinearelasticity.derivedsoltype = (Elasticity::DerivedSolutionType)derived_sol_elas_type;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Type of derived solutions for the elasticity process is not given. " <<  isotropiclinearelasticity.DerivedSolutionTypes[0] << " is assumed." << std::endl;
    }
  }

  derived_sol_heat_type = this->transientthermal.derivedsoltype;

  PetscCall(PetscOptionsEList("-derived_solution_thermal_type", "Type of derived solutions", NULL, this->transientthermal.DerivedSolutionTypes, transientthermal.NUM_DERIVED_SOLUTION_TYPES, this->transientthermal.DerivedSolutionTypes[this->transientthermal.derivedsoltype], &derived_sol_heat_type, &flg_heat));
  if (flg_heat)
  {
    this->transientthermal.derivedsoltype = (Thermal::DerivedSolutionType)derived_sol_heat_type;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Type of derived solution for the heat transfer process is not given. " << transientthermal.DerivedSolutionTypes[0] << " is assumed." << std::endl;}
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);

}

PetscErrorCode TransientThermoPoroElasticity::SetupTS(DM dm, Vec *sol, const char name[], TS* ts) 
{
    char           view_call[PETSC_MAX_PATH_LEN]; 

    PetscFunctionBegin;

    if (this->expTS) {
      std::cout << "hahaha" << std::endl;
        PetscCall(DMTSSetRHSFunctionLocal(dm, DMPlexTSComputeRHSFunctionFEM, NULL));
        if (this->lumped) PetscCall(DMTSCreateRHSMassMatrixLumped(dm));
        else            PetscCall(DMTSCreateRHSMassMatrix(dm));
    } else {
      std::cout << "hoha" << std::endl;
        PetscCall(DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, NULL));
        PetscCall(DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, NULL));
    }
    PetscCall(TSSetExactFinalTime(*ts, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(*ts));
    PetscCall(TSSetComputeInitialCondition(*ts, SetInitialConditions));


    PetscCall(DMCreateGlobalVector(dm, sol));

    PetscCall(DMTSCheckFromOptions(*ts, *sol));
    PetscCall(SetInitialConditions(*ts, *sol));

    PetscCall(PetscObjectSetName((PetscObject) *sol, name));

    PetscFunctionReturn(0); 
}



PetscErrorCode TransientThermoPoroElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscFunctionBegin;   

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientThermoPoroElasticity::SetMaterial(MPI_Comm comm) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);
    PetscScalar ther_cond_temp(3.5);
    PetscScalar fluid_density_temp(1000.0), fluid_heat_capa_temp(4182.0); // flow_vel_temp(0.0);
    PetscScalar lin_ther_expan_c_temp(0.0);
    PetscScalar density_temp(2700.0), heat_capa_temp(800.0), spec_heat_temp(1.0);
    PetscScalar permeability_darcy(1e-15); // /* m^2 */ // fluid_density_darcy(1000.0)/* kg/m^3 */
    PetscScalar fluid_viscosity_darcy(0.001)/* Pa s */;
    PetscScalar gravity_darcy(0.0) ;/* m/s^2 */
    // PetscScalar compressibility_darcy(5e-10) /* 1/Pa */, porosity_darcy(0.05) /*  */ ;
    PetscScalar biot_modulus_hm(1.0e10) /* Pa */, biot_coefficient_hm(1.0) ;

    char option_call[PETSC_MAX_PATH_LEN];
    char material_name[128];
    char material_description[128];
    PetscBool flg_material_name, flg_material_description;
    std::string name(""), description("");
    PetscInt label_number(0);

    PetscFunctionBegin;   

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-youngs_modulus"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &youngs_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-lames_first_parameter"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lambda_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-shear_modulus"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &shear_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-poissons_ratio"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &poisson_r_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-bulk_modulus"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &bulk_m_temp, NULL);

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-name"));
    PetscOptionsGetString(NULL, NULL, option_call, material_name, sizeof(material_name), &flg_material_name);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-description"));
    PetscOptionsGetString(NULL, NULL, option_call, material_description, sizeof(material_description), &flg_material_description);

    if (flg_material_name) name = material_name; 
    if (flg_material_description) description = material_description;

    switch (this->isotropiclinearelasticity.inputpara) {
    case IsotropicLinearElasticity::InputParameterType::YOUNGSMODULUS_POISSONSRATIO:
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      shear_m_temp = youngs_m_temp / (2.0*(1.0+poisson_r_temp));
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      break;
    case IsotropicLinearElasticity::InputParameterType::YOUNGSMODULUS_SHEARMODULUS:
      poisson_r_temp = youngs_m_temp / (2.0*shear_m_temp) -1.0 ;
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      break;
    case IsotropicLinearElasticity::InputParameterType::BULKMODULUS_SHAERMODULUS:
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      lambda_temp = bulk_m_temp - (2.0/3.0*shear_m_temp);
      break;
    case IsotropicLinearElasticity::InputParameterType::LAMESPARAMETERS:
      bulk_m_temp = lambda_temp + (2.0 / 3.0 * shear_m_temp);
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid input parameter type: %s", isotropiclinearelasticity.inputParameterTypes[PetscMin(this->isotropiclinearelasticity.inputpara, isotropiclinearelasticity.NUM_INPUT_PARA_TYPES)]);
    }

    
    Property youngs_m(youngs_m_temp, "Pa", "Young's modulus, E" ); /* Young's modulus, E */
    Property shear_m(shear_m_temp, "Pa", "Shear modulus, G (Lame's second parameter, mu)" ); /* G = Lame's second parameter, mu */
    Property lambda(lambda_temp, "Pa", "Lame's first parameter, lambda" ); /*Lame's first parameter, lambda*/
    Property poisson_r(poisson_r_temp, " ", "Poisson's ratio, nu" ); /* Poisson's ratio, nu  */
    Property bulk_m(bulk_m_temp, "Pa", "bulk modulus, K" ); /* bulk modulus, K */
    
    material.SetName(name);
    material.SetLabel_number(label_number); 
    material.SetDescription(description); 
    material.AddProperty(shear_m); // 0
    material.AddProperty(lambda); // 1
    material.AddProperty(poisson_r); // 2
    material.AddProperty(youngs_m); // 3
    material.AddProperty(bulk_m); // 4


    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-thermal_conductivity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &ther_cond_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-density"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-heat_capacity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-specific_heat"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &spec_heat_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_density"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_heat_capacity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_heat_capa_temp, NULL);

    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" );  
    Property density(density_temp, "kg/m3", "Density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property spec_heat(spec_heat_temp, " ", "Specific Heat" );
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" );  
    Property fluid_heat_capa(fluid_heat_capa_temp, "J/(kg K)", "Fluid heat capacity" );  

    material.AddProperty(ther_cond); // 5
    material.AddProperty(fluid_density); // 6
    material.AddProperty(fluid_heat_capa); // 7

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-lin_ther_expan_c"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lin_ther_expan_c_temp, NULL);

    Property lin_ther_expan_c(lin_ther_expan_c_temp, "1/K", "Linear thermal expansion coefficient" ); 
    material.AddProperty(lin_ther_expan_c); // 8

    material.AddProperty(density); // 9
    material.AddProperty(heat_capa); // 10
    material.AddProperty(spec_heat); // 11

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-permeability"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &permeability_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_viscosity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_viscosity_darcy, NULL);
    if (transientdarcysflow.useGravity == PETSC_TRUE) {gravity_darcy = 9.8;}

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-biot_modulus"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_modulus_hm, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-biot_coefficient"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_coefficient_hm, NULL);
    
    Property permeability(permeability_darcy, "m2", "Permeability" ); 
    Property fluid_viscosity(fluid_viscosity_darcy, "Pa s", "Fluid viscosity" ); 
    Property gravity(gravity_darcy, "m/s^2", "Gravity" ); 
    Property biot_modulus(biot_modulus_hm, "1/Pa", "Biot Modulus" ); 
    Property biot_coefficient(biot_coefficient_hm, "-", "Biot coefficient" ); 
       

    material.AddProperty(permeability); //12
    material.AddProperty(fluid_viscosity); // 13
    material.AddProperty(gravity); // 14
    material.AddProperty(biot_modulus);  //15
    material.AddProperty(biot_coefficient); //16

    PetscFunctionReturn(0); 

}


PetscErrorCode TransientThermoPoroElasticity::SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);

    char option_call[PETSC_MAX_PATH_LEN];
    char material_name[128];
    char material_description[128];
    PetscBool flg_material_name, flg_material_description; 
    std::string name("");
    Material mtrl;

    PetscFunctionBegin;  

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

    mtrl.SetName(name);
    mtrl.SetLabel_number(label_number); 
    mtrl.SetDescription(description); 
 

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-youngs_modulus_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &youngs_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-lames_first_parameter_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lambda_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-shear_modulus_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &shear_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-poissons_ratio_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &poisson_r_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-bulk_modulus_%d", label_number));
    PetscOptionsGetScalar(NULL, NULL, option_call, &bulk_m_temp, NULL);


    switch (this->isotropiclinearelasticity.inputpara) {
    case IsotropicLinearElasticity::InputParameterType::YOUNGSMODULUS_POISSONSRATIO:
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      shear_m_temp = youngs_m_temp / (2.0*(1.0+poisson_r_temp));
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      break;
    case IsotropicLinearElasticity::InputParameterType::YOUNGSMODULUS_SHEARMODULUS:
      poisson_r_temp = youngs_m_temp / (2.0*shear_m_temp) -1.0 ;
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      break;
    case IsotropicLinearElasticity::InputParameterType::BULKMODULUS_SHAERMODULUS:
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      lambda_temp = bulk_m_temp - (2.0/3.0*shear_m_temp);
      break;
    case IsotropicLinearElasticity::InputParameterType::LAMESPARAMETERS:
      bulk_m_temp = lambda_temp + (2.0 / 3.0 * shear_m_temp);
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid input parameter type: %s", isotropiclinearelasticity.inputParameterTypes[PetscMin(this->isotropiclinearelasticity.inputpara, isotropiclinearelasticity.NUM_INPUT_PARA_TYPES)]);
    }

    Property youngs_m(youngs_m_temp, "Pa", "Young's modulus, E" ); /* Young's modulus, E */
    Property shear_m(shear_m_temp, "Pa", "Shear modulus, G (Lame's second parameter, mu)" ); /* G = Lame's second parameter, mu */
    Property lambda(lambda_temp, "Pa", "Lame's first parameter, lambda" ); /*Lame's first parameter, lambda*/
    Property poisson_r(poisson_r_temp, " ", "Poisson's ratio, nu" ); /* Poisson's ratio, nu  */
    Property bulk_m(bulk_m_temp, "Pa", "bulk modulus, K" ); /* bulk modulus, K */
    
    mtrl.AddProperty(shear_m); // 0
    mtrl.AddProperty(lambda); // 1
    mtrl.AddProperty(poisson_r); // 2
    mtrl.AddProperty(youngs_m); // 3
    mtrl.AddProperty(bulk_m); // 4

    Property ther_cond;
    Property fluid_density;
    Property fluid_heat_capa;
    // Property flow_vel;
    Property lin_ther_expan_c;
    Property density;
    Property heat_capa;
    Property spec_heat;
    Property heat_source; // volumetric (in 3D) or areal (in 2D) haet source
    Property ref_temp; 

    PetscBool heatsource_pwr_bool(PETSC_FALSE); //for high-level nuclear waste disposal modeling

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-heat_source_pwr_%d", label_number));
    PetscOptionsGetBool(NULL, NULL, option_call, &heatsource_pwr_bool, NULL);

    PetscCall(ther_cond.SetPropertyUser(comm, "W/(m K)", "Thermal conductivity", "thermal_conductivity", label_number));
    PetscCall(fluid_density.SetPropertyUser(comm, "kg/m3", "Fluid density", "fluid_density", label_number));
    PetscCall(fluid_heat_capa.SetPropertyUser(comm, "J/(kg K)", "Fluid heat capacity", "fluid_heat_capacity", label_number));
    // PetscCall(flow_vel.SetPropertyUser(comm, "m/s", "Flow velocity", "flow_velocity", label_number));
    PetscCall(lin_ther_expan_c.SetPropertyUser(comm, "1/K", "Linear thermal expansion coefficient", "lin_ther_expan_c", label_number));
    PetscCall(density.SetPropertyUser(comm, "kg/m3", "Density", "density", label_number));
    PetscCall(heat_capa.SetPropertyUser(comm, "J/(kg K)", "Heat capacity", "heat_capacity", label_number));
    PetscCall(spec_heat.SetPropertyUser(comm, " ", "Specific Heat", "specific_heat", label_number));
    if (heatsource_pwr_bool) PetscCall(heat_source.SetPropertyUser(comm, "W/m3", "Heat Source PWR", "heat_source_pwr_num_canisters_per_volume", label_number));
    else PetscCall(heat_source.SetPropertyUser(comm, "W/m3", "Heat Source", "heat_source", label_number));
    PetscCall(ref_temp.SetPropertyUser(comm, "degC", "Reference temperature for thermal stress analysis", "ref_temp", label_number));
    
    Property permeability; 
    Property fluid_viscosity; 
    Property gravity; 
    Property biot_modulus; 
    Property biot_coefficient; 
    Property ref_pres; 
    Property porosity;
    Property vol_flu_ther_expan_c;  
    
    PetscCall(permeability.SetPropertyUser(comm, "m2", "Permeability", "permeability", label_number));
    PetscCall(fluid_viscosity.SetPropertyUser(comm, "Pa s", "Fluid viscosity", "fluid_viscosity", label_number));
    PetscCall(gravity.SetPropertyUser(comm, "m/s^2", "Gravity", "gravity", label_number));
    PetscCall(biot_modulus.SetPropertyUser(comm, "1/Pa", "Biot Modulus", "biot_modulus", label_number));
    PetscCall(biot_coefficient.SetPropertyUser(comm, "-", "Biot coefficient", "biot_coefficient", label_number));
    PetscCall(ref_pres.SetPropertyUser(comm, "Pa", "Reference pressure for HM analysis", "ref_pres", label_number));
    PetscCall(porosity.SetPropertyUser(comm, "-", "Porosity for THM analysis", "porosity", label_number));
    PetscCall(vol_flu_ther_expan_c.SetPropertyUser(comm, "1/K", "Volumetric thermal expansion coefficient of fluid", "vol_flu_ther_expan_c", label_number));
        
    mtrl.AddProperty(ther_cond); // 5
    mtrl.AddProperty(fluid_density); // 6
    mtrl.AddProperty(fluid_heat_capa); // 7
    mtrl.AddProperty(lin_ther_expan_c); // 8
    mtrl.AddProperty(density); // 9
    mtrl.AddProperty(heat_capa); // 10
    mtrl.AddProperty(spec_heat); // 11
    mtrl.AddProperty(permeability); //12
    // mtrl.AddProperty(fluid_density);
    mtrl.AddProperty(fluid_viscosity); // 13
    mtrl.AddProperty(gravity); // 14
    mtrl.AddProperty(biot_modulus);  //15
    mtrl.AddProperty(biot_coefficient); //16
    mtrl.AddProperty(heat_source); // 17
    mtrl.AddProperty(ref_temp); // 18
    mtrl.AddProperty(ref_pres); // 19
    mtrl.AddProperty(porosity); // 20
    mtrl.AddProperty(vol_flu_ther_expan_c); // 21

    AddMaterial(mtrl);

    PetscFunctionReturn(0); 

}


PetscErrorCode TransientThermoPoroElasticity::SetMaterials(MPI_Comm comm) {
  
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

PetscErrorCode TransientThermoPoroElasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
{
    PetscFunctionBeginUser;
   
    isotropiclinearelasticity.BCSetup(dm, ds, wf, 0);

    if (material_fields) {transientthermal.BCSetup_TM_mat_fields(dm, ds, wf, 1);}
    else transientthermal.BCSetup_TM_mat_const(dm, ds, wf, 1);

    transientdarcysflow.BCSetup_HM(dm, ds, wf, 2);

    PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoPoroElasticity::DomainSetup(PetscDS ds) {
 
  PetscFunctionBeginUser;
  if (!material_fields) {
    PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_thermoporoelas_u));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_thermoporoelas_uu));
    PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_thermoporoelas_ut, NULL)); // 1
    PetscCall(PetscDSSetJacobian(ds, 0, 2, NULL, NULL, Pw_Functions::g2_thermoporoelas_up, NULL)); 
    PetscCall(PetscDSSetJacobian(ds, 0, 3, NULL, NULL, Pw_Functions::g2_thermoporoelas_ustrain, NULL)); 

    PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_thermoporoelas_t_transient, Pw_Functions::f1_thermoporoelas_t_transient));
    PetscCall(PetscDSSetJacobian(ds, 1, 1, Pw_Functions::g0_thermoporoelas_tt_transient, NULL, NULL, Pw_Functions::g3_thermoporoelas_tt));
    if (transientthermal.flg_iniv_equation) {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature_user, (void*)(&transientthermal.init_temperature_equation)));}
    else {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature, (void*)(&transientthermal.init_temperature)));}
    PetscCall(PetscDSSetExactSolutionTimeDerivative(ds, 1, Pw_Functions::initial_temperature_t, NULL));

    PetscCall(PetscDSSetResidual(ds, 2, Pw_Functions::f0_thermoporoelas_pressure_transient, Pw_Functions::f1_thermoporoelas_pressure));
    PetscCall(PetscDSSetJacobian(ds, 2, 1, Pw_Functions::g0_thermoporoelas_pressure_temperature_transient, NULL, NULL, NULL));
    PetscCall(PetscDSSetJacobian(ds, 2, 2, Pw_Functions::g0_thermoporoelas_pressure_pressure_transient, NULL, NULL, Pw_Functions::g3_thermoporoelas_pressure_pressure));
    PetscCall(PetscDSSetJacobian(ds, 2, 3, Pw_Functions::g0_thermoporoelas_pressure_strain_transient, NULL, NULL, NULL));
    PetscCall(PetscDSSetExactSolution(ds, 2, Pw_Functions::initial_pressure_thermoporoelas, (void*)&transientdarcysflow.init_pressure));
    PetscCall(PetscDSSetResidual(ds, 3, Pw_Functions::f0_thermoporoelas_strain, NULL));
    PetscCall(PetscDSSetJacobian(ds, 3, 0, NULL, Pw_Functions::g1_thermoporoelas_strain_u, NULL, NULL));
    PetscCall(PetscDSSetJacobian(ds, 3, 3, Pw_Functions::g0_thermoporoelas_strain_strain, NULL, NULL, NULL));
  }
  else {
    PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_thermoporoelas_u_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_thermoporoelas_uu_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_thermoporoelas_ut_mat_fields, NULL)); // 1
    PetscCall(PetscDSSetJacobian(ds, 0, 2, NULL, NULL, Pw_Functions::g2_thermoporoelas_up_mat_fields, NULL)); 
    PetscCall(PetscDSSetJacobian(ds, 0, 3, NULL, NULL, Pw_Functions::g2_thermoporoelas_ustrain_mat_fields, NULL)); 
    PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_thermoporoelas_t_transient_mat_fields, Pw_Functions::f1_thermoporoelas_t_transient_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 1, 1, Pw_Functions::g0_thermoporoelas_tt_transient_mat_fields, NULL, NULL, Pw_Functions::g3_thermoporoelas_tt_mat_fields));
    
    if (transientthermal.flg_iniv_equation) {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature_user, (void*)(&transientthermal.init_temperature_equation)));}
    else {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature, (void*)(&transientthermal.init_temperature)));}
    PetscCall(PetscDSSetExactSolutionTimeDerivative(ds, 1, Pw_Functions::initial_temperature_t, NULL));

    PetscCall(PetscDSSetResidual(ds, 2, Pw_Functions::f0_thermoporoelas_pressure_transient_mat_fields, Pw_Functions::f1_thermoporoelas_pressure_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 2, 1, Pw_Functions::g0_thermoporoelas_pressure_temperature_transient_mat_fields, NULL, NULL, NULL));
    PetscCall(PetscDSSetJacobian(ds, 2, 2, Pw_Functions::g0_thermoporoelas_pressure_pressure_transient_mat_fields, NULL, NULL, Pw_Functions::g3_thermoporoelas_pressure_pressure_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 2, 3, Pw_Functions::g0_thermoporoelas_pressure_strain_transient_mat_fields, NULL, NULL, NULL));
    PetscCall(PetscDSSetExactSolution(ds, 2, Pw_Functions::initial_pressure_thermoporoelas, (void*)&transientdarcysflow.init_pressure));
    
    PetscCall(PetscDSSetResidual(ds, 3, Pw_Functions::f0_thermoporoelas_strain, NULL));
    PetscCall(PetscDSSetJacobian(ds, 3, 0, NULL, Pw_Functions::g1_thermoporoelas_strain_u, NULL, NULL));
    PetscCall(PetscDSSetJacobian(ds, 3, 3, Pw_Functions::g0_thermoporoelas_strain_strain, NULL, NULL, NULL));
}
  PetscFunctionReturn(0);
}


PetscErrorCode TransientThermoPoroElasticity::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));

  PetscFunctionReturn(0);
}


PetscErrorCode TransientThermoPoroElasticity::SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe)
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

  /* displacement */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[0]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, dim, simplex, name[0] ? prefix : NULL, -1, &fe[0]));
  PetscCall(PetscObjectSetName((PetscObject) fe[0], name[0]));
 
 
  /* temperature */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[1]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[1] ? prefix : NULL, -1, &fe[1]));
  PetscCall(PetscFECopyQuadrature(fe[0], fe[1]));
  PetscCall(PetscObjectSetName((PetscObject) fe[1], name[1]));

  /* pressure */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[2]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[2] ? prefix : NULL, -1, &fe[2]));
  PetscCall(PetscFECopyQuadrature(fe[1], fe[2]));
  PetscCall(PetscObjectSetName((PetscObject) fe[2], name[2]));

  /* volumetric strain */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[3]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[3] ? prefix : NULL, -1, &fe[3]));
  PetscCall(PetscFECopyQuadrature(fe[2], fe[3]));
  PetscCall(PetscObjectSetName((PetscObject) fe[3], name[3]));

   /* Create finite element for material properties */
  if (this->material_fields) {
    properties_num = this->materials[0].GetNumofProperties();
  }
  PetscFE  feAux[properties_num];

//option bool을 넣고, true 이면 아래 수정해서 property 수만큼 수행. q
  if (this->material_fields) {
    for (int i = 0; i < properties_num ; i++) {
      PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, "Material_", PETSC_DEFAULT, &feAux[i]));
      PetscCall(PetscFECopyQuadrature(fe[0], feAux[i]));
      PetscCall(PetscObjectSetName((PetscObject) feAux[i], const_cast<char*>(this->materials[0].GetPropertyDescription(i).c_str())));
    } 
  }
  
  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe[0]));
  PetscCall(DMSetField(dm, 1, NULL, (PetscObject) fe[1]));
  PetscCall(DMSetField(dm, 2, NULL, (PetscObject) fe[2]));
  PetscCall(DMSetField(dm, 3, NULL, (PetscObject) fe[3]));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));

  if (this->expTS) {
    for (PetscInt f = 0; f < 4; ++f) {PetscCall(PetscDSSetImplicit(ds, f, PETSC_FALSE));}
  }

  PetscCall(SetupPrimalProblem(dm, ds));
  

  if (!this->material_fields) {
    PetscCall(ParameterSetup(ds));
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if ((this->material_view) & (rank == 0)) PetscCall(material.ViewMaterial());
    while (cdm)
    {
      PetscCall(DMCopyDisc(dm, cdm));
      if (this->isotropiclinearelasticity.useNearNullspace)
        PetscCall(DMSetNearNullSpaceConstructor(cdm, 0, CreateElasticityNullSpace));
      /* TODO: Check whether the boundary of coarse meshes is marked */
      PetscCall(DMGetCoarseDM(cdm, &cdm));
    }
  }
  else {
    while (cdm)
    {
      PetscCall(SetupAuxDM(cdm, properties_num, feAux, 0.0, 0));
      PetscCall(DMCopyDisc(dm, cdm));
      if (this->isotropiclinearelasticity.useNearNullspace)
        PetscCall(DMSetNearNullSpaceConstructor(cdm, 0, CreateElasticityNullSpace));
      /* TODO: Check whether the boundary of coarse meshes is marked */
      PetscCall(DMGetCoarseDM(cdm, &cdm));
    }
  }
  
  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFEDestroy(&feAux[i]));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoPoroElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


PetscErrorCode TransientThermoPoroElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
{
  PetscInt       para = 0; 
  std::string option_name = "-";
  PetscBool      flg;
  PetscMPIInt    rank;

  PetscFunctionBeginUser;
  option_name.append(name);
  option_name.append("_inputpara_type");
  // std::cout << option_name << std::endl;
  PetscOptionsBegin(comm, "", "Type of two input elastic parameters", "");
  para = this->isotropiclinearelasticity.inputpara;

  PetscCall(PetscOptionsEList(option_name.c_str(), "Type of two input elastic parameters", NULL, this->isotropiclinearelasticity.inputParameterTypes, isotropiclinearelasticity.NUM_INPUT_PARA_TYPES, this->isotropiclinearelasticity.inputParameterTypes[this->isotropiclinearelasticity.inputpara], &para, &flg));
  if (flg)
  {
    this->isotropiclinearelasticity.inputpara = (IsotropicLinearElasticity::InputParameterType)para;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Input type for two elastic parameters is not given. LAMESPARAMETERS are assumed to be given." << std::endl;
    }
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

TransientThermoPoroElasticity::TransientThermoPoroElasticity()
{
  this->expTS = PETSC_FALSE;
  this->lumped = PETSC_FALSE;

}
