#include "TransientThermoElasticity.hh"


PetscErrorCode TransientThermoElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) {
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

PetscErrorCode TransientThermoElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step, PetscFE feAux[]) {
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

PetscErrorCode TransientThermoElasticity::SolutionMonitor_TM(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
  PetscBool      simplex;
  PetscInt       dim, comp;
  PetscDS        ds;
  PetscFE        fe[2];  //, fe_derived;
  DM             dm;//, dmAux;

  TransientThermoElasticity* tte = (TransientThermoElasticity*)ctx; 
  PetscInt properties_num(0);

/* displacement, temperature, heat flux, total stress, effective stress, strain  */

  PetscFunctionBeginUser;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMGetDS(dm, &ds));

  PetscCall(DMGetField(dm, 0, NULL, (PetscObject *)&fe[0])); /* displacement */
  PetscCall(DMGetField(dm, 1, NULL, (PetscObject *)&fe[1])); /* temperature */

  
  if (tte->material_fields) {
  properties_num = tte->materials[0].GetNumofProperties();
  }
  PetscFE  feAux[properties_num];
  PetscCall(DMGetDimension(dm, &dim));
  DMPlexIsSimplex(dm, &simplex);

  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, "Material_", PETSC_DEFAULT, &feAux[i]));
    PetscCall(PetscFECopyQuadrature(fe[0], feAux[i]));
    PetscCall(PetscObjectSetName((PetscObject) feAux[i], const_cast<char*>(tte->materials[0].GetPropertyDescription(i).c_str())));
  } 
    PetscCall(VecViewFromOptions(u, NULL, "-vec_sol_view_2"));  
    PetscCall(tte->SetupAuxDM(dm, properties_num, feAux, time, step, &u));
  
    PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::displacement_TM, &u, "displacement", time, step, dim, dim));
    PetscCall(CreateSolveDerivedField(dm, &fe[1], Pw_Functions::temperature_TM, &u, "temperature", time, step, 1, dim));
    if (tte->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::heatflux_TM_mat_fields, &u, "heat_flux", time, step, dim, dim));
    else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::heatflux_TM, &u, "heat_flux", time, step, dim, dim));

    if (dim == 2) {comp = 4;}
    else if (dim == 3) {comp = 6;}

    switch (tte->isotropiclinearelasticity.derivedsoltype) {
      case Elasticity::DerivedSolutionType::PLANESTRAIN_2D:
      if (tte->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_planestrain_TM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_planestrain_TM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystrain_planestrain_TM, &u, "strain", time, step, comp, dim));
      break;
      case Elasticity::DerivedSolutionType::CAUCHY_3D:
      if (tte->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_3D_TM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystress_3D_TM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::cauchystrain_3D_TM, &u, "strain", time, step, comp, dim));
      break;
      case Elasticity::DerivedSolutionType::AXISYMMETRIC_2D_PLANESTRAIN:
      if (tte->material_fields) PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_stress_planestrain_TM_mat_fields, &u, "stress", time, step, comp, dim));
      else PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_stress_planestrain_TM, &u, "stress", time, step, comp, dim));
      PetscCall(CreateSolveDerivedField(dm, &fe[0], Pw_Functions::axisymmetric_2d_strain_planestrain_TM, &u, "strain", time, step, comp, dim));
      break;
      default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type for elasticity");
    }   
 
  for (int i = 0; i < properties_num ; i++) {
    PetscCall(PetscFEDestroy(&feAux[i]));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoElasticity::SolveFE(DM dm, Vec *sol, const char name[], TS* ts)
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;

  PetscCall(TSSolve(*ts, *sol));
  PetscCall(DMTSCheckFromOptions(*ts, *sol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*sol, NULL, view_call));

  PetscFunctionReturn(0); 
}


PetscErrorCode TransientThermoElasticity::ProcessOptions(MPI_Comm comm) {
    
  PetscInt       derived_sol_elas_type = 0, derived_sol_heat_type = 0, rank; 
  PetscBool      flg_elas, flg_heat;
  
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-near_nullspace", "Use the rigid body modes as an AMG near nullspace", NULL, this->isotropiclinearelasticity.useNearNullspace, &this->isotropiclinearelasticity.useNearNullspace, NULL)); 

  PetscCall(PetscOptionsBool("-explicit", "Use explicit timestepping", NULL, this->transientthermal.expTS, &this->transientthermal.expTS, NULL));
  PetscCall(PetscOptionsBool("-lumped", "Lump the mass matrix", NULL, this->transientthermal.lumped, &this->transientthermal.lumped, NULL));

  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));
  PetscCall(PetscOptionsBool("-material_fields", "Material properties are assigned to auxiliary fields", NULL, this->material_fields, &this->material_fields, NULL));

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

PetscErrorCode TransientThermoElasticity::SetupTS(DM dm, Vec *sol, const char name[], TS* ts) 
{
    char           view_call[PETSC_MAX_PATH_LEN]; 

    PetscFunctionBegin;

    if (this->transientthermal.expTS) {
        PetscCall(DMTSSetRHSFunctionLocal(dm, DMPlexTSComputeRHSFunctionFEM, NULL));
        if (this->transientthermal.lumped) PetscCall(DMTSCreateRHSMassMatrixLumped(dm));
        else            PetscCall(DMTSCreateRHSMassMatrix(dm));
    } else {
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


// not used
PetscErrorCode TransientThermoElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);
    PetscScalar ther_cond_temp(3.5);
    PetscScalar fluid_density_temp(1000.0), fluid_heat_capa_temp(4182.0), flow_vel_temp(0.0);
    PetscScalar lin_ther_expan_c_temp(0.0);
    PetscScalar density_temp(2700.0), heat_capa_temp(800.0), spec_heat_temp(1.0);
    
    char option_call[PETSC_MAX_PATH_LEN];

    PetscFunctionBegin;   

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_youngs_modulus", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &youngs_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_lames_first_parameter", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lambda_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_shear_modulus", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &shear_m_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_poissons_ratio", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &poisson_r_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_bulk_modulus", const_cast<char*>(name.c_str())));
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
    
    material.SetName(name);
    material.SetLabel_number(label_number); 
    material.SetDescription(description); 
    material.AddProperty(shear_m); // 0
    material.AddProperty(lambda); // 1
    material.AddProperty(poisson_r); // 2
    material.AddProperty(youngs_m); // 3
    material.AddProperty(bulk_m); // 4


    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_thermal_conductivity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &ther_cond_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_density", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_heat_capacity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_specific_heat", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &spec_heat_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_density", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_heat_capacity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_flow_velocity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);

    
    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" );  
    Property density(density_temp, "kg/m3", "Density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property spec_heat(spec_heat_temp, " ", "Specific Heat" );
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" );  
    Property fluid_heat_capa(fluid_heat_capa_temp, "J/(kg K)", "Fluid heat capacity" );  
    Property flow_vel(flow_vel_temp, "m/s", "Flow velocity" );  

    material.AddProperty(ther_cond); // 5
    material.AddProperty(fluid_density); // 6
    material.AddProperty(fluid_heat_capa); // 7
    material.AddProperty(flow_vel); // 8

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_lin_ther_expan_c", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lin_ther_expan_c_temp, NULL);

    Property lin_ther_expan_c(lin_ther_expan_c_temp, "1/K", "Linear thermal expansion coefficient" ); 
    material.AddProperty(lin_ther_expan_c); // 9

    material.AddProperty(density); // 10
    material.AddProperty(heat_capa); // 11
    material.AddProperty(spec_heat); // 12

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientThermoElasticity::SetMaterial(MPI_Comm comm) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);
    PetscScalar ther_cond_temp(3.5);
    PetscScalar fluid_density_temp(1000.0), fluid_heat_capa_temp(4182.0), flow_vel_temp(0.0);
    PetscScalar lin_ther_expan_c_temp(0.0);
    PetscScalar density_temp(2700.0), heat_capa_temp(800.0), spec_heat_temp(1.0);
    
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
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-flow_velocity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);

    
    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" );  
    Property density(density_temp, "kg/m3", "Density" ); 
    Property heat_capa(heat_capa_temp, "J/(kg K)", "Heat capacity" ); 
    Property spec_heat(spec_heat_temp, " ", "Specific Heat" );
    Property fluid_density(fluid_density_temp, "kg/m3", "Fluid density" );  
    Property fluid_heat_capa(fluid_heat_capa_temp, "J/(kg K)", "Fluid heat capacity" );  
    Property flow_vel(flow_vel_temp, "m/s", "Flow velocity" );  
     
    material.AddProperty(ther_cond); // 5
    material.AddProperty(fluid_density); // 6
    material.AddProperty(fluid_heat_capa); // 7
    material.AddProperty(flow_vel); // 8

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-lin_ther_expan_c"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &lin_ther_expan_c_temp, NULL);

    Property lin_ther_expan_c(lin_ther_expan_c_temp, "1/K", "Linear thermal expansion coefficient" ); 
    material.AddProperty(lin_ther_expan_c); // 9

    material.AddProperty(density); // 10
    material.AddProperty(heat_capa); // 11
    material.AddProperty(spec_heat); // 12

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientThermoElasticity::SetMaterial(MPI_Comm comm, PetscInt label_number, std::string description) {

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
    Property flow_vel;
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
    PetscCall(flow_vel.SetPropertyUser(comm, "m/s", "Flow velocity", "flow_velocity", label_number));
    PetscCall(lin_ther_expan_c.SetPropertyUser(comm, "1/K", "Linear thermal expansion coefficient", "lin_ther_expan_c", label_number));
    PetscCall(density.SetPropertyUser(comm, "kg/m3", "Density", "density", label_number));
    PetscCall(heat_capa.SetPropertyUser(comm, "J/(kg K)", "Heat capacity", "heat_capacity", label_number));
    PetscCall(spec_heat.SetPropertyUser(comm, " ", "Specific Heat", "specific_heat", label_number));
    if (heatsource_pwr_bool) PetscCall(heat_source.SetPropertyUser(comm, "W/m3", "Heat Source PWR", "heat_source_pwr_num_canisters_per_volume", label_number));
    else PetscCall(heat_source.SetPropertyUser(comm, "W/m3", "Heat Source", "heat_source", label_number));
    PetscCall(ref_temp.SetPropertyUser(comm, "degC", "Reference temperature for thermal stress analysis", "ref_temp", label_number));
    
    mtrl.AddProperty(ther_cond); // 5
    mtrl.AddProperty(fluid_density); // 6
    mtrl.AddProperty(fluid_heat_capa); // 7
    mtrl.AddProperty(flow_vel); // 8

    mtrl.AddProperty(lin_ther_expan_c); // 9
    mtrl.AddProperty(density); // 10
    mtrl.AddProperty(heat_capa); // 11
    mtrl.AddProperty(spec_heat); // 12
    mtrl.AddProperty(heat_source); // 13
    mtrl.AddProperty(ref_temp); // 14

    AddMaterial(mtrl);

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientThermoElasticity::SetMaterials(MPI_Comm comm) {
  
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

PetscErrorCode TransientThermoElasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
{
    PetscFunctionBeginUser;
   
    isotropiclinearelasticity.BCSetup(dm, ds, wf, 0);

    if (material_fields) {transientthermal.BCSetup_TM_mat_fields(dm, ds, wf, 1);}
    else transientthermal.BCSetup_TM_mat_const(dm, ds, wf, 1);
    

    PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoElasticity::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  if (!material_fields) {
    PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_thermoelas_u));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_thermoelas_uu));
    PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_thermoelas_ut, NULL)); // 1
    PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_thermoelas_t_transient, Pw_Functions::f1_thermoelas_t_transient));
    PetscCall(PetscDSSetJacobian(ds, 1, 1, Pw_Functions::g0_thermoelas_tt_transient, Pw_Functions::g1_thermoelas_tt, NULL, Pw_Functions::g3_thermoelas_tt));
    if (transientthermal.flg_iniv_equation) {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature_user, (void*)(&transientthermal.init_temperature_equation)));}
    else {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature, (void*)(&transientthermal.init_temperature)));}
    PetscCall(PetscDSSetExactSolutionTimeDerivative(ds, 1, Pw_Functions::initial_temperature_t, NULL));
  }
  else {
    PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_thermoelas_u_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_thermoelas_uu_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_thermoelas_ut_mat_fields, NULL)); // 1
    PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_thermoelas_t_transient_mat_fields, Pw_Functions::f1_thermoelas_t_transient_mat_fields));
    PetscCall(PetscDSSetJacobian(ds, 1, 1, Pw_Functions::g0_thermoelas_tt_transient_mat_fields, Pw_Functions::g1_thermoelas_tt_mat_fields, NULL, Pw_Functions::g3_thermoelas_tt_mat_fields));
    if (transientthermal.flg_iniv_equation) {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature_user, (void*)(&transientthermal.init_temperature_equation)));}
    else {PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_temperature, (void*)(&transientthermal.init_temperature)));}
    PetscCall(PetscDSSetExactSolutionTimeDerivative(ds, 1, Pw_Functions::initial_temperature_t, NULL));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TransientThermoElasticity::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));
  PetscFunctionReturn(0);
}


PetscErrorCode TransientThermoElasticity::SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe)
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
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));

  if (this->transientthermal.expTS) {
    for (PetscInt f = 0; f < 2; ++f) {PetscCall(PetscDSSetImplicit(ds, f, PETSC_FALSE));}
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
      // PetscCall(SetupAuxDM(dm, properties_num, feAux));
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

PetscErrorCode TransientThermoElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


PetscErrorCode TransientThermoElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
{
  PetscInt       para = 0; 
  std::string option_name = "-";
  PetscBool      flg;
  PetscMPIInt    rank;

  PetscFunctionBeginUser;
  option_name.append(name);
  option_name.append("_inputpara_type");
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

TransientThermoElasticity::TransientThermoElasticity()
{

}
