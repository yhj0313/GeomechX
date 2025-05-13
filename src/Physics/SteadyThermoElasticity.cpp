#include "SteadyThermoElasticity.hh"

PetscErrorCode SteadyThermoElasticity::SolveFE(DM dm, Vec *sol, const char name[], SNES* snes)
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
  std::cout << view_call << std::endl;
  PetscCall(VecViewFromOptions(*sol, NULL, view_call));

  PetscFunctionReturn(0);  
}

PetscErrorCode SteadyThermoElasticity::ProcessOptions(MPI_Comm comm)
{
  PetscInt       derived_sol_elas_type = 0, derived_sol_heat_type = 0, rank; 
  PetscBool      flg_elas, flg_heat;
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-near_nullspace", "Use the rigid body modes as an AMG near nullspace", NULL, this->isotropiclinearelasticity.useNearNullspace, &this->isotropiclinearelasticity.useNearNullspace, NULL));
  // PetscCall(PetscOptionsBool("-radial_tangential_2D", "Calculate radial and tangential stress and strain", NULL, this->radial_tangential_2D, &this->radial_tangential_2D, NULL));
  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));

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

  
  derived_sol_heat_type = this->steadythermal.derivedsoltype;

  PetscCall(PetscOptionsEList("-derived_solution_thermal_type", "Type of derived solutions", NULL, this->steadythermal.DerivedSolutionTypes, steadythermal.NUM_DERIVED_SOLUTION_TYPES, this->steadythermal.DerivedSolutionTypes[this->steadythermal.derivedsoltype], &derived_sol_heat_type, &flg_heat));
  if (flg_heat)
  {
    // this->steadythermal.derivedsoltype = (steadythermal.DerivedSolutionType)derived_sol_heat_type;
    this->steadythermal.derivedsoltype = (Thermal::DerivedSolutionType)derived_sol_heat_type;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Type of derived solution for the heat transfer process is not given. " << steadythermal.DerivedSolutionTypes[0] << " is assumed." << std::endl;}
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermoElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);
    PetscScalar ther_cond_temp(3.5), fluid_density_temp(1000.0), fluid_heat_capa_temp(4182.0), flow_vel_temp(0.0);
    PetscScalar lin_ther_expan_c_temp(0.0);
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
    Property lambda(lambda_temp, "MPa", "Lame's first parameter, lambda" ); /*Lame's first parameter, lambda*/
    Property poisson_r(poisson_r_temp, "MPa", "Poisson's ratio, nu" ); /* Poisson's ratio, nu  */
    Property bulk_m(bulk_m_temp, "MPa", "bulk modulus, K" ); /* bulk modulus, K */
    
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
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_density", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_heat_capacity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_heat_capa_temp, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_flow_velocity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &flow_vel_temp, NULL);

    
    Property ther_cond(ther_cond_temp, "W/(m K)", "Thermal conductivity" ); 
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

    PetscFunctionReturn(0); 

}

PetscErrorCode SteadyThermoElasticity::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_thermoelas_u));
  PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_thermoelas_uu));
  PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_thermoelas_ut, NULL)); // 1
  // PetscCall(PetscDSSetJacobian(ds, 1, 0, NULL, NULL, Pw_Functions::g2_thermoelas_ut, NULL)); // 2  // 3 no use of the two
  PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_thermoelas_t, Pw_Functions::f1_thermoelas_t));
  PetscCall(PetscDSSetJacobian(ds, 1, 1, NULL, Pw_Functions::g1_thermoelas_tt, NULL, Pw_Functions::g3_thermoelas_tt));
  PetscFunctionReturn(0);
}


PetscErrorCode SteadyThermoElasticity::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));

  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermoElasticity::SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe)
{
  DM             cdm  = dm;
  char           prefix[PETSC_MAX_PATH_LEN];
  PetscBool      simplex;
  PetscInt       dim;  
  PetscMPIInt    rank;

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
 std::cout << name[0] << std::endl;
 std::cout << name[1] << std::endl;
  
  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe[0]));
  PetscCall(DMSetField(dm, 1, NULL, (PetscObject) fe[1]));

  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(SetupPrimalProblem(dm, ds));
  PetscCall(ParameterSetup(ds));

  CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  if ((this->material_view) & (rank == 0)) PetscCall(material.ViewMaterial());
  
  while (cdm) {
    PetscCall(DMCopyDisc(dm, cdm));
    if (this->isotropiclinearelasticity.useNearNullspace) PetscCall(DMSetNearNullSpaceConstructor(cdm, 0, CreateElasticityNullSpace));
    /* TODO: Check whether the boundary of coarse meshes is marked */
    PetscCall(DMGetCoarseDM(cdm, &cdm));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SteadyThermoElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


PetscErrorCode SteadyThermoElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
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

PetscErrorCode SteadyThermoElasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
{
    PetscFunctionBeginUser;
   
    isotropiclinearelasticity.BCSetup(dm, ds, wf, 0);
    steadythermal.BCSetup(dm, ds, wf, 1);

    PetscFunctionReturn(0);
}
    
PetscErrorCode SteadyThermoElasticity::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) 
{
  PetscInt       n;
  
  PetscFunctionBeginUser;

  switch (isotropiclinearelasticity.derivedsoltype) {
    case Elasticity::DerivedSolutionType::PLANESTRAIN_2D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_planestrain, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_planestrain_TM, sol, &derivedsols[1], "stress");
    break;
    case Elasticity::DerivedSolutionType::CAUCHY_3D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_3D, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_3D_TM, sol, &derivedsols[1], "stress");
    break;
    case Elasticity::DerivedSolutionType::AXISYMMETRIC_2D_PLANESTRAIN:
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_strain_planestrain, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_stress_planestrain_TM, sol, &derivedsols[1], "stress");
    break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type");

  }

  PetscFunctionReturn(0);
   
}


