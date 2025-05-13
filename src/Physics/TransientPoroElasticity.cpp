#include "TransientPoroElasticity.hh"


PetscErrorCode TransientPoroElasticity::SolutionMonitor_HM(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
  PetscBool      simplex;
  PetscInt       dim, comp;
  PetscDS        ds;
  PetscFE        fe[2], fe_derived;
  DM             dm, dmAux;
  Vec            derived_solutions[6];
//   TransientDarcysFlow t_darcy ;
  TransientPoroElasticity* tpe = (TransientPoroElasticity*)ctx; 

/* displacement, darcy velocity, pressure, total stress, effective stress, strain  */

  PetscFunctionBeginUser;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(DMGetField(dm, 0, NULL, (PetscObject *)&fe[0])); /* displacement */
  // PetscCall(DMGetField(dm, 1, NULL, (PetscObject *)&fe[1])); /* darcy velocity */
  PetscCall(DMGetField(dm, 1, NULL, (PetscObject *)&fe[1])); /* pore pressure */

  PetscCall(DMGetDimension(dm, &dim));

  PetscCall(tpe->CreateDerivedField(dm, &dmAux, "displacement", &fe[0], &fe_derived, dim, dim)); /* displacement y(# of component = dim) */
  PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::scale_displacement_poroelas, &u, &derived_solutions[0], "displacement", time, step));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));

  PetscCall(tpe->CreateDerivedField(dm, &dmAux, "scaled_pressure", &fe[1], &fe_derived, 1, dim)); /* pressure(# of component = 1) */
  PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::scale_pressure_poroelas, &u, &derived_solutions[1], "scaled_pressure", time, step));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));
  PetscCall(tpe->CreateDerivedField(dm, &dmAux, "scaled_darcyvelocity", &fe[0], &fe_derived, dim, dim)); /* darcy velocity(# of component = dim) */
  PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::scale_darcyvelocity_poroelas, &u, &derived_solutions[2], "scaled_darcyvelocity", time, step));
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));
  if (dim == 2) {comp = 4;}
  else if (dim == 3) {comp = 6;}
  PetscCall(tpe->CreateDerivedField(dm, &dmAux, "derived_stress_strain", &fe[0], &fe_derived, comp, dim)); /* stress / strain (# of component = 4 or 6) */
  switch (tpe->isotropiclinearelasticity.derivedsoltype) {
    case Elasticity::DerivedSolutionType::PLANESTRAIN_2D:
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_planestrain_poroelas, &u, &derived_solutions[3], "total_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::effectivestress_planestrain_poroelas, &u, &derived_solutions[4], "effective_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_planestrain_poroelas, &u, &derived_solutions[5], "strain", time, step));
    break;
    case Elasticity::DerivedSolutionType::CAUCHY_3D:
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_3D_poroelas, &u, &derived_solutions[3], "total_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::effectivestress_3D_poroelas, &u, &derived_solutions[4], "effective_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_3D_poroelas, &u, &derived_solutions[5], "strain", time, step));
    break;
    case Elasticity::DerivedSolutionType::AXISYMMETRIC_2D_PLANESTRAIN:
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_stress_planestrain_poroelas, &u, &derived_solutions[3], "total_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_effectivestress_planestrain_poroelas, &u, &derived_solutions[4], "effective_stress", time, step));
    PetscCall(tpe->SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_strain_planestrain_poroelas, &u, &derived_solutions[5], "strain", time, step));
    break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type for elasticity");
  }
  PetscCall(DMDestroy(&dmAux)); 
  PetscCall(PetscFEDestroy(&fe_derived));
  PetscCall(VecDestroy(&derived_solutions[0]));
  PetscCall(VecDestroy(&derived_solutions[1]));
  PetscCall(VecDestroy(&derived_solutions[2]));
  PetscCall(VecDestroy(&derived_solutions[3]));
  PetscCall(VecDestroy(&derived_solutions[4]));
  PetscCall(VecDestroy(&derived_solutions[5]));

  PetscFunctionReturn(0);
}


PetscErrorCode TransientPoroElasticity::ProcessOptions(MPI_Comm comm) {
    
  PetscInt       derived_sol_elas_type = 0, rank;
  PetscBool      flg_elas;
  
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-gravity", "Use gravity in the y direction in 2D and in the z direction in 3D", NULL, this->transientdarcysflow.useGravity, &this->transientdarcysflow.useGravity, NULL));

  PetscCall(PetscOptionsBool("-explicit", "Use explicit timestepping", NULL, this->transientdarcysflow.expTS, &this->transientdarcysflow.expTS, NULL));
  PetscCall(PetscOptionsBool("-lumped", "Lump the mass matrix", NULL, this->transientdarcysflow.lumped, &this->transientdarcysflow.lumped, NULL));
  
  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));
  PetscCall(PetscOptionsBool("-material_fields", "Material properties are assigned to auxiliary fields", NULL, this->material_fields, &this->material_fields, NULL));
  PetscCall(PetscOptionsBool("-cryer", "cryer model", NULL, this->isotropiclinearelasticity.cryer, &this->isotropiclinearelasticity.cryer, NULL));
  
  PetscCall(PetscOptionsBool("-near_nullspace", "Use the rigid body modes as an AMG near nullspace", NULL, this->isotropiclinearelasticity.useNearNullspace, &this->isotropiclinearelasticity.useNearNullspace, NULL));
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

  PetscOptionsEnd();
  PetscFunctionReturn(0);

}

PetscErrorCode TransientPoroElasticity::SetupTS(DM dm, Vec *sol, const char name[], TS* ts) 
{
    char           view_call[PETSC_MAX_PATH_LEN]; 

    PetscFunctionBegin;

    if (this->transientdarcysflow.expTS) {
        PetscCall(DMTSSetRHSFunctionLocal(dm, DMPlexTSComputeRHSFunctionFEM, NULL));
        if (this->transientdarcysflow.lumped) PetscCall(DMTSCreateRHSMassMatrixLumped(dm));
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

PetscErrorCode TransientPoroElasticity::SolveFE(DM dm, Vec *sol, const char name[], TS* ts)
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 
  
  PetscFunctionBegin;

  PetscCall(TSSolve(*ts, *sol));
  PetscCall(DMTSCheckFromOptions(*ts, *sol));

  PetscFunctionReturn(0); 
}

PetscErrorCode TransientPoroElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);

    PetscScalar permeability_darcy(1e-15) /* m^2 */, fluid_density_darcy(1000.0)/* kg/m^3 */, fluid_viscosity_darcy(0.001)/* Pa s */;
    PetscScalar gravity_darcy(0.0) ;/* m/s^2 */
    // PetscScalar compressibility_darcy(5e-10) /* 1/Pa */, porosity_darcy(0.05) /*  */ ;
    PetscScalar biot_modulus_hm(1.0e10) /* Pa */, biot_coefficient_hm(1.0) ;
    
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
    material.AddProperty(shear_m);
    material.AddProperty(lambda);
    material.AddProperty(poisson_r);
    material.AddProperty(youngs_m);
    material.AddProperty(bulk_m);

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_permeability", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &permeability_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_density", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_fluid_viscosity", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_viscosity_darcy, NULL);
    if (transientdarcysflow.useGravity == PETSC_TRUE) {gravity_darcy = 9.8;}

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_biot_modulus", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_modulus_hm, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_biot_coefficient", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_coefficient_hm, NULL);
    
    Property permeability(permeability_darcy, "m2", "Permeability" ); 
    Property fluid_density(fluid_density_darcy, "kg/m3", "Fluid density" ); 
    Property fluid_viscosity(fluid_viscosity_darcy, "Pa s", "Fluid viscosity" ); 
    Property gravity(gravity_darcy, "m/s^2", "Gravity" ); 
    Property biot_modulus(biot_modulus_hm, "1/Pa", "Biot Modulus" ); 
    Property biot_coefficient(biot_coefficient_hm, "-", "Biot coefficient" ); 
       

    material.AddProperty(permeability);
    material.AddProperty(fluid_density);
    material.AddProperty(fluid_viscosity);
    material.AddProperty(gravity);
    material.AddProperty(biot_modulus);
    material.AddProperty(biot_coefficient);

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientPoroElasticity::SetMaterial(MPI_Comm comm) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);

    PetscScalar permeability_darcy(1e-15) /* m^2 */, fluid_density_darcy(1000.0)/* kg/m^3 */, fluid_viscosity_darcy(0.001)/* Pa s */;
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
    Property lambda(lambda_temp, "MPa", "Lame's first parameter, lambda" ); /*Lame's first parameter, lambda*/
    Property poisson_r(poisson_r_temp, "MPa", "Poisson's ratio, nu" ); /* Poisson's ratio, nu  */
    Property bulk_m(bulk_m_temp, "MPa", "bulk modulus, K" ); /* bulk modulus, K */
    
    material.SetName(name);
    material.SetLabel_number(label_number);
    material.SetDescription(description);
    material.AddProperty(shear_m);
    material.AddProperty(lambda);
    material.AddProperty(poisson_r);
    material.AddProperty(youngs_m);
    material.AddProperty(bulk_m);

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-permeability"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &permeability_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_density"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_density_darcy, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-fluid_viscosity"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &fluid_viscosity_darcy, NULL);
    if (transientdarcysflow.useGravity == PETSC_TRUE) {gravity_darcy = 9.8;}

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-biot_modulus"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_modulus_hm, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-biot_coefficient"));
    PetscOptionsGetScalar(NULL, NULL, option_call, &biot_coefficient_hm, NULL);
    
    Property permeability(permeability_darcy, "m2", "Permeability" ); 
    Property fluid_density(fluid_density_darcy, "kg/m3", "Fluid density" ); 
    Property fluid_viscosity(fluid_viscosity_darcy, "Pa s", "Fluid viscosity" ); 
    Property gravity(gravity_darcy, "m/s^2", "Gravity" ); 
    Property biot_modulus(biot_modulus_hm, "1/Pa", "Biot Modulus" ); 
    Property biot_coefficient(biot_coefficient_hm, "-", "Biot coefficient" ); 
       

    material.AddProperty(permeability);
    material.AddProperty(fluid_density);
    material.AddProperty(fluid_viscosity);
    material.AddProperty(gravity);
    material.AddProperty(biot_modulus);
    material.AddProperty(biot_coefficient);

    PetscFunctionReturn(0); 

}

PetscErrorCode TransientPoroElasticity::ParameterSetup(PetscDS ds)
{
  const PetscInt num = material.GetNumofProperties();
  PetscScalar constants[num];

  PetscFunctionBeginUser;

  for (int i = 0; i < num; i++)
  {
      constants[i] = material.GetProperty(i);
  }

  PetscCall(PetscDSSetConstants(ds, num, constants));

  PetscFunctionReturn(0);
}    

PetscErrorCode TransientPoroElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
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

PetscErrorCode TransientPoroElasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
{
    PetscFunctionBeginUser;
   
    isotropiclinearelasticity.BCSetup(dm, ds, wf, 0);
    transientdarcysflow.BCSetup_HM(dm, ds, wf, 1);

    PetscFunctionReturn(0);
}

PetscErrorCode TransientPoroElasticity::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_poroelas_u));
  PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_poroelas_uu));
  PetscCall(PetscDSSetJacobian(ds, 0, 1, NULL, NULL, Pw_Functions::g2_poroelas_up, NULL)); 
  PetscCall(PetscDSSetJacobian(ds, 0, 2, NULL, NULL, Pw_Functions::g2_poroelas_ustrain, NULL)); 

  PetscCall(PetscDSSetResidual(ds, 1, Pw_Functions::f0_poroelas_pressure_transient, Pw_Functions::f1_poroelas_pressure));
  PetscCall(PetscDSSetJacobian(ds, 1, 1, Pw_Functions::g0_poroelas_pressure_pressure_transient, NULL, NULL, Pw_Functions::g3_poroelas_pressure_pressure));
  PetscCall(PetscDSSetJacobian(ds, 1, 2, Pw_Functions::g0_poroelas_pressure_strain_transient, NULL, NULL, NULL));
  PetscCall(PetscDSSetExactSolution(ds, 1, Pw_Functions::initial_pressure_poroelas, (void*)&transientdarcysflow.init_pressure));

  PetscCall(PetscDSSetResidual(ds, 2, Pw_Functions::f0_poroelas_strain, NULL));
  PetscCall(PetscDSSetJacobian(ds, 2, 0, NULL, Pw_Functions::g1_poroelas_strain_u, NULL, NULL));
  PetscCall(PetscDSSetJacobian(ds, 2, 2, Pw_Functions::g0_poroelas_strain_strain, NULL, NULL, NULL));
  
  
  PetscFunctionReturn(0);
}

PetscErrorCode TransientPoroElasticity::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf));

  PetscFunctionReturn(0);
}

PetscErrorCode TransientPoroElasticity::SetupFE(DM dm, PetscDS ds, const char *name[], PetscFE *fe)
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
 
  /* pressure */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[1]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[1] ? prefix : NULL, -1, &fe[1]));
  PetscCall(PetscFECopyQuadrature(fe[0], fe[1]));
  PetscCall(PetscObjectSetName((PetscObject) fe[1], name[1]));

  /* volumetric strain */
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name[2]));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, simplex, name[2] ? prefix : NULL, -1, &fe[2]));
  PetscCall(PetscFECopyQuadrature(fe[1], fe[2]));
  PetscCall(PetscObjectSetName((PetscObject) fe[2], name[2]));
  
  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe[0]));
  PetscCall(DMSetField(dm, 1, NULL, (PetscObject) fe[1]));
  PetscCall(DMSetField(dm, 2, NULL, (PetscObject) fe[2]));

  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  // added in Transient
  // if (this->transientdarcysflow.expTS) {PetscCall(PetscDSSetImplicit(ds, 2, PETSC_FALSE));}
  if (this->transientdarcysflow.expTS) {
    for (PetscInt f = 0; f < 3; ++f) {PetscCall(PetscDSSetImplicit(ds, f, PETSC_FALSE));}
  }
  // added in Transient
  PetscCall(SetupPrimalProblem(dm, ds));
  // PetscCall(SetupPrimalProblem(dm, ds, fe, febd)); //added 230508
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

PetscErrorCode TransientPoroElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

TransientPoroElasticity::TransientPoroElasticity()
{

}