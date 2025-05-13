#include "IsotropicLinearElasticity.hh"


PetscErrorCode IsotropicLinearElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    PetscScalar shear_m_temp(1.0), lambda_temp(1.0), poisson_r_temp(1.0), youngs_m_temp(1.0), bulk_m_temp(1.0);
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

    
    switch (this->inputpara) {
    case YOUNGSMODULUS_POISSONSRATIO:
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      shear_m_temp = youngs_m_temp / (2.0*(1.0+poisson_r_temp));
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      break;
    case YOUNGSMODULUS_SHEARMODULUS:
      poisson_r_temp = youngs_m_temp / (2.0*shear_m_temp) -1.0 ;
      bulk_m_temp = youngs_m_temp / (3.0*(1.0-2.0*poisson_r_temp));
      lambda_temp = youngs_m_temp * poisson_r_temp / ((1.0+poisson_r_temp)*(1.0-2.0*poisson_r_temp));
      break;
    case BULKMODULUS_SHAERMODULUS:
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      lambda_temp = bulk_m_temp - (2.0/3.0*shear_m_temp);
      break;
    case LAMESPARAMETERS:
      bulk_m_temp = lambda_temp + (2.0 / 3.0 * shear_m_temp);
      youngs_m_temp = (9.0 * shear_m_temp * bulk_m_temp ) / (3.0 * bulk_m_temp + shear_m_temp);
      poisson_r_temp = (3.0 * bulk_m_temp - 2.0 * shear_m_temp)/(6.0 * bulk_m_temp + 2.0 * shear_m_temp); 
      break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid input parameter type: %s", inputParameterTypes[PetscMin(this->inputpara, NUM_INPUT_PARA_TYPES)]);
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

    PetscFunctionReturn(0); 

}

PetscErrorCode IsotropicLinearElasticity::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_elas_u));
  PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_elas_uu));
  PetscFunctionReturn(0);
}


PetscErrorCode IsotropicLinearElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
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
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, dim, simplex, name ? prefix : NULL, -1, fe));
  PetscCall(PetscObjectSetName((PetscObject) *fe, name));
  /* Set discretization and boundary conditions for each mesh */
  PetscCall(DMSetField(dm, 0, NULL, (PetscObject) *fe));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(SetupPrimalProblem(dm, ds));
  PetscCall(ParameterSetup(ds));

  CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  if ((this->material_view) & (rank == 0)) PetscCall(material.ViewMaterial());
  
  while (cdm) {
    PetscCall(DMCopyDisc(dm, cdm));
    if (this->useNearNullspace) PetscCall(DMSetNearNullSpaceConstructor(cdm, 0, CreateElasticityNullSpace));
    /* TODO: Check whether the boundary of coarse meshes is marked */
    PetscCall(DMGetCoarseDM(cdm, &cdm));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode IsotropicLinearElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
{
  PetscInt       para = 0; 
  std::string option_name = "-";
  PetscBool      flg;
  PetscMPIInt    rank;

  PetscFunctionBeginUser;
  option_name.append(name);
  option_name.append("_inputpara_type");
  PetscOptionsBegin(comm, "", "Type of two input elastic parameters", "");
  para = this->inputpara;

  PetscCall(PetscOptionsEList(option_name.c_str(), "Type of two input elastic parameters", NULL, this->inputParameterTypes, IsotropicLinearElasticity::InputParameterType::NUM_INPUT_PARA_TYPES, this->inputParameterTypes[this->inputpara], &para, &flg));
  if (flg)
  {
    this->inputpara = (IsotropicLinearElasticity::InputParameterType)para;
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


PetscErrorCode IsotropicLinearElasticity::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols)
{

  PetscInt       n;
  
  PetscFunctionBeginUser;

  switch (this->derivedsoltype) {
    case PLANESTRAIN_2D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_planestrain, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_planestrain, sol, &derivedsols[1], "stress");
    break;
    case CAUCHY_3D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_3D, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_3D, sol, &derivedsols[1], "stress");
    break;
    case AXISYMMETRIC_2D_PLANESTRAIN:
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_strain_planestrain, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_2d_stress_planestrain, sol, &derivedsols[1], "stress");
    break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type");

  }
  
  PetscFunctionReturn(0);

  
}
