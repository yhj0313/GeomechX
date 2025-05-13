#include "TransverselyIsotropicLinearElasticity.hh"


PetscErrorCode TransverselyIsotropicLinearElasticity::SetMaterial(MPI_Comm comm, const std::string name, PetscInt label_number, std::string description) {

    Mat       SMat, SMat_rotation, CMat_rotation, IMat, TMat, tTMat; // SMat: Compliance matrix, CMat: Stiffness matrix, IMat: identity matrix, TMat: Matrix for transformation (A_epsilon, https://en.wikipedia.org/wiki/Transverse_isotropy), tTMat: transpose of TMat
    PetscInt i, i1[3] = {0, 1, 2}, j1[3] = {0, 1, 2}, i2[3] = {3, 4, 5}, j2[3] = {3, 4, 5};
    PetscScalar Identity_d = 1.0; 
    /* Transversely isotropic material with the (x, y) plane of isotropy, 
    E_x = E_y = E, E_z = E', mu_xy = mu_yx = mu, mu_zx = mu_zy = mu', mu_xz = mu_yz = mu' * E/E'
    G_xz = G_yz = G', G_xy = G */  

    PetscScalar youngs_m(1.0), youngs_m_(1.0), poisson_r(0.25), poisson_r_(0.25), shear_m_(1.0), theta(0.0);
    PetscReal radian;

    char option_call[PETSC_MAX_PATH_LEN];

    PetscFunctionBegin;   

    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_youngs_modulus", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &youngs_m, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_youngs_modulus_", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &youngs_m_, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_poissons_ratio", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &poisson_r, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_poissons_ratio_", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &poisson_r_, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_shear_modulus_", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &shear_m_, NULL);
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_theta", const_cast<char*>(name.c_str())));
    PetscOptionsGetScalar(NULL, NULL, option_call, &theta, NULL);

    radian = theta / 180.0 * PETSC_PI;

    PetscReal S_upperleft[9] = {1.0/youngs_m, -poisson_r_/youngs_m_, -poisson_r/youngs_m, -poisson_r_/youngs_m_, 1.0/youngs_m_, -poisson_r_/youngs_m_, -poisson_r/youngs_m, -poisson_r_/youngs_m_, 1.0/youngs_m};
    PetscReal S_lowerright[9] = {1.0/shear_m_, 0.0, 0.0, 0.0, 2.0*(1.0 + poisson_r) / youngs_m, 0.0, 0.0, 0.0, 1.0/shear_m_};

    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, 6, 6, NULL, &SMat);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, 6, 6, NULL, &SMat_rotation);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, 6, 6, NULL, &CMat_rotation);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, 6, 6, NULL, &IMat);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, 6, 6, NULL, &TMat);
    MatSetValues(SMat, 3, i1, 3, j1, S_upperleft, INSERT_VALUES);
    MatSetValues(SMat, 3, i2, 3, j2, S_lowerright, INSERT_VALUES);

    for (i=0; i<6; i++) {
      MatSetValue(IMat, i, i, Identity_d, INSERT_VALUES);
    }
    MatAssemblyBegin(SMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(SMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(IMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(IMat, MAT_FINAL_ASSEMBLY);
    MatViewFromOptions(SMat, NULL, "-SMat_view");
    MatViewFromOptions(IMat, NULL, "-IMat_view");

    PetscReal T_upperleft[9] = {PetscPowReal(PetscCosReal(radian),2), PetscPowReal(PetscSinReal(radian),2), 0.0, PetscPowReal(PetscSinReal(radian),2), PetscPowReal(PetscCosReal(radian),2), 0.0, 0.0, 0.0, 1.0};
    PetscReal T_upperright[9] = {0.0, 0.0, PetscCosReal(radian)*PetscSinReal(radian), 0.0, 0.0, -PetscCosReal(radian)*PetscSinReal(radian), 0.0, 0.0, 0.0};
    PetscReal T_lowerleft[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 * PetscCosReal(radian)*PetscSinReal(radian), 2.0 * PetscCosReal(radian)*PetscSinReal(radian), 0.0};
    PetscReal T_lowerright[9] = {PetscCosReal(radian), -PetscSinReal(radian), 0.0, PetscSinReal(radian), PetscCosReal(radian), 0.0, 0.0, 0.0, PetscPowReal(PetscCosReal(radian),2)-PetscPowReal(PetscSinReal(radian),2)};

    MatSetValues(TMat, 3, i1, 3, j1, T_upperleft, INSERT_VALUES);
    MatSetValues(TMat, 3, i1, 3, j2, T_upperright, INSERT_VALUES);
    MatSetValues(TMat, 3, i2, 3, j1, T_lowerleft, INSERT_VALUES);
    MatSetValues(TMat, 3, i2, 3, j2, T_lowerright, INSERT_VALUES);
    MatAssemblyBegin(TMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(TMat, MAT_FINAL_ASSEMBLY);
    MatViewFromOptions(TMat, NULL, "-TMat_view");
    MatTranspose(TMat, MAT_INITIAL_MATRIX, &tTMat);
    MatViewFromOptions(tTMat, NULL, "-tTMat_view");

    MatMatMatMult(TMat, SMat, tTMat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &SMat_rotation);
    MatAssemblyBegin(SMat_rotation, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(SMat_rotation, MAT_FINAL_ASSEMBLY);
    MatViewFromOptions(SMat_rotation, NULL, "-SMat_rotation_view"); // compliance tensor

    MatLUFactor(SMat_rotation, NULL, NULL, NULL);
    MatMatSolve(SMat_rotation, IMat, CMat_rotation);  // calculate the inverse of SMat_rotation (= CMat_rotation)


    MatViewFromOptions(CMat_rotation, NULL, "-CMat_rotation_view"); // stiffness tensor

    PetscScalar C00, C01, C02, C05, C11, C12, C15, C22, C25, C33, C34, C44, C55;
    MatGetValue(CMat_rotation, 0, 0, &C00); 
    MatGetValue(CMat_rotation, 0, 1, &C01);
    MatGetValue(CMat_rotation, 0, 2, &C02);
    MatGetValue(CMat_rotation, 0, 5, &C05);
    MatGetValue(CMat_rotation, 1, 1, &C11);
    MatGetValue(CMat_rotation, 1, 2, &C12);
    MatGetValue(CMat_rotation, 1, 5, &C15);
    MatGetValue(CMat_rotation, 2, 2, &C22);
    MatGetValue(CMat_rotation, 2, 5, &C25);
    MatGetValue(CMat_rotation, 3, 3, &C33); //yz component
    MatGetValue(CMat_rotation, 3, 4, &C34);
    MatGetValue(CMat_rotation, 4, 4, &C44); //xz component
    MatGetValue(CMat_rotation, 5, 5, &C55); //xy component

    MatDestroy(&SMat); MatDestroy(&SMat_rotation); MatDestroy(&IMat); MatDestroy(&CMat_rotation); MatDestroy(&TMat); MatDestroy(&tTMat);

    Property D11(C00, "Pa", "D11 component of stiffness tensor" );
    Property D12(C01, "Pa", "D12 component of stiffness tensor" );
    Property D13(C02, "Pa", "D13 component of stiffness tensor" );
    Property D16(C05, "Pa", "D16 component of stiffness tensor" );
    Property D22(C11, "Pa", "D22 component of stiffness tensor" );
    Property D23(C12, "Pa", "D23 component of stiffness tensor" );
    Property D26(C15, "Pa", "D26 component of stiffness tensor" );
    Property D33(C22, "Pa", "D33 component of stiffness tensor" );
    Property D36(C25, "Pa", "D36 component of stiffness tensor" );
    Property D44(C33, "Pa", "D44 component of stiffness tensor" );
    Property D45(C34, "Pa", "D45 component of stiffness tensor" );
    Property D55(C44, "Pa", "D55 component of stiffness tensor" );
    Property D66(C55, "Pa", "D66 component of stiffness tensor" );
    
    material.SetName(name);
    material.SetLabel_number(label_number);
    material.SetDescription(description);
    material.AddProperty(D11);
    material.AddProperty(D12);
    material.AddProperty(D13);
    material.AddProperty(D16);
    material.AddProperty(D22);
    material.AddProperty(D23);
    material.AddProperty(D26);
    material.AddProperty(D33);
    material.AddProperty(D36);
    material.AddProperty(D44);
    material.AddProperty(D45);
    material.AddProperty(D55);
    material.AddProperty(D66);

    PetscFunctionReturn(0); 

}

PetscErrorCode TransverselyIsotropicLinearElasticity::DomainSetup(PetscDS ds) {

  PetscFunctionBeginUser;
  PetscCall(PetscDSSetResidual(ds, 0, NULL, Pw_Functions::f1_TI));
  PetscCall(PetscDSSetJacobian(ds, 0, 0, NULL, NULL, NULL, Pw_Functions::g3_TI));
  PetscFunctionReturn(0);
}

PetscErrorCode TransverselyIsotropicLinearElasticity::SetupFE(DM dm, PetscDS ds, const char name[], PetscFE *fe)
{
  // AppCtx        *user = (AppCtx *) ctx;
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

PetscErrorCode TransverselyIsotropicLinearElasticity::InputParameterOptions(MPI_Comm comm, const std::string name)
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
  para = this->inputpara;

  PetscCall(PetscOptionsEList(option_name.c_str(), "Type of two input elastic parameters", NULL, this->inputParameterTypes, TransverselyIsotropicLinearElasticity::InputParameterType::NUM_INPUT_PARA_TYPES, this->inputParameterTypes[this->inputpara], &para, &flg));
  if (flg)
  {
    this->inputpara = (TransverselyIsotropicLinearElasticity::InputParameterType)para;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Input type for two elastic parameters is not given. Five parameters are assumed to be given." << std::endl;
    }
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

PetscErrorCode TransverselyIsotropicLinearElasticity::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols)
{
  // // add 
  // PetscInt       n;
  
  PetscFunctionBeginUser;

  switch (this->derivedsoltype) {
    case CAUCHY_3D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystrain_3D, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::cauchystress_3d_TI, sol, &derivedsols[1], "stress");
    break;
    case AXISYMMETRIC_3D:
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_3d_strain_TI, sol, &derivedsols[0], "strain");
    SolveDerivedField(dm, &dmAux, Pw_Functions::axisymmetric_3d_stress_TI, sol, &derivedsols[1], "stress");
    break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type");

  }

  PetscFunctionReturn(0);

  
}

