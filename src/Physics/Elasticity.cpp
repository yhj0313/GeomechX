#include "Elasticity.hh"

PetscErrorCode Elasticity::CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace)
{
  PetscFunctionBegin;
  PetscCall(DMPlexCreateRigidBody(dm, origField, nullspace));
  PetscFunctionReturn(0);
}

PetscErrorCode Elasticity::SolveFE(DM dm, Vec *sol, const char name[], SNES* snes)
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

PetscErrorCode Elasticity::CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim)
{
  PetscBool      simplex;
  DM             coordDM;
  char           prefix[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;

  DMPlexIsSimplex(dm, &simplex);
  PetscCall(PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name));
  PetscCall(PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, numComp, simplex, name ? prefix : NULL, -1, fe_derived));
  PetscCall(PetscObjectSetName((PetscObject) *fe_derived, name));
  PetscCall(PetscFECopyQuadrature(*fe, *fe_derived));
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  PetscCall(DMClone(dm, dmAux));
  PetscCall(DMSetCoordinateDM(*dmAux, coordDM));
  PetscCall(DMSetField(*dmAux, 0, NULL, (PetscObject) *fe_derived));
  PetscCall(DMCreateDS(*dmAux));

  PetscFunctionReturn(0);  
}

PetscErrorCode Elasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[])
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 
  // // add 
  // PetscInt       n;

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, 0.0, *sol, &func, INSERT_VALUES, *derivedsol));
  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, *derivedsol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  

  // // add
  // PetscCall(DMGetNumAuxiliaryVec(dm, &n));
  // // std::cout << n << "DMGetNumAuxiliaryVec" << std::endl;
  

  PetscFunctionReturn(0);  
}

PetscErrorCode Elasticity::ProcessOptions(MPI_Comm comm)
{
  PetscInt       derived_sol_type = 0, rank; 
  PetscBool      flg;
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-near_nullspace", "Use the rigid body modes as an AMG near nullspace", NULL, this->useNearNullspace, &this->useNearNullspace, NULL));
  // PetscCall(PetscOptionsBool("-radial_tangential_2D", "Calculate radial and tangential stress and strain", NULL, this->radial_tangential_2D, &this->radial_tangential_2D, NULL));
  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));
  PetscCall(PetscOptionsBool("-TI_stripload", "Strip load test of transversely isotropic rock core sample", NULL, this->TIStripLineLoad, &this->TIStripLineLoad, NULL));
  PetscCall(PetscOptionsBool("-cryermodel", "3D cryer model", NULL, this->cryer, &this->cryer, NULL));

  derived_sol_type = this->derivedsoltype;

  PetscCall(PetscOptionsEList("-derived_solution_type", "Type of derived solutions", NULL, this->DerivedSolutionTypes, NUM_DERIVED_SOLUTION_TYPES, this->DerivedSolutionTypes[this->derivedsoltype], &derived_sol_type, &flg));
  if (flg)
  {
    this->derivedsoltype = (DerivedSolutionType)derived_sol_type;
  }
  else
  {
    CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0)
    {
      std::cout << "ERROR: Type of derived solutions is not given. " <<  DerivedSolutionTypes[0] << " is assumed." << std::endl;
    }
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

PetscErrorCode Elasticity::ParameterSetup(PetscDS ds)
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


PetscErrorCode Elasticity::SetupPrimalProblem(DM dm, PetscDS ds)
{
  PetscWeakForm    wf;
  PetscFunctionBeginUser;

  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetWeakForm(ds, &wf));
  PetscCall(DomainSetup(ds));
  PetscCall(BCSetup(dm, ds, wf, 0));

  PetscFunctionReturn(0);
}

PetscErrorCode Elasticity::CreateAuxiliaryVecv22(DM dm, DM *auxdm, Vec *la, std::vector<BCElasticity> bc)
{
  PetscErrorCode (**afuncs)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar *, void *);
  
  PetscBool simplex;
  PetscInt  dim, Nf, f;

  PetscFE  fe;
  MPI_Comm comm;

  const PetscInt id_2d[] = {1, 2, 3, 4};
  const PetscInt id_3d[] = {1, 2, 5, 6, 3, 4};

  DMLabel label;
  PetscInt val[3] = {1, 2, 3};
  PetscInt val_1[1] = {1};
  PetscInt val_2[1] = {2};
  PetscInt val_3[1] = {3};
  PetscInt numAuxLabels;
  const char *labelName;

  PetscInt vec_size;

  PetscFunctionBeginUser;
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMPlexIsSimplex(dm, &simplex));
  PetscCall(DMGetNumFields(dm, &Nf));
  PetscCall(PetscMalloc1(Nf, &afuncs));
  for (f = 0; f < Nf; ++f) afuncs[f] = Pw_Functions::boundary_normalstress; //take this fuction just to test
  PetscCall(DMClone(dm, auxdm));
  // PetscCall(SetupDiscretization(*auxdm, dim, simplex, user));
  PetscCall(PetscObjectGetComm((PetscObject)*auxdm, &comm));
  PetscCall(PetscFECreateDefault(comm, dim, 1, simplex, "normalstress_boundary_", -1, &fe));
  PetscCall(PetscObjectSetName((PetscObject)fe, "normalstress_boundary"));
  PetscCall(DMSetField(*auxdm, 0, NULL, (PetscObject)fe));
  PetscCall(PetscFEDestroy(&fe));
  PetscCall(DMCreateDS(*auxdm));

  PetscCall(PetscObjectSetName((PetscObject)*auxdm, "aux_dm"));

  PetscCall(DMGetLabel(*auxdm, "Faces", &label));
  PetscCall(DMPlexLabelAddFaceCells(*auxdm, label));
  
  


  PetscCall(DMGetNumLabels(*auxdm, &numAuxLabels));
  PetscPrintf(PETSC_COMM_WORLD, "numLabels is %d \n", numAuxLabels);
  for (int i = 0 ; i < numAuxLabels ; i++) {
    PetscCall(DMGetLabelName(*auxdm, i, &labelName));
    PetscPrintf(PETSC_COMM_WORLD, "LabelName is %s \n", labelName);
  }

  PetscCall(DMCreateGlobalVector(*auxdm, la));
  // PetscCall(DMGetLabel(*auxdm, "aux", &label));
  PetscCall(DMGetLabel(*auxdm, "Faces", &label));

  PetscCall(DMViewFromOptions(*auxdm, NULL, "-aux_dm_view"));
  
  vec_size = bc.size();
  auto beg = bc.begin();
  PetscScalar *bdvalue;
 
  for (int i = 0 ; i < vec_size ; ++i) {
    PetscPrintf(PETSC_COMM_WORLD, "bc[i].GetBdValue() is %f, and i is %d, and vecsize is %d. \n", bc[i].GetBdValue(), i, vec_size);
      
    if (bc[i].GetBCType() == BCElasticity::BoundaryConditionType::NORMAL_STRESS) {
      *bdvalue = bc[i].GetBdValue(); 
      PetscPrintf(PETSC_COMM_WORLD, "bdvalue (bc[i].GetBdValue()) is %f, and i is %d, and vecsize is %d. \n", *bdvalue, i, vec_size);
 
      if (dim == 2) PetscCall(DMProjectFunctionLabel(*auxdm, 0.0, label, 1, &id_2d[i], 0, NULL, afuncs, (void**)(&bdvalue), INSERT_VALUES, *la));
      else if (dim == 3) PetscCall(DMProjectFunctionLabel(*auxdm, 0.0, label, 1, &id_3d[i], 0, NULL, afuncs, (void**)(&bdvalue), INSERT_VALUES, *la));
    }
    
  }

  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, *la));
  PetscCall(VecViewFromOptions(*la, NULL, "-local_aux_view"));
  PetscCall(PetscFree(afuncs));
  PetscFunctionReturn(0);
}

PetscErrorCode Elasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf)
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
    const char *bctype_call[6] = {"-bctype_elas_left", "-bctype_elas_right", "-bctype_elas_bottom", "-bctype_elas_top", "-bctype_elas_front", "-bctype_elas_back"};
    const char *bdvalue_call[6] = {"-bdvalue_elas_left", "-bdvalue_elas_right", "-bdvalue_elas_bottom", "-bdvalue_elas_top", "-bdvalue_elas_front", "-bdvalue_elas_back"};
    const char *bdvalues_call[6] = {"-bdvalues_elas_left", "-bdvalues_elas_right", "-bdvalues_elas_bottom", "-bdvalues_elas_top", "-bdvalues_elas_front", "-bdvalues_elas_back"};
     

    void (*normal_stress[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_normal_stress_bd_0, Pw_Functions::f0_normal_stress_bd_1, Pw_Functions::f0_normal_stress_bd_2, Pw_Functions::f0_normal_stress_bd_3, Pw_Functions::f0_normal_stress_bd_4, Pw_Functions::f0_normal_stress_bd_5};
    void (*traction[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_traction_bd_0, Pw_Functions::f0_traction_bd_1, Pw_Functions::f0_traction_bd_2, Pw_Functions::f0_traction_bd_3, Pw_Functions::f0_traction_bd_4, Pw_Functions::f0_traction_bd_5};



    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCElasticity(0)); /* Default boundary condition is the first option. (FREE in Isotropic linear elascticity) */
    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");

    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}

    for (int i = 0; i < n_bd; i++)
    {
      PetscCall(BoundaryConditions[i].SetBdIndex(i));
      PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCElasticity::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
      if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        
        if (bc == 3) {  // displacement boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;

          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 4) { // normal stress boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            std::cout << BoundaryConditions[0].GetBdValue() << "is BoundaryConditions[0].GetBdValue()." << std::endl;
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            std::cout << BoundaryConditions[i].GetBdValue() << "is value and i is " << i << std::endl;
            std::cout << BoundaryConditions[0].GetBdValue() << "is BoundaryConditions[0].GetBdValue()." << std::endl;
            flg_bdvalue = PETSC_FALSE;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}

        }
        else if (bc == 5) { // traction boundary

          PetscOptionsGetScalarArray(NULL, NULL, bdvalues_call[i], bdvs, &dim, &flg_bdvalues);
          if (flg_bdvalues) {
            PetscCall(BoundaryConditions[i].SetBdValues(bdvs, dim));
            flg_bdvalues = PETSC_FALSE;

            }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary values '%s' are not given. \n", bdvalues_call[i]);}
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
        case BCElasticity::BoundaryConditionType::FREE:
        break;
        case BCElasticity::BoundaryConditionType::ROLLER:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::FIXED:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, dim, components_3d, (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::DISPLACEMENT:
        // PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::disp_bd, NULL, NULL, NULL));
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_displacement, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCElasticity::BoundaryConditionType::NORMAL_STRESS:
        *BDVALUE = BoundaryConditions[i].GetBdValue();
        {PetscPrintf(PETSC_COMM_WORLD, "BDVALUE is %f \n", *BDVALUE);}
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        // PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, normal_stress[i], 0, NULL));
        // PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_normal_stress_bd, 0, NULL));
        // PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, BoundaryConditions[i].f0_normal_stress_bd, 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, Pw_Functions::f0_normal_stress_bd_aux, 0, NULL));
        break;
        case BCElasticity::BoundaryConditionType::TRACTION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], 0, 0, 0, traction[i], 0, NULL));
        break;
        default:
        break;
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCElasticity::BoundaryConditionType::FREE:
        break;
        case BCElasticity::BoundaryConditionType::ROLLER:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::FIXED:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, dim, components_2d, (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::DISPLACEMENT:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_displacement, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCElasticity::BoundaryConditionType::NORMAL_STRESS:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        // PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, normal_stress[i], 0, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, Pw_Functions::f0_normal_stress_bd_aux, 0, NULL));
        break;
        case BCElasticity::BoundaryConditionType::TRACTION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], 0, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], 0, 0, 0, traction[i], 0, NULL));
        break;
        default:
        break;
        }
        }
    }

    for (int i = 0 ; i < 8 ; i++) {
    PetscPrintf(PETSC_COMM_WORLD, "BoundaryConditions[i].GetBdValue() is %f, and, BoundaryConditions[i].GetBCType() is %d , and i is %d. \n", BoundaryConditions[i].GetBdValue(), BoundaryConditions[i].GetBCType(), i );
    }

    DM dmAuxiliary;
    Vec       la;
    PetscCall(CreateAuxiliaryVecv22(dm, &dmAuxiliary, &la, BoundaryConditions));

    PetscCall(DMDestroy(&dmAuxiliary));
    PetscCall(VecDestroy(&la));
    

    PetscFunctionReturn(0);
}

PetscErrorCode Elasticity::BCSetup(DM dm, PetscDS ds, PetscWeakForm wf, PetscInt field)
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
    // const char *boundary_call[6] = {"-bc_type_left", "-bc_type_right", "-bc_type_front", "-bc_type_back", "-bc_type_bottom", "-bc_type_top"};
    const char *bctype_call[6] = {"-bctype_elas_left", "-bctype_elas_right", "-bctype_elas_bottom", "-bctype_elas_top", "-bctype_elas_front", "-bctype_elas_back"};
    const char *bdvalue_call[6] = {"-bdvalue_elas_left", "-bdvalue_elas_right", "-bdvalue_elas_bottom", "-bdvalue_elas_top", "-bdvalue_elas_front", "-bdvalue_elas_back"};
    const char *bdvalues_call[6] = {"-bdvalues_elas_left", "-bdvalues_elas_right", "-bdvalues_elas_bottom", "-bdvalues_elas_top", "-bdvalues_elas_front", "-bdvalues_elas_back"};
    
    void (*normal_stress[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_normal_stress_bd_0, Pw_Functions::f0_normal_stress_bd_1, Pw_Functions::f0_normal_stress_bd_2, Pw_Functions::f0_normal_stress_bd_3, Pw_Functions::f0_normal_stress_bd_4, Pw_Functions::f0_normal_stress_bd_5};
    void (*traction[6])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[]) = {Pw_Functions::f0_traction_bd_0, Pw_Functions::f0_traction_bd_1, Pw_Functions::f0_traction_bd_2, Pw_Functions::f0_traction_bd_3, Pw_Functions::f0_traction_bd_4, Pw_Functions::f0_traction_bd_5};

    
    // std::vector<BoundaryConditionType> bctypes;

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "Faces", &label));
    PetscCall(PetscDSGetSpatialDimension(ds, &dim));

    BoundaryConditions.assign(6, BCElasticity(0)); /* Default boundary condition is the first option. (FREE in Isotropic linear elascticity) */

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Boundary Condition Options", "DMPLEX");

    PetscCall(PetscMalloc1(dim, &bdvs));

    if (dim == 2) {n_bd = 4;}
    else if (dim == 3) {n_bd = 6;}

    for (int i = 0; i < n_bd; i++)
    {
      PetscCall(BoundaryConditions[i].SetBdIndex(i));
      PetscCall(PetscOptionsEList(bctype_call[i], "Type of boundary condition", NULL, BCTypes, BCElasticity::BoundaryConditionType::NUM_BC_TYPES, BCTypes[0], &bc, &flg_bctype));
      if (flg_bctype) {
        PetscCall(BoundaryConditions[i].SetBCType(bc));
        // std::cout << BoundaryConditions[i].GetBCType() << std::endl;
        
        if (bc == 3) {  // displacement boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
            // bdv = 0.0;
            // std::cout << BoundaryConditions[i].GetBdValue() << std::endl;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
        }
        else if (bc == 4) { // normal stress boundary
          PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          if (flg_bdvalue) {
            PetscCall(BoundaryConditions[i].SetBdValue(bdv));
            flg_bdvalue = PETSC_FALSE;
          }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
          // PetscOptionsGetScalarArray(NULL, NULL, bdvalues_call[i], bdvs, &dim, &flg_bdvalues);
          //  if (flg_bdvalues) {
          //   PetscCall(BoundaryConditions[i].SetBdValues(bdvs, dim));
          //   flg_bdvalues = PETSC_FALSE;
          //   // bdv = 0.0;
          //   // std::cout << BoundaryConditions[i].GetBdValues() << std::endl;
        }
        else if (bc == 5) { // traction boundary
          // PetscOptionsGetScalar(NULL, NULL, bdvalue_call[i], &bdv, &flg_bdvalue);
          // if (flg_bdvalue) {
          //   PetscCall(BoundaryConditions[i].SetBdValue(bdv));
          //   flg_bdvalue = PETSC_FALSE;
          // }
          // else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary value '%s' is not given. \n", bdvalue_call[i]);}
          PetscOptionsGetScalarArray(NULL, NULL, bdvalues_call[i], bdvs, &dim, &flg_bdvalues);
          if (flg_bdvalues) {
            PetscCall(BoundaryConditions[i].SetBdValues(bdvs, dim));
            flg_bdvalues = PETSC_FALSE;
            // bdv = 0.0;
            // std::cout << BoundaryConditions[i].GetBdValues() << std::endl;
            }
          else {PetscPrintf(PETSC_COMM_WORLD, "ERROR: Boundary values '%s' are not given. \n", bdvalues_call[i]);}
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
        case BCElasticity::BoundaryConditionType::FREE:
        break;
        case BCElasticity::BoundaryConditionType::ROLLER:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::FIXED:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, dim, components_3d, (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::DISPLACEMENT:
        // PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::disp_bd, NULL, NULL, NULL));
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 1, &cmp_3d[i], (void (*)(void))Pw_Functions::boundary_displacement, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCElasticity::BoundaryConditionType::NORMAL_STRESS:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, normal_stress[i], 0, NULL));
        break;
        case BCElasticity::BoundaryConditionType::TRACTION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_3d[i]], label, 1, &id_3d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_3d[i], field, 0, 0, traction[i], 0, NULL));
        break;
        default:
        break;
        }
        }
    }
    else if (dim == 2) {
        for (int i = 0; i < 4; i++)
        {
        switch (BoundaryConditions[i].GetBCType())
        {
        case BCElasticity::BoundaryConditionType::FREE:
        break;
        case BCElasticity::BoundaryConditionType::ROLLER:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::FIXED:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, dim, components_2d, (void (*)(void))Pw_Functions::zero_disp, NULL, NULL, NULL));
        break;
        case BCElasticity::BoundaryConditionType::DISPLACEMENT:
        PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 1, &cmp_2d[i], (void (*)(void))Pw_Functions::boundary_displacement, NULL, (void*)(&BoundaryConditions[i]), NULL));
        break;
        case BCElasticity::BoundaryConditionType::NORMAL_STRESS:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, normal_stress[i], 0, NULL));
        break;
        case BCElasticity::BoundaryConditionType::TRACTION:
        PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, boundary_name[id_2d[i]], label, 1, &id_2d[i], field, 0, NULL, (void (*)(void))NULL, NULL, NULL, &bd));
        PetscCall(PetscDSGetBoundary(ds, bd, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
        PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, id_2d[i], field, 0, 0, traction[i], 0, NULL));
        break;
        default:
        break;
        }
        }
    }

    PetscFunctionReturn(0);
}
   

Elasticity::DerivedSolutionType Elasticity::GetDerivedSolType() {
  return this->derivedsoltype;
}

Elasticity::Elasticity() {
  this->useNearNullspace = PETSC_TRUE;
  // this->radial_tangential_2D = PETSC_FALSE;
  this->material_view = PETSC_FALSE;
  this->TIStripLineLoad = PETSC_FALSE;
  this->cryer = PETSC_FALSE;
}
  
