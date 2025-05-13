#include "Thermal.hh"

void Thermal::AddMaterial(Material& mtrl){
    this->materials.push_back(mtrl);
}

PetscInt Thermal::GetNumofMaterials(){
    return this->materials.size();
}

PetscErrorCode Thermal::SetupMaterial(DM dm, DM dmAux)
{
  Vec   A;
  DMLabel material;
  const PetscInt numofmaterial = this->GetNumofMaterials();
  const PetscInt numofproperties = this->materials[0].GetNumofProperties();
  PetscErrorCode (*matFuncs[numofproperties])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar u[], void *ctx);
  void *ctxs[numofproperties];
  PetscInt id;
  PetscScalar properties[numofproperties];

  PetscFunctionBegin;
  for (int i= 0; i < numofproperties; i++) {
    matFuncs[i] = Pw_Functions_common::material_property;
  }
  PetscCall(DMGetLabel(dmAux, "Cell Sets", &material));
  PetscCall(DMPlexLabelComplete(dmAux, material));
  PetscCall(DMCreateLocalVector(dmAux, &A));
  for (id = 1; id < numofmaterial+1; id++) {
    for (int j = 0; j < numofproperties; j++){
      // ctxs[j] = (void*)(materials[id].GetProperty(j));
      properties[j] = materials[id-1].GetProperty(j);
      ctxs[j] = &(properties[j]);
      }
    PetscCall(DMProjectFunctionLabelLocal(dmAux, 0.0, material, 1, &id, PETSC_DETERMINE, NULL, matFuncs, ctxs, INSERT_ALL_VALUES, A));
    }
  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, A));
  PetscCall(VecDestroy(&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Thermal::SetupMaterial(DM dm, DM dmAux, PetscScalar time, PetscInt step, Vec *sol)
{
  Vec A;
  DMLabel material;
  const PetscInt numofmaterial = this->GetNumofMaterials();
  const PetscInt numofproperties = this->materials[0].GetNumofProperties();
  PetscErrorCode (*matFuncs[numofproperties])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar u[], void *ctx);
  void *ctxs[numofproperties];
  PetscInt id;
  PetscScalar properties[numofproperties];
  // PetscPointFunc pointfuncs[numofproperties];
  std::string prop_equations[numofproperties];
  // 추가
  void (*pointfuncs_user[numofproperties])(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar property[]);

  PetscInt num_funcs(0);

  PetscFunctionBegin;

  PetscCall(DMGetLabel(dmAux, "Cell Sets", &material));
  PetscCall(DMPlexLabelComplete(dmAux, material));
  PetscCall(DMCreateGlobalVector(dmAux, &A));
  if (time == 0.0)
    PetscCall(VecSet(A, 0.0));
  PetscCall(PetscObjectSetName((PetscObject)A, "Material"));


  for (id = 1; id < numofmaterial + 1; id++)
  {
    for (int j = 0; j < numofproperties; j++)
    {
      if (j == HEAT_SOURCE_T) // heat source
      {
        if ((materials[id - 1].GetPropertyDescription(j)) == "Heat Source PWR")
        { // for high-level nuclear waste disposal modeling
          num_funcs++;
          properties[j] = materials[id - 1].GetProperty(j);
          ctxs[j] = &(properties[j]); // Number of PWR canisters per unit volume (m3)
          pointfuncs_user[j] = Pw_Functions_common::pwr_heatdecay; // to update decay heat
          matFuncs[j] = Pw_Functions_common::material_property_nochange;
        }
        else if ((materials[id - 1].GetPropertyEquation(j)) == "")
        {
          properties[j] = materials[id - 1].GetProperty(j);
          // prop_equations[j] = "";
          ctxs[j] = &(properties[j]);
          // pointfuncs[j] = Pw_Functions::nochange_prop ;
          // 추가
          pointfuncs_user[j] = Pw_Functions_common::nochange_prop_user;
          matFuncs[j] = Pw_Functions_common::material_property;
        }
        else
        {
          num_funcs++;
          // pointfuncs[j] = Pw_Functions::lin_thermal_expan_c;
          // pointfuncs[j] = materials[id-1].GetPropertyFunction(j);
          prop_equations[j] = materials[id - 1].GetPropertyEquation(j);
          // 추가
          // ctxs[j] = &(properties[j]);
          ctxs[j] = &(prop_equations[j]);
          pointfuncs_user[j] = Pw_Functions::property_user_T;
          matFuncs[j] = Pw_Functions_common::material_property_nochange;
        }
      }
      // if (!(materials[id-1].GetPropertyFunction(j))) {
      else
      {
        if ((materials[id - 1].GetPropertyEquation(j)) == "")
        {
          properties[j] = materials[id - 1].GetProperty(j);
          // prop_equations[j] = "";
          ctxs[j] = &(properties[j]);
          // pointfuncs[j] = Pw_Functions::nochange_prop ;
          // 추가
          pointfuncs_user[j] = Pw_Functions_common::nochange_prop_user;
          matFuncs[j] = Pw_Functions_common::material_property;
        }
        else
        {
          num_funcs++;
          // pointfuncs[j] = Pw_Functions::lin_thermal_expan_c;
          // pointfuncs[j] = materials[id-1].GetPropertyFunction(j);
          prop_equations[j] = materials[id - 1].GetPropertyEquation(j);
          // 추가
          // ctxs[j] = &(properties[j]);
          ctxs[j] = &(prop_equations[j]);
          pointfuncs_user[j] = Pw_Functions::property_user_T;
          matFuncs[j] = Pw_Functions_common::material_property_nochange;
        }
      }
    }
    // PetscCall(DMSetApplicationContext(dmAux, &ctxs)); /// 추가
    // PetscCall(DMProjectFunctionLabel(dmAux, time, material, 1, &id, PETSC_DETERMINE, NULL, matFuncs, ctxs, INSERT_ALL_VALUES, A));
    // DMProjectFunctionLabel 에서 vecset 0.0 을 지움
    {
      Vec localX;
      PetscCall(DMGetLocalVector(dmAux, &localX));
      // PetscCall(VecSet(localX, 0.));
      PetscCall(DMProjectFunctionLabelLocal(dmAux, time, material, 1, &id, PETSC_DETERMINE, NULL, matFuncs, ctxs, INSERT_VALUES, localX));
      PetscCall(DMLocalToGlobalBegin(dmAux, localX, INSERT_VALUES, A));
      PetscCall(DMLocalToGlobalEnd(dmAux, localX, INSERT_VALUES, A));
      PetscCall(DMRestoreLocalVector(dmAux, &localX));
    }

    if (num_funcs > 0)
    {
      // PetscCall(DMProjectFieldLabel(dmAux, time, material, 1, &id, PETSC_DETERMINE, NULL, *sol, pointfuncs, INSERT_VALUES, A));
      // DMProjectFieldLabel에서 VecSet 0.0 을 지움
      DM dmIn;
      Vec localU, localX;

      PetscCall(VecGetDM(*sol, &dmIn));
      PetscCall(DMGetLocalVector(dmIn, &localU));
      PetscCall(DMGetLocalVector(dmAux, &localX));
      // PetscCall(VecSet(localX, 0.));
      PetscCall(DMGlobalToLocalBegin(dmIn, *sol, INSERT_VALUES, localU));
      PetscCall(DMGlobalToLocalEnd(dmIn, *sol, INSERT_VALUES, localU));
      // PetscCall(DMProjectFieldLabelLocal(dmAux, time, material, 1, &id, PETSC_DETERMINE, NULL, localU, pointfuncs, INSERT_VALUES, localX));
      PetscCall(DMProjectFieldLabelLocal_User(dmAux, time, material, 1, &id, PETSC_DETERMINE, NULL, localU, pointfuncs_user, ctxs, INSERT_VALUES, localX));
      PetscCall(DMLocalToGlobalBegin(dmAux, localX, INSERT_VALUES, A));
      PetscCall(DMLocalToGlobalEnd(dmAux, localX, INSERT_VALUES, A));
      PetscCall(DMRestoreLocalVector(dmAux, &localX));
      PetscCall(DMRestoreLocalVector(dmIn, &localU));
    }
    num_funcs = 0;
  }

  PetscCall(DMSetOutputSequenceNumber(dmAux, step, time));
  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, A));
  PetscCall(VecViewFromOptions(A, NULL, "-materials_properties_view"));
  PetscCall(VecDestroy(&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}



PetscErrorCode Thermal::SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[])
{
  DM       dmAux, coordDM;
  PetscInt f;

  PetscFunctionBegin;
  /* MUST call DMGetCoordinateDM() in order to get p4est setup if present */
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  if (!feAux) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DMClone(dm, &dmAux));
  PetscCall(DMSetCoordinateDM(dmAux, coordDM));
  for (f = 0; f < NfAux; ++f) PetscCall(DMSetField(dmAux, f, NULL, (PetscObject)feAux[f]));
  PetscCall(DMCreateDS(dmAux));
  PetscCall(SetupMaterial(dm, dmAux));
  PetscCall(DMDestroy(&dmAux));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Thermal::SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step)
{
  DM dmAux, coordDM;
  PetscInt f;
  Vec *vec_sol;

  PetscFunctionBegin;
  /* MUST call DMGetCoordinateDM() in order to get p4est setup if present */
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  if (!feAux)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DMClone(dm, &dmAux));
  PetscCall(DMSetCoordinateDM(dmAux, coordDM));
  for (f = 0; f < NfAux; ++f) PetscCall(DMSetField(dmAux, f, NULL, (PetscObject)feAux[f]));
  PetscCall(DMCreateDS(dmAux));
  // PetscCall(SetupMaterial(dm, dmAux, time));
  PetscCall(DMGetGlobalVector(dm, vec_sol));
  PetscCall(SetupMaterial(dm, dmAux, time, step, vec_sol));
  PetscCall(DMDestroy(&dmAux));
  // PetscCall(VecDestroy(vec_sol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Thermal::CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim)
{
  PetscBool      simplex;
  DM             coordDM;
  char           prefix[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;

  PetscCall(DMPlexIsSimplex(dm, &simplex));
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

PetscErrorCode Thermal::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[])
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  // PetscCall(DMProjectField(*dmAux, 0.0, *sol, &func, INSERT_VALUES, *derivedsol));
  // PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, *derivedsol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  
  
  PetscFunctionReturn(0);  
}

PetscErrorCode Thermal::ProcessOptions(MPI_Comm comm) {
    
  PetscInt       derived_sol_type = 0, rank; 
  PetscBool      flg;
  
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));
  PetscCall(PetscOptionsBool("-material_fields", "Material properties are assigned to auxiliary fields", NULL, this->material_fields, &this->material_fields, NULL));
  
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
      std::cout << "ERROR: Type of derived solutions is not given. " <<  DerivedSolutionTypes[0] << " is assumed." << std::endl;}
  }

  PetscOptionsEnd();
  PetscFunctionReturn(0);

}

PetscErrorCode Thermal::ParameterSetup(PetscDS ds)
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

PetscErrorCode Thermal::InputParameterOptions(MPI_Comm comm, const std::string name)
{
  PetscFunctionBeginUser;
//   option_name.append(name);
//   option_name.append("_inputpara_type");
//   // std::cout << option_name << std::endl;
//   ierr = PetscOptionsBegin(comm, "", "Type of two input elastic parameters", "");PetscCall(ierr);
//   para = this->inputpara;

//   PetscCall(PetscOptionsEList(option_name.c_str(), "Type of two input elastic parameters", NULL, this->inputParameterTypes, IsotropicLinearElasticity::InputParameterType::NUM_INPUT_PARA_TYPES, this->inputParameterTypes[this->inputpara], &para, &flg));
//   if (flg)
//   {
//     this->inputpara = (IsotropicLinearElasticity::InputParameterType)para;
//   }
//   else
//   {
//     CHKERRMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
//     if (rank == 0)
//     {
//       std::cout << "ERROR: Input type for two elastic parameters is not given. LAMESPARAMETERS are assumed to be given." << std::endl;
//     }
//   }

//   ierr = PetscOptionsEnd();PetscCall(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode Thermal::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols)
{
  PetscFunctionBeginUser;

  switch (this->derivedsoltype) {
    case HEATFLUX:
    SolveDerivedField(dm, &dmAux, Pw_Functions::heatflux, sol, &derivedsols[0], "heatflux");
    break;
    default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid derived solution type");
  }
  PetscFunctionReturn(0);

}

Thermal::Thermal() {
  this->material_view = PETSC_FALSE;
}

