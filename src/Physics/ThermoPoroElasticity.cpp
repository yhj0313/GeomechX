#include "ThermoPoroElasticity.hh"

PetscErrorCode ThermoPoroElasticity::CreateElasticityNullSpace(DM dm, PetscInt origField, PetscInt field, MatNullSpace *nullspace)
{
  PetscFunctionBegin;
  PetscCall(DMPlexCreateRigidBody(dm, origField, nullspace));
  PetscFunctionReturn(0);
}

void ThermoPoroElasticity::AddMaterial(Material& mtrl){
    this->materials.push_back(mtrl);
}

PetscInt ThermoPoroElasticity::GetNumofMaterials(){
    return this->materials.size();
}

PetscErrorCode ThermoPoroElasticity::SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step)
{
  DM       dmAux, coordDM;
  PetscInt f;
  Vec      *vec_sol;

  PetscFunctionBegin;
  /* MUST call DMGetCoordinateDM() in order to get p4est setup if present */
  PetscCall(DMGetCoordinateDM(dm, &coordDM));
  if (!feAux) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DMClone(dm, &dmAux));
  PetscCall(DMSetCoordinateDM(dmAux, coordDM));
  for (f = 0; f < NfAux; ++f) PetscCall(DMSetField(dmAux, f, NULL, (PetscObject)feAux[f]));
  PetscCall(DMCreateDS(dmAux));
  PetscCall(DMGetGlobalVector(dm, vec_sol));
  PetscCall(SetupMaterial(dm, dmAux, time, step, vec_sol));
  PetscCall(DMDestroy(&dmAux));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ThermoPoroElasticity::SetupAuxDM(DM dm, PetscInt NfAux, PetscFE feAux[], PetscScalar time, PetscInt step, Vec *sol)
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
  PetscCall(SetupMaterial(dm, dmAux, time, step, sol));
  PetscCall(DMDestroy(&dmAux));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ThermoPoroElasticity::SetupMaterial(DM dm, DM dmAux, PetscScalar time, PetscInt step, Vec *sol)
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
      if (j == HEAT_SOURCE_THM) // heat source 17
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
          pointfuncs_user[j] = Pw_Functions::property_user_THM;
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
          pointfuncs_user[j] = Pw_Functions::property_user_THM;
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

PetscErrorCode ThermoPoroElasticity::CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim)
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

PetscErrorCode ThermoPoroElasticity::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[])
{ 
  char           view_call[PETSC_MAX_PATH_LEN]; 
  // add 
  PetscInt       n;

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, 0.0, *sol, &func, INSERT_VALUES, *derivedsol));
  PetscCall(DMSetAuxiliaryVec(dm, NULL, 0, 0, *derivedsol));
  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  

  // add
  PetscCall(DMGetNumAuxiliaryVec(dm, &n));
  std::cout << n << std::endl;
  

  PetscFunctionReturn(0);  
}

PetscErrorCode ThermoPoroElasticity::ParameterSetup(PetscDS ds)
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

PetscErrorCode ThermoPoroElasticity::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) 
{
  PetscFunctionBeginUser;
  PetscFunctionReturn(0);
}

ThermoPoroElasticity::ThermoPoroElasticity() {
  this->material_view = PETSC_FALSE;
}
  
