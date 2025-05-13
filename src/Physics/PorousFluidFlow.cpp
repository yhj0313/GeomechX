#include "PorousFluidFlow.hh"

PetscErrorCode PorousFluidFlow::ProcessOptions(MPI_Comm comm) {
    
  PetscFunctionBeginUser;
  PetscCall(PetscStrncpy(this->dmType, DMPLEX, 256));

  PetscOptionsBegin(comm, "", "Problem Options", "");
  PetscCall(PetscOptionsBool("-gravity", "Use gravity in the y direction in 2D and in the z direction in 3D", NULL, this->useGravity, &this->useGravity, NULL));

  PetscCall(PetscOptionsBool("-material_view", "View material and its properties", NULL, this->material_view, &this->material_view, NULL));
  PetscCall(PetscOptionsFList("-dm_type", "Convert DMPlex to another format", NULL, DMList, this->dmType, this->dmType, 256, NULL));

  PetscOptionsEnd();
  PetscFunctionReturn(0);

}

PetscErrorCode PorousFluidFlow::ParameterSetup(PetscDS ds)
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


PetscErrorCode PorousFluidFlow::InputParameterOptions(MPI_Comm comm, const std::string name)
{
//   PetscInt       para = 0; 
//   PetscErrorCode ierr;
//   std::string option_name = "-";
//   PetscBool      flg;
//   PetscMPIInt    rank;

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

PetscErrorCode PorousFluidFlow::CreateDerivedField(DM dm, DM *dmAux, const char name[], PetscFE *fe, PetscFE *fe_derived, PetscInt numComp, PetscInt dim) {
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

PetscErrorCode PorousFluidFlow::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[]) {
  char           view_call[PETSC_MAX_PATH_LEN]; 
  // add 
  PetscInt       n;

  PetscFunctionBegin;
  
  PetscCall(DMCreateGlobalVector(*dmAux, derivedsol));
  PetscCall(VecSet(*derivedsol, 0.0));
  PetscCall(PetscObjectSetName((PetscObject) *derivedsol, name));
  PetscCall(DMProjectField(*dmAux, 0.0, *sol, &func, INSERT_VALUES, *derivedsol));

  PetscCall(PetscSNPrintf(view_call, PETSC_MAX_PATH_LEN, "-%s_view", name));
  PetscCall(VecViewFromOptions(*derivedsol, NULL, view_call));  


  PetscFunctionReturn(0);  
}

PetscErrorCode PorousFluidFlow::SolveDerivedField(DM dm, DM *dmAux, PetscPointFunc func, Vec *sol, Vec* derivedsol, const char name[], PetscReal time, PetscInt step) {
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

PetscErrorCode PorousFluidFlow::SolveDerivedFields(DM dm, DM &dmAux, Vec *sol, Vec *derivedsols) 
{
  PetscInt       n;
  
  PetscFunctionBeginUser;
  SolveDerivedField(dm, &dmAux, Pw_Functions::scale_darcyvelocity, sol, &derivedsols[0], "darcyvelocity");
  SolveDerivedField(dm, &dmAux, Pw_Functions::scale_pressure, sol, &derivedsols[1], "pressure");

  PetscFunctionReturn(0);

}

PorousFluidFlow::PorousFluidFlow() {
  this->material_view = PETSC_FALSE;
  this->useGravity = PETSC_FALSE;
}

