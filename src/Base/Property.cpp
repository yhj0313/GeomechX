#include "Property.hh"

Property::Property() 
    : value(0.0), unit(""), description(""), pfunc(NULL), prop_equation("") {} 

Property::Property(PetscScalar value, std::string unit, std::string description) {
    this->value = value;
    this->unit = unit;
    this->description = description;
    this->pfunc = NULL;
    this->prop_equation = "";
}

Property::Property(PetscPointFunc pfunc, std::string unit, std::string description) {
    this->value = 0.0;
    this->unit = unit;
    this->description = description;
    this->pfunc = pfunc;
    this->prop_equation = "";
}

void Property::SetProperty(PetscScalar value, std::string unit, std::string description) {
    this->value = value;
    this->unit = unit;
    this->description = description;
    this->pfunc = NULL;
    this->prop_equation = "";
}

void Property::SetProperty(PetscPointFunc pfunc, std::string unit, std::string description) {
    this->value = 0.0;
    this->unit = unit;
    this->description = description;
    this->pfunc = pfunc;
    this->prop_equation = "";
}

void Property::SetProperty(std::string prop_equation, std::string unit, std::string description) {
    this->value = 0.0;
    this->unit = unit;
    this->description = description;
    this->pfunc = NULL;
    this->prop_equation = prop_equation;
}

PetscScalar Property::GetValue() {
    return this->value;
}

std::string Property::GetUnit() {
    return this->unit;
}

std::string Property::GetDescription(){
    return this->description;
}

PetscPointFunc Property::GetFunction() {
    return this->pfunc;
}

std::string Property::GetEquation() {
    return this->prop_equation;
}

PetscErrorCode Property::SetPropertyUser(MPI_Comm comm, std::string unit, std::string description, const char property_option_call[], PetscInt label_number) {
    
    char option_call[PETSC_MAX_PATH_LEN];
    char prop_equation[PETSC_MAX_PATH_LEN]; 
    PetscBool flg_prop_equation, flg_prop_value;
    PetscScalar prop_value_temp(0.0);

    PetscFunctionBegin; 
    
    PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_equation_%d", property_option_call, label_number));
    PetscOptionsGetString(NULL, NULL, option_call, prop_equation, sizeof(prop_equation), &flg_prop_equation);
    
    if (flg_prop_equation) {
      // PetscPointFunc lin_ther_expan_c_func = Pw_Functions::lin_thermal_expan_c;
      // lin_ther_expan_c.SetProperty(Pw_Functions::lin_thermal_expan_c, "1/K", "Linear thermal expansion coefficient" );
      this->SetProperty(prop_equation, unit, description);
    }
    else {
      PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-%s_%d", property_option_call, label_number));
      PetscOptionsGetScalar(NULL, NULL, option_call, &prop_value_temp, &flg_prop_value);
      if (!flg_prop_value) PetscCall(PetscPrintf(comm, "A value for %s is not given. It is assumed to be 0.0.\n", const_cast<char*>(description.c_str())));
      this->SetProperty(prop_value_temp, unit, description); 
    }

    PetscFunctionReturn(0); 
}