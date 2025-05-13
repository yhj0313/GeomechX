#include "Material.hh"


Material::Material() 
    : name(""), label_number(0), description(""), properties() {}

Material::Material(std::string name, PetscInt label_number, std::string description, std::vector <Property> properties) {
    this->name = name;
    this->label_number = label_number;
    this->description = description;
    this->properties = properties;
}

Material::Material(std::string name, PetscInt label_number, std::string description) {
    this->name = name;
    this->label_number = label_number;
    this->description = description;
}

void Material::SetName(std::string name) {
    this->name = name;
}

void Material::SetLabel_number(PetscInt label_number){
    this->label_number = label_number;
}

void Material::SetDescription(std::string description){
    this->description = description;
}

void Material::AddProperty(Property& prop){
    this->properties.push_back(prop);

}

std::string Material::GetName() {
    return this->name;
}

PetscInt Material::GetLabel_number() {
    return this->label_number;
};

std::string Material::GetDescription() {
    return this->description;
};

PetscScalar Material::GetProperty(PetscInt index){
    return this->properties[index].GetValue();
}

PetscInt Material::GetNumofProperties(){
    return this->properties.size();
}

std::string Material::GetPropertyUnit(PetscInt index){
    return this->properties[index].GetUnit();
}

std::string Material::GetPropertyDescription(PetscInt index){
    return this->properties[index].GetDescription();
}

PetscPointFunc Material::GetPropertyFunction(PetscInt index){
    return this->properties[index].GetFunction();
}

std::string Material::GetPropertyEquation(PetscInt index){
    return this->properties[index].GetEquation();
}

PetscErrorCode Material::ViewMaterial()
{
  const PetscInt num = GetNumofProperties();
  PetscFunctionBeginUser;

  std::ofstream fout;
  fout.open("Material.txt");

  if (fout)
  {
      fout << "Material Name: " << GetName() << std::endl;
      fout << "Label Number: " << GetLabel_number() << std::endl;
      fout << "Description: " << GetDescription() << std::endl;
  }
  else {
      std::cout << "cannot open the file to write Material." << std::endl;
  }

  for (int i = 0; i < num; i++)
  {
      fout << GetPropertyDescription(i) << " = " << GetProperty(i) << " " << GetPropertyUnit(i)<< std::endl;
  }

  fout.close();

  PetscFunctionReturn(0);
}
