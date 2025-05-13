#ifndef MATERIAL_HH
#define MATERIAL_HH

#include "petscheaders.hh"

#include "Property.hh"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

class Material{
    private: 
        std::string name;
        PetscInt label_number;
        std::string description;
        
    protected:   
        std::vector <Property> properties;
        
        
    public:    
        Material();
        Material(std::string name, PetscInt label_number, std::string description, std::vector <Property> properties);
        Material(std::string name, PetscInt label_number, std::string description);
        void SetName(std::string name);
        void SetLabel_number(PetscInt label_number);
        void SetDescription(std::string description);
        void AddProperty(Property& prop);
        std::string GetName();
        PetscInt GetLabel_number();
        std::string GetDescription();
        PetscScalar GetProperty(PetscInt index);
        PetscInt GetNumofProperties();
        std::string GetPropertyUnit(PetscInt index);
        std::string GetPropertyDescription(PetscInt index);
        PetscPointFunc GetPropertyFunction(PetscInt index);
        std::string GetPropertyEquation(PetscInt index);
        PetscErrorCode ViewMaterial();
        
};

#endif