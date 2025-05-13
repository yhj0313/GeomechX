#ifndef PROPERTY_HH
#define PROPERTY_HH

#include "petscheaders.hh"

#include <string>
// #include "Material.hh"

class Property {
    protected:
        PetscScalar value;
        PetscPointFunc pfunc; 

    private: 
        std::string unit;
        std::string description;
        std::string prop_equation;

    public:
        Property();
        Property(PetscScalar value, std::string unit, std::string description);
        Property(PetscPointFunc pfunc, std::string unit, std::string description);
        void SetProperty(PetscScalar value, std::string unit, std::string description);
        void SetProperty(PetscPointFunc pfunc, std::string unit, std::string description);
        void SetProperty(std::string prop_equation, std::string unit, std::string description);
        PetscScalar GetValue();
        std::string GetUnit();
        std::string GetDescription();
        PetscPointFunc GetFunction();
        std::string GetEquation();
        PetscErrorCode SetPropertyUser(MPI_Comm comm, std::string unit, std::string description, const char property_option_call[], PetscInt label_number);


};

#endif
