#  CMakeLists.txt
#  Modified from $PETSC_DIR/share/petsc/CMakeLists.txt
#
#  Usage:  To build GeomechX from the source code using CMake
#    rm -fr build
#    mkdir build
#    cd build
#    cmake -S ../ -B ./
#    make all
#
#  The 'GeomechX' executable will be built in the 'build' directory.


#  By default, the compiler and library information is retrieved from PETSc
#
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set the root directory to locate PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

# Remove the lines below if you do not want to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)

# Define the list of THM source files
## Physics-related source files
file(GLOB THM_PHYSICS_SOURCES
    ./src/Physics/*.cpp
)

## PointwiseFunctions-related source files
file(GLOB THM_PF_SOURCES
    ./src/PointwiseFunctions/*.cpp
)

## Base-related source files
file(GLOB THM_BASE_SOURCES
    ./src/Base/*.cpp
)

## BoundaryCondition-related source files
file(GLOB THM_BC_SOURCES
    ./src/BoundaryCondition/*.cpp
)

## Main-related source files
file(GLOB THM_MAIN_SOURCES
    ./src/Mains/*.cpp
)

## Merge all source files
set(THM_SOURCE_FILES
    ${THM_PHYSICS_SOURCES}
    ${THM_PF_SOURCES}
    ${THM_BASE_SOURCES}
    ${THM_BC_SOURCES}
    ${THM_MAIN_SOURCES}
)


# set(THM_SOURCE_FILES
# ./src/Physics/Problem.cpp
# ./src/Physics/Elasticity.cpp
# ./src/Physics/IsotropicLinearElasticity.cpp  
# ./src/Physics/TransverselyIsotropicLinearElasticity.cpp 
# ./src/Physics/Thermal.cpp 
# ./src/Physics/SteadyThermal.cpp 
# ./src/Physics/TransientThermal.cpp 
# ./src/Physics/ThermoElasticity.cpp 
# ./src/Physics/SteadyThermoElasticity.cpp 
# ./src/Physics/TransientThermoElasticity.cpp 
# ./src/Physics/PorousFluidFlow.cpp 
# ./src/Physics/SteadyDarcysFlow.cpp 
# ./src/Physics/TransientDarcysFlow.cpp 
# ./src/Physics/PoroElasticity.cpp
# ./src/Physics/TransientPoroElasticity.cpp
# ./src/Physics/ThermoPoroElasticity.cpp
# ./src/Physics/TransientThermoPoroElasticity.cpp
# ./src/PointwiseFunctions/PF_constants.cpp 
# ./src/PointwiseFunctions/PF_Problem.cpp 
# ./src/PointwiseFunctions/PF_Elasticity.cpp 
# ./src/PointwiseFunctions/PF_Thermal.cpp 
# ./src/PointwiseFunctions/PF_ThermoElasticity.cpp 
# ./src/PointwiseFunctions/PF_PorousFluidFlow.cpp 
# ./src/PointwiseFunctions/PF_PoroElasticity.cpp
# ./src/PointwiseFunctions/PF_ThermoPoroElasticity.cpp
# ./src/Base/Mesh.cpp 
# ./src/Base/Material.cpp 
# ./src/Base/Property.cpp 
# ./src/Base/Projection.cpp
# ./src/BoundaryCondition/BoundaryCondition.cpp 
# ./src/BoundaryCondition/BCThermal.cpp 
# ./src/BoundaryCondition/BCElasticity.cpp 
# ./src/BoundaryCondition/BCPorousFluidFlow.cpp 
# ./src/Mains/main.cpp 
# ./src/Mains/IsotropicLinearElasticity_main.cpp 
# ./src/Mains/TransverselyIsotropicLinearElasticity_main.cpp 
# ./src/Mains/SteadyThermal_main.cpp 
# ./src/Mains/TransientThermal_main.cpp 
# ./src/Mains/SteadyThermoElasticity_main.cpp 
# ./src/Mains/TransientThermoElasticity_main.cpp 
# ./src/Mains/SteadyDarcysFlow_main.cpp
# ./src/Mains/TransientDarcysFlow_main.cpp
# ./src/Mains/TransientPoroElasticity_main.cpp
# ./src/Mains/TransientThermoPoroElasticity_main.cpp
# )

project(GeomechX)
add_executable(GeomechX ${THM_SOURCE_FILES})

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
set(muparser_DIR ${PETSC}/lib/cmake/muparser/)
find_package(muparser 2.0 REQUIRED)
# target_link_libraries(GeomechX PkgConfig::PETSC) # not to use muparser
target_link_libraries(GeomechX PkgConfig::PETSC muparser::muparser)


