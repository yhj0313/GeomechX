#include "PF_PorousFluidFlow.hh"
const char *bdvalue_call_H[6] = {"-bdvalue_darcy_left", "-bdvalue_darcy_right", "-bdvalue_darcy_bottom", "-bdvalue_darcy_top", "-bdvalue_darcy_front", "-bdvalue_darcy_back"};
const char *bdvalues_call_H[6] = {"-bdvalues_darcy_left", "-bdvalues_darcy_right", "-bdvalues_darcy_bottom", "-bdvalues_darcy_top", "-bdvalues_darcy_front", "-bdvalues_darcy_back"};

void Getbdvalue_H(PetscScalar &boundaryvalue, PetscInt bdlocation)
{
  PetscScalar bdv = 0.0;
  std::cout << "bd value call: " <<  bdvalue_call_H[bdlocation] << std::endl;
  PetscOptionsGetScalar(NULL, NULL, bdvalue_call_H[bdlocation], &bdv, NULL);
//   boundaryvalue = bdv * 1.0e-6; // Pore pressure
//   boundaryvalue = bdv * 1.0e-9; // Darcy velocity
  boundaryvalue = bdv;
}

void inline Getbdvalue_pp_HM(PetscScalar &boundaryvalue, PetscInt bdlocation)
{
  PetscScalar bdv = 0.0;
  PetscOptionsGetScalar(NULL, NULL, bdvalue_call_H[bdlocation], &bdv, NULL);
  boundaryvalue = bdv;
}

namespace Pw_Functions
{
        
    void f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {

        PetscInt c;
        for (c = 0 ; c < dim; ++c) {
                f0[0] += u_x[c * dim + c];
        }
    }

    void f0_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        const PetscReal compr = PetscRealPart(constants[COMPRESSIBILITY_H]);
        const PetscReal poro = PetscRealPart(constants[POROSITY_H]);

        PetscInt c;
        for (c = 0 ; c < dim; ++c) {
                f0[0] += u_x[c * dim + c];
        }
        f0[0] += compr * poro * u_t[uOff[1]];
    }

   void f0_pressure_transient_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        const PetscReal compr = PetscRealPart(constants[COMPRESSIBILITY_H]);
        const PetscReal poro = PetscRealPart(constants[POROSITY_H]);

        f0[0] = compr * poro * u_t[uOff[0]];
    }

    void f1_pressure_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        const PetscReal grav = PetscRealPart(constants[GRAVITY_H]);
        const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_H]);

        for (c = 0; c < dim; ++c) {
                f1[c] = perm_over_visc * u_x[c];
        }
        f1[dim - 1] -= perm_over_visc * grav * dens;
    }

    void g0_pressure_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
    {
        PetscInt c;   
        const PetscReal compr = PetscRealPart(constants[COMPRESSIBILITY_H]);
        const PetscReal poro = PetscRealPart(constants[POROSITY_H]); 
  
        for (c = 0 ; c < dim; ++c) {
                g0[c] = compr * poro * u_tShift;
        }

    }

   void g0_pressure_pressure_transient_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
    {
        PetscInt c;   
        const PetscReal compr = PetscRealPart(constants[COMPRESSIBILITY_H]);
        const PetscReal poro = PetscRealPart(constants[POROSITY_H]); 
  
        for (c = 0 ; c < dim; ++c) {
                g0[c] = compr * poro * u_tShift;
        }

    }

    void g3_pressure_pressure_transient_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
        const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
        const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
        PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
        {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        
        PetscInt d;
        for (d = 0; d < dim; ++d)
        g3[d * dim + d] = perm_over_visc;
    }

  void f0_darcy_velocity_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 0); // value on the left boundary

    f0[0] = -q0; // An inflow q0 is positive.
  }

  void f0_darcy_velocity_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 1); // value on the right boundary

    f0[0] = -q0; // An inflow q0 is positive.
  }

  void f0_darcy_velocity_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 2); // value on the bottom boundary    

    f0[0] = -q0; // An inflow q0 is positive.
  }

  void f0_darcy_velocity_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 3); // value on the top boundary        
 
    f0[0] = -q0;  // An inflow q0 is positive.
  }

  void f0_darcy_velocity_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 4); // value on the front boundary    

    f0[0] = -q0; // An inflow q0 is positive.
  }

  void f0_darcy_velocity_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue_H(q0, 5); // value on the back boundary        

    f0[0] = -q0; // An inflow q0 is positive.
  }


    void g1_pressure_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]){
        PetscInt c;
        for (c = 0; c < dim; ++c) {
            g1[c*dim + c] = 1.0;
        }
        
    }

    void f0_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        const PetscReal grav = PetscRealPart(constants[GRAVITY_H]);
        const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_H]);
  
        for (c = 0, f0[0] = 0.0 ; c < dim; ++c) {
                f0[c] = u[c];
        }

        f0[dim-1] -= grav * dens * perm_over_visc;
    }

    
    void f1_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;

        for (c = 0; c < dim; ++c) {
                f1[c*dim + c] = - perm_over_visc * u[uOff[1]];
        }
    }

    void g0_fluidflux_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
    {
        PetscInt c;    
  
        for (c = 0 ; c < dim; ++c) {
                g0[c * dim + c] = 1.0;
        }

    }
    
    void g2_fluidflux_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[]) {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
  
        for (c = 0 ; c < dim; ++c) {
                g2[c * dim + c] = - perm_over_visc;
        }
    }

    void f0_noflow_pressure_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        for (c = 0 ; c < dim; ++c) {
                f0[0] += u[c] * n[c];
        }
    }

    void g0_noflow_pressureflux_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
    {
        PetscInt c;
        for (c = 0 ; c < dim; ++c) {
                g0[c] = n[c];
        }
    }

    PetscErrorCode boundary_value(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
    {
        BoundaryCondition *p = (BoundaryCondition *)ctx;

        *u = p->GetBdValue();

        return 0;
    }

    PetscErrorCode boundary_values(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
    {
        BoundaryCondition *p = (BoundaryCondition *)ctx;

        *u = *p->GetBdValues();

        return 0;
    }

    PetscErrorCode zero(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
    {
       PetscInt c;
        for (c = 0; c < Nc; ++c) u[c] = 0.0;
        return 0;
    }

    void f0_pressureboundary_fluidflux_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
  
        for (c = 0; c < dim; ++c) {
            f0[c] +=  perm_over_visc * x[0] * n[c];

        }
    }

    void f0_pressureboundary_fluidflux_bd_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;
        Getbdvalue_H(P, 0);
        for (c = 0; c < dim; ++c) {
            f0[c] +=  perm_over_visc * P * n[c];

        }


    }
    void f0_pressureboundary_fluidflux_bd_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;

        Getbdvalue_H(P, 1);
        for (c = 0; c < dim; ++c) {
                f0[c] +=  perm_over_visc * P * n[c];

        }
    }
    void f0_pressureboundary_fluidflux_bd_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;

        Getbdvalue_H(P, 2);
        for (c = 0; c < dim; ++c) {
                f0[c] +=  perm_over_visc * P * n[c];
        }

    }
    void f0_pressureboundary_fluidflux_bd_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;
        Getbdvalue_H(P, 3);
        for (c = 0; c < dim; ++c) {
                f0[c] +=  perm_over_visc * P * n[c];

        }

    }
    void f0_pressureboundary_fluidflux_bd_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;

        Getbdvalue_H(P, 4);
        for (c = 0; c < dim; ++c) {
                f0[c] +=  perm_over_visc * P * n[c];
        }

    }
    void f0_pressureboundary_fluidflux_bd_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        PetscReal P;

        Getbdvalue_H(P, 5);
        for (c = 0; c < dim; ++c) {
                f0[c] +=  perm_over_visc * P * n[c];

        }

    }

    void f0_massflux_pressure_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_H]);
        for (c = 0 ; c < dim; ++c) {
                f0[0] += dens * u[c] * n[c];
        }

        
    }

    void f0_massflux_fluidflux_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
        const PetscReal grav = PetscRealPart(constants[GRAVITY_H]);
        const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_H]);
  
        for (c = 0, f0[0] = 0.0 ; c < dim; ++c) {
                f0[c] += n[c] * 1.0;
        }

        f0[dim-1] += grav * dens * perm_over_visc;
        
    }

    void f1_massflux_fluidflux_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;

        for (c = 0; c < dim; ++c) {
                f1[c*dim + c] = - perm_over_visc * u[uOff[1]];
        }
    }

    void g2_massflux_fluxpressure_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
    {
        PetscInt c;
        const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_H]);
        const PetscReal perm = PetscRealPart(constants[PERMEABILITY_H]);
        PetscReal perm_over_visc = perm / fluid_visc;
  
        for (c = 0 ; c < dim; ++c) {
                g2[c * dim + c] = - perm_over_visc;
        }
    }
    

    void g0_massflux_pressureflux_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
    {
        PetscInt c;
        const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_H]);
        for (c = 0 ; c < dim; ++c) {
                g0[c] = dens * n[c];
        }
    }

    void scale_darcyvelocity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar darcyvelocity[])
    {
    PetscInt i;
    const PetscScalar* darcyv = &u[uOff[0]];
    for (i = 0; i < dim; i++) {darcyvelocity[i] = darcyv[i] * 1.0e-9;}
    }

    void scale_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar pressure[])
    {
    pressure[0] = u[uOff[0]] ;
    }

    PetscErrorCode initial_pressure(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
    {
        PetscScalar *p = (PetscScalar *)ctx;

        *u = *p;

        return 0;
    }

    PetscErrorCode initial_pressure_t(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
    {
        *u = 0.0;
        return 0;
    }

    void f0_pressureboundary_fluidflux_HM_bd_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;
    Getbdvalue_pp_HM(P, 0);
    for (c = 0; c < dim; ++c) {
        f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_HM_bd_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;
    Getbdvalue_pp_HM(P, 1);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_HM_bd_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 2);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_HM_bd_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 3);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_HM_bd_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 4);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_HM_bd_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 5);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }    

    
    void f0_pressureboundary_fluidflux_THM_bd_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 0);
    for (c = 0; c < dim; ++c) {
        f0[c] =  perm_over_visc * P * n[c];

    }
    }
    void f0_pressureboundary_fluidflux_THM_bd_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 1);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 2);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 3);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 4);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 5);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }  

    void f0_pressureboundary_fluidflux_THM_bd_0_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 0);
    for (c = 0; c < dim; ++c) {
        f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_1_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 1);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_2_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 2);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }
    }
    void f0_pressureboundary_fluidflux_THM_bd_3_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 3);
    std::cout << P << std::endl; 
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_4_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 4);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];

    }

    }
    void f0_pressureboundary_fluidflux_THM_bd_5_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
    {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal P;

    Getbdvalue_pp_HM(P, 5);
    for (c = 0; c < dim; ++c) {
            f0[c] =  perm_over_visc * P * n[c];
    }

    }    
  

}

