#include "PF_PoroElasticity.hh"

namespace Pw_Functions {

void f1_poroelas_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
  const PetscReal bulk_m = PetscRealPart(constants[BULK_M_HM]);  
  const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);  
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);

    }
    f1[c*dim+c] += lambda * u[uOff[2]] ;

    f1[c*dim+c] -= biot_coeff * u[uOff[1]] ; /* pore pressure */
  }
}

void g3_poroelas_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_HM]); 

  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
    }
  }
}

void g2_poroelas_up(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] = -biot_coeff;
  }
}

void g2_poroelas_ustrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;

  const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {

      g2[c*dim+c] = lambda;
  }
}

void f0_poroelas_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    const PetscReal grav = PetscRealPart(constants[GRAVITY_HM]);
    const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_HM]);

    for (c = 0, f0[0] = 0.0 ; c < dim; ++c) {
            f0[c] = u[uOff[1] + c];
    }

    f0[dim-1] -= grav * dens * perm_over_visc;


}

void f1_poroelas_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0; c < dim; ++c) {
            f1[c*dim + c] = - perm_over_visc * u[uOff[2]];
    }
}

void g0_poroelas_fluidflux_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;    

    for (c = 0 ; c < dim; ++c) {
            g0[c * dim + c] = 1.0;
    }

}

void g2_poroelas_fluidflux_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[]) {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0 ; c < dim; ++c) {
            g2[c * dim + c] = - perm_over_visc;
    }
}

void f0_poroelas_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{

    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);
    const PetscReal biot_m = PetscRealPart(constants[BIOT_MODULUS_HM]);

    f0[0] = u_t[uOff[1]] / biot_m ;

    f0[0] += biot_alpha * u_t[uOff[2]]  ; 
}

void f1_poroelas_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    const PetscReal grav = PetscRealPart(constants[GRAVITY_HM]);
    const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_HM]);

    for (c = 0; c < dim; ++c) {
            f1[c] = perm_over_visc * u_x[uOff_x[1]+c];
    }
    f1[dim - 1] -= perm_over_visc * grav * dens;
}

void g3_poroelas_pressure_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    
    PetscInt d;
    for (d = 0; d < dim; ++d) {
      g3[d * dim + d] = perm_over_visc;
    }    
}

void g1_poroelas_pressure_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]){
    PetscInt c;
    for (c = 0; c < dim; ++c) {
        g1[c*dim + c] = 1.0;
    }
    
}

void g0_poroelas_pressure_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    const PetscReal biot_m = PetscRealPart(constants[BIOT_MODULUS_HM]);

    g0[0] = u_tShift / biot_m;

}

void g0_poroelas_pressure_strain_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);
    
    g0[0] =  u_tShift * biot_alpha;

}

PetscErrorCode initial_pressure_poroelas(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
    PetscScalar *p = (PetscScalar *)ctx;

    *u = *p;

    return 0;
}

PetscErrorCode initial_pressure_t_poroelas(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
    *u = 0.0;
    return 0;
}

void f0_poroelas_strain(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;

  for (d = 0; d < dim; ++d) f0[0] += u_x[uOff_x[0]+ d * dim + d];
  f0[0] -= u[uOff[2]];
}

void g1_poroelas_strain_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) g1[d * dim + d] = 1.0;
}

void g0_poroelas_strain_strain(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = -1.0;
}

void scale_displacement_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;
const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = disp[i] ;}
}

void scale_darcyvelocity_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar darcyvelocity[])
{
PetscInt i;
const PetscScalar* pp_x = &u_x[uOff_x[1]];
const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_HM]);
const PetscReal perm = PetscRealPart(constants[PERMEABILITY_HM]);
PetscReal perm_over_visc = perm / fluid_visc;
const PetscReal grav = PetscRealPart(constants[GRAVITY_HM]);
const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_HM]);
for (i = 0; i < dim; i++) {darcyvelocity[i] = -perm_over_visc * pp_x[i] ;}
darcyvelocity[dim - 1] += perm_over_visc * grav * dens;

}

void scale_pressure_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar pressure[])
{
pressure[0] = u[uOff[1]];
}

void cauchystrain_3D_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar strain[])
{
  const PetscScalar* disp_x = &u_x[uOff_x[0]]; 
  const PetscInt Nc = dim;
  
  const PetscScalar strain_xx = disp_x[0*Nc+0];  
  const PetscScalar strain_yy = disp_x[1*Nc+1];
  const PetscScalar strain_zz = disp_x[2*Nc+2];
  const PetscScalar strain_xy = 0.5*(disp_x[0*Nc+1] + disp_x[1*Nc+0]);
  const PetscScalar strain_yz = 0.5*(disp_x[1*Nc+2] + disp_x[2*Nc+1]);
  const PetscScalar strain_xz = 0.5*(disp_x[0*Nc+2] + disp_x[2*Nc+0]);
  strain[0] = strain_xx;
  strain[1] = strain_yy;
  strain[2] = strain_zz;
  strain[3] = strain_xy;
  strain[4] = strain_yz;
  strain[5] = strain_xz;
}

void cauchystress_3D_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscScalar biot_pp = biot_alpha * u[uOff[1]] ; 

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace - biot_pp;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace - biot_pp;
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);
    const PetscScalar stress_yz = shear_m * (disp_x[1*Nc+2] + disp_x[2*Nc+1]);
    const PetscScalar stress_xz = shear_m * (disp_x[0*Nc+2] + disp_x[2*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_zz; // stress_zz
    stress[3] = stress_xy; // stress_xy
    stress[4] = stress_yz; // stress_yz
    stress[5] = stress_xz; // stress_xz
} // stress

void effectivestress_3D_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar effective_stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace;
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);
    const PetscScalar stress_yz = shear_m * (disp_x[1*Nc+2] + disp_x[2*Nc+1]);
    const PetscScalar stress_xz = shear_m * (disp_x[0*Nc+2] + disp_x[2*Nc+0]);

    effective_stress[0] = stress_xx; // stress_xx
    effective_stress[1] = stress_yy; // stress_yy
    effective_stress[2] = stress_zz; // stress_zz
    effective_stress[3] = stress_xy; // stress_xy
    effective_stress[4] = stress_yz; // stress_yz
    effective_stress[5] = stress_xz; // stress_xz
} // effective stress

void cauchystrain_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar strain[])
{
  const PetscScalar* disp_x = &u_x[uOff_x[0]]; 
  const PetscInt Nc = dim;
  
  const PetscScalar strain_xx = disp_x[0*Nc+0];  
  const PetscScalar strain_yy = disp_x[1*Nc+1];
  const PetscScalar strain_zz = 0.0;
  const PetscScalar strain_xy = 0.5*(disp_x[0*Nc+1] + disp_x[1*Nc+0]);

  strain[0] = strain_xx;
  strain[1] = strain_yy;
  strain[2] = strain_xy;
  strain[3] = strain_zz;

}

void cauchystress_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_HM]);
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscScalar biot_pp = biot_alpha * u[uOff[1]] ;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace - biot_pp;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy + 2 * biot_pp) - biot_pp; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void effectivestress_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar effective_stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_HM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) ;
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    effective_stress[0] = stress_xx; // stress_xx
    effective_stress[1] = stress_yy; // stress_yy
    effective_stress[2] = stress_xy; // stress_xy
    effective_stress[3] = stress_zz; // stress_zz
} // effective stress

void axisymmetric_2d_stress_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_HM]);
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscScalar biot_pp = biot_alpha * u[uOff[1]] ;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace - biot_pp; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace - biot_pp; 
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy + 2 * biot_pp) - biot_pp; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    //(PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sy + 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
    const PetscScalar stress_rr = PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * stress_xx + PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * stress_yy + 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * stress_xy;
    // (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sy - 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
    const PetscScalar stress_thetatheta = PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * stress_xx + PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * stress_yy - 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * stress_xy;
    // - (PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x))) * solid.sx +(PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x)))  * solid.sy +((PetscCosReal(PetscAtan2Real(y,x)))^2.0 - (PetscSinReal(PetscAtan2Real(y,x)))^2.0 ) * solid.sxy
    const PetscScalar stress_rtheta = - (PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0]))) * stress_xx +(PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0])))  * stress_yy +(PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) - PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) ) * stress_xy; 

    stress[0] = stress_rr; // stress_rr (radial stress)
    stress[1] = stress_thetatheta; // stress_thetatheta (tangential stress)
    stress[2] = stress_rtheta; // stress_rtheta
    stress[3] = stress_zz; // stress_zz

} // stress

void axisymmetric_2d_effectivestress_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar effective_stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_HM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_HM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_HM]);  

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace ; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace; 
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) ; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    //(PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sy + 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
    const PetscScalar stress_rr = PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * stress_xx + PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * stress_yy + 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * stress_xy;
    // (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sy - 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
    const PetscScalar stress_thetatheta = PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * stress_xx + PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * stress_yy - 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * stress_xy;
    // - (PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x))) * solid.sx +(PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x)))  * solid.sy +((PetscCosReal(PetscAtan2Real(y,x)))^2.0 - (PetscSinReal(PetscAtan2Real(y,x)))^2.0 ) * solid.sxy
    const PetscScalar stress_rtheta = - (PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0]))) * stress_xx +(PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0])))  * stress_yy +(PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) - PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) ) * stress_xy; 

    effective_stress[0] = stress_rr; // stress_rr (radial stress)
    effective_stress[1] = stress_thetatheta; // stress_thetatheta (tangential stress)
    effective_stress[2] = stress_rtheta; // stress_rtheta
    effective_stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_strain_planestrain_poroelas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar strain[])
{
  const PetscScalar* disp_x = &u_x[uOff_x[0]]; 
  const PetscInt Nc = dim;
  
  const PetscScalar strain_xx = disp_x[0*Nc+0];  
  const PetscScalar strain_yy = disp_x[1*Nc+1];
  const PetscScalar strain_zz = 0.0;
  const PetscScalar strain_xy = 0.5*(disp_x[0*Nc+1] + disp_x[1*Nc+0]);

  //(PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sy + 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar strain_rr = PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * strain_xx + PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * strain_yy + 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * strain_xy;
  // (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sy - 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar strain_thetatheta = PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) * strain_xx + PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) * strain_yy - 2*PetscSinReal(PetscAtan2Real(x[1],x[0])) * PetscCosReal(PetscAtan2Real(x[1],x[0])) * strain_xy;
  // - (PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x))) * solid.sx +(PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x)))  * solid.sy +((PetscCosReal(PetscAtan2Real(y,x)))^2.0 - (PetscSinReal(PetscAtan2Real(y,x)))^2.0 ) * solid.sxy
  const PetscScalar strain_rtheta = - (PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0]))) * strain_xx +(PetscCosReal(PetscAtan2Real(x[1],x[0]))) * (PetscSinReal(PetscAtan2Real(x[1],x[0])))  * strain_yy +(PetscPowReal(PetscCosReal(PetscAtan2Real(x[1],x[0])), 2) - PetscPowReal(PetscSinReal(PetscAtan2Real(x[1],x[0])), 2) ) * strain_xy; 

  strain[0] = strain_rr; // strain_rr (radial strain)
  strain[1] = strain_thetatheta; // strain_thetatheta (tangential strain)
  strain[2] = strain_rtheta; // strain_rtheta 
  strain[3] = strain_zz; // strain_zz

}

}

