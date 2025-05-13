#include "PF_ThermoElasticity.hh"

namespace Pw_Functions {

void f1_thermoelas_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]);
  const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
  const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]);  
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);
      f1[c*dim+c] += lambda*u_x[d*dim+d] ;
    }
    f1[c*dim+c] -= 3.0 * ther_expan * bulk_m * u[uOff[1]] ; /* Thermal stress */
  }
}

void g3_thermoelas_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_TM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]); 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
      g3[((c*Nc + d)*dim + c)*dim + d] += lambda;
    }
  }
}

void g2_thermoelas_ut(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
  const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]);  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 3.0 * ther_expan * bulk_m;
  }
}

//f0_temp
void f0_thermoelas_t(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_TM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_TM]); 
  const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_TM]); 

  f0[0] = fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)
//   f0[1] = density * heat_capa * flow_vel * u_x[uOff_x[1]+1];
//   f0[2] = density * heat_capa * flow_vel * u_x[uOff_x[1]+2];
}

//f1_temp
void f1_thermoelas_t(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = ther_cond * u_x[uOff_x[1] + d];
}

//g1_temp
void g1_thermoelas_tt(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_TM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_TM]); 
  const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_TM]); 
  g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
}

//g3_temp
void g3_thermoelas_tt(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond;

}

void f0_thermoelas_t_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_TM]);
  const PetscReal density = PetscRealPart(constants[DENSITY_TM]);
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_TM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_TM]); 
  const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_TM]); 

  
  f0[0] = density * heat_capa * u_t[uOff[1]+0];
  f0[0] += fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)

}

void g0_thermoelas_tt_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_TM]);
  const PetscReal density = PetscRealPart(constants[DENSITY_TM]);
  g0[0] = u_tShift * heat_capa * density;
}

void f1_thermoelas_t_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond * u_x[uOff_x[1]+d];
}

void g3_thermoelas_tt_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d)
    g3[d * dim + d] = ther_cond;
}

//mat_fields

void f1_thermoelas_u_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;

  PetscReal shear_m = a[SHEAR_M_TM]; 
  PetscReal lambda = a[LAMBDA_TM];
  PetscReal bulk_m = a[BULK_M_TM];  
  PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM]; 
  PetscReal ref_temp = a[REF_TEMP_TM]; // reference temperature for thermal stress analysis
  PetscInt  c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);
      f1[c*dim+c] += lambda*u_x[d*dim+d] ;
    }
    f1[c*dim+c] -= 3.0 * ther_expan * bulk_m * (u[uOff[1]] - ref_temp) ; /* Thermal stress */
  }

}

void g3_thermoelas_uu_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;

  PetscReal shear_m  = a[SHEAR_M_TM]; 
  PetscReal lambda = a[LAMBDA_TM]; 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
      g3[((c*Nc + d)*dim + c)*dim + d] += lambda;
    }
  }
}

void g2_thermoelas_ut_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  PetscReal bulk_m = a[BULK_M_TM];  
  PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM];  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 3.0 * ther_expan * bulk_m;
  }
}

//f0_temp
void f0_thermoelas_t_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_TM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_TM]; 
  PetscReal flow_vel = a[FLOW_VEL_TM]; 

  f0[0] = fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)
//   f0[1] = density * heat_capa * flow_vel * u_x[uOff_x[1]+1];
//   f0[2] = density * heat_capa * flow_vel * u_x[uOff_x[1]+2];
}

//f1_temp
void f1_thermoelas_t_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = ther_cond * u_x[uOff_x[1] + d];

}

//g1_temp
void g1_thermoelas_tt_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_TM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_TM]; 
  PetscReal flow_vel = a[FLOW_VEL_TM]; 
  g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
}

//g3_temp
void g3_thermoelas_tt_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond;

}

void f0_thermoelas_t_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal heat_capa = a[HEAT_CAPA_TM];
  PetscReal density = a[DENSITY_TM];
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_TM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_TM]; 
  PetscReal flow_vel = a[FLOW_VEL_TM]; 
  PetscReal heat_source = a[HEAT_SOURCE_TM];

  f0[0] = density * heat_capa * u_t[uOff[1]+0];
  f0[0] += fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)
  f0[0] -= heat_source; 

}

void g0_thermoelas_tt_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  PetscReal heat_capa = a[HEAT_CAPA_TM];
  PetscReal density = a[DENSITY_TM];
  g0[0] = u_tShift * heat_capa * density;
}

void f1_thermoelas_t_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond * u_x[uOff_x[1]+d];
}

void g3_thermoelas_tt_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d)
    g3[d * dim + d] = ther_cond;
}

//end of mat_fields

void displacement_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;
const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = disp[i] ;}
}

void temperature_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar temperature[])
{
  temperature[0] = u[uOff[1]];
}

void heatflux_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar heatflux[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;

  for (d = 0; d < dim; ++d)
    heatflux[d] = -ther_cond * u_x[uOff_x[1]+d];
}

void heatflux_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar heatflux[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;

  for (d = 0; d < dim; ++d)
    heatflux[d] = -ther_cond * u_x[uOff_x[1]+d];
}

void cauchystrain_planestrain_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_planestrain_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_TM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]); 

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) + thermal_stress; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_strain_planestrain_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void axisymmetric_2d_stress_planestrain_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_TM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]);  

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    PetscScalar stress_zz;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress; 
    if (dim == 2) {stress_zz = poisson_r * (stress_xx + stress_yy);}
    else if (dim == 3) {stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress;}
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

void cauchystrain_3D_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_3D_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]);  

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress;
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

void cauchystress_planestrain_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_TM]; 
    PetscReal lambda = a[LAMBDA_TM];
    PetscReal poisson_r = a[POISSON_R_TM];
    PetscReal bulk_m = a[BULK_M_TM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM]; 

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) + thermal_stress; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_stress_planestrain_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_TM]; 
    PetscReal lambda = a[LAMBDA_TM];
    PetscReal poisson_r = a[POISSON_R_TM];
    PetscReal bulk_m = a[BULK_M_TM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM];  

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    PetscScalar stress_zz;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress; 
    if (dim == 2) {stress_zz = poisson_r * (stress_xx + stress_yy);}
    else if (dim == 3) {stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress;}
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

void cauchystress_3D_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_TM]; 
    PetscReal lambda = a[LAMBDA_TM];
    PetscReal bulk_m = a[BULK_M_TM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM];  

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress;
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

void f1_thermoelas_u_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]);

  PetscInt        c, d;

  const PetscReal gamma = lambda/shear_m ; 

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += (u_x[c*dim+d] + u_x[d*dim+c]);
      f1[c*dim+c] += gamma*u_x[d*dim+d] ;
    }
    f1[c*dim+c] -= u[uOff[1]] ; /* Thermal stress */
  }
}

void g3_thermoelas_uu_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_TM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_TM]); 
  PetscInt        c, d;

  const PetscReal gamma = lambda/shear_m; 

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += 1;
      g3[((c*Nc + d)*dim + d)*dim + c] += 1;
      g3[((c*Nc + d)*dim + c)*dim + d] += gamma;
    }
  }
}

void g2_thermoelas_ut_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;

  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 1;
  }
}

void f0_thermoelas_t_transient_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_TM]);
  const PetscReal density = PetscRealPart(constants[DENSITY_TM]);
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_TM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_TM]); 
  const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_TM]); 
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);

  
  f0[0] = u_t[uOff[1]+0];
  f0[0] += fluid_density * fluid_heat_capa * flow_vel / ther_cond * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)

}

void f1_thermoelas_t_transient_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond  * u_x[uOff_x[1]+d];
}

void g0_thermoelas_tt_transient_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{

  g0[0] = u_tShift ;
}

void g1_thermoelas_tt_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_TM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_TM]); 
  const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_TM]); 
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);

  g1[0] = fluid_density * fluid_heat_capa * flow_vel / ther_cond ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
}

void g3_thermoelas_tt_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond ;

}

void displacement_TM_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;

const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_TM]); 
const PetscReal bulk_m = PetscRealPart(constants[BULK_M_TM]);  
const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_TM]); 
const PetscReal coeff = 3.0 * ther_expan * bulk_m / shear_m ;

const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = coeff * disp[i] ;}
}


void f1_thermoelas_u_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  PetscReal shear_m = a[SHEAR_M_TM]; 
  PetscReal lambda = a[LAMBDA_TM];
  PetscInt        c, d;
  PetscReal gamma = lambda/shear_m ; 

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += (u_x[c*dim+d] + u_x[d*dim+c]);
      f1[c*dim+c] += gamma*u_x[d*dim+d] ;
    }
    f1[c*dim+c] -= u[uOff[1]] ; /* Thermal stress */
  }
}

void g3_thermoelas_uu_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;

  PetscReal shear_m  = a[SHEAR_M_TM]; 
  PetscReal lambda = a[LAMBDA_TM]; 
  PetscInt        c, d;

  PetscReal gamma = lambda/shear_m; 

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += 1;
      g3[((c*Nc + d)*dim + d)*dim + c] += 1;
      g3[((c*Nc + d)*dim + c)*dim + d] += gamma;
    }
  }
}

void g2_thermoelas_ut_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 1;
  }
}

void f0_thermoelas_t_transient_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal heat_capa = a[HEAT_CAPA_TM];
  PetscReal density = a[DENSITY_TM];
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_TM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_TM]; 
  PetscReal flow_vel = a[FLOW_VEL_TM]; 
  PetscReal ther_cond = a[THER_COND_TM];
  f0[0] = u_t[uOff[1]+0];
  f0[0] += fluid_density * fluid_heat_capa * flow_vel / ther_cond * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)

}

void f1_thermoelas_t_transient_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond * u_x[uOff_x[1]+d];
}

void g0_thermoelas_tt_transient_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = u_tShift ;
}

void g1_thermoelas_tt_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_TM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_TM]; 
  PetscReal flow_vel = a[FLOW_VEL_TM]; 
  PetscReal ther_cond = a[THER_COND_TM];

  g1[0] = fluid_density * fluid_heat_capa * flow_vel / ther_cond ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
}

void g3_thermoelas_tt_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscReal ther_cond = a[THER_COND_TM];
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond ;

}

void displacement_TM_dimensionless_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;

PetscReal shear_m = a[SHEAR_M_TM]; 
PetscReal bulk_m = a[BULK_M_TM];  
PetscReal ther_expan = a[LIN_THER_EXPAN_C_TM]; 
PetscReal coeff = 3.0 * ther_expan * bulk_m / shear_m ;

const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = coeff * disp[i] ;}
}

// 아래 두 함수에 void *ctx 추가함. 
void property_user_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar property[])
{
  PetscScalar u1 = u[uOff[0]];
  PetscScalar u2= u[uOff[0]+1];
  PetscScalar u3;
  if (dim == 3) u3 = u[uOff[0]+2];
  PetscScalar T = u[uOff[1]];  // temperature in degC
  PetscScalar x1= x[0];
  PetscScalar x2= x[1];
  PetscScalar x3;
  if (dim == 3) x3 = x[2];

  std::string *p = (std::string *)ctx;

  const std::string pro_eq = *p;
  try {
      mu::Parser parser;

      parser.SetExpr(pro_eq);
      parser.DefineVar("T", &T);
      parser.DefineVar("u1", &u1);
      parser.DefineVar("u2", &u2);
      if (dim ==3) parser.DefineVar("u3", &u3);
      parser.DefineVar("x", &x1);
      parser.DefineVar("y", &x2);
      if (dim ==3) parser.DefineVar("z", &x3);
      // parser.Compile();
      property[0] = parser.Eval();
  } catch (mu::Parser::exception_type& e) {
      std::cout << e.GetMsg() << std::endl;
  }


}

}

