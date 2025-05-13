#include "PF_ThermoPoroElasticity.hh"

namespace Pw_Functions {
 
void f1_thermoporoelas_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);
  const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
  const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]);  
  const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);
      // f1[c*dim+c] += lambda*u_x[d*dim+d] ;
    }
    f1[c*dim+c] += lambda * u[uOff[3]] ;
    f1[c*dim+c] -= biot_coeff * u[uOff[2]] ; /* pore pressure */
    f1[c*dim+c] -= 3.0 * ther_expan * bulk_m * u[uOff[1]] ; /* Thermal stress */
  }
}

void g3_thermoporoelas_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_THM]); 
  // const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]); 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
      // g3[((c*Nc + d)*dim + c)*dim + d] += lambda;
    }
  }
}

void g2_thermoporoelas_ut(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
  const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]);  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 3.0 * ther_expan * bulk_m;
  }
}

void g2_thermoporoelas_up(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
  const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      // g2[c*dim+c] += biot_coeff;
      g2[c*dim+c] -= biot_coeff;
  }
}

void g2_thermoporoelas_ustrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal bulk_m = PetscRealPart(constants[BULK_M_HM]);  
  // const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);  
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      // g2[c*dim+c] += biot_coeff;
      g2[c*dim+c] = lambda;
  }
}

void f1_thermoporoelas_t(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_THM]);
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = ther_cond * u_x[uOff_x[1] + d];

}

void g3_thermoporoelas_tt(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_THM]);
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond;
  
}

void f0_thermoporoelas_t_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_THM]);
  const PetscReal density = PetscRealPart(constants[DENSITY_THM]);
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_THM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_THM]); 
  const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
  const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
  PetscReal perm_over_visc = perm / fluid_visc;
  // const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_THM]); 
  PetscInt  d;
  PetscReal flow_vel_delta_temp(0.0); /* inner product of flow velocity and delta temperature  */
  
  f0[0] = density * heat_capa * u_t[uOff[1]+0];
  // for (d = 0; d < dim; ++d){
  //   // flow_vel_delta_temp += u[uOff[2]+d] * u_x[uOff_x[1]+d];
  //   flow_vel_delta_temp += - perm_over_visc * u_x[uOff_x[2]+d] * u_x[uOff_x[1]+d];
  // }
  // f0[0] += fluid_density * fluid_heat_capa * flow_vel_delta_temp; //convection (x, y, z)
  // f0[0] += fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)

}

void g0_thermoporoelas_tt_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_THM]);
  const PetscReal density = PetscRealPart(constants[DENSITY_THM]);
  g0[0] = u_tShift * heat_capa * density;
}

void g0_thermoporoelas_tv_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_THM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_THM]); 
  PetscInt  d;
  
  for (d = 0; d < dim; ++d) {
    g0[d] = fluid_density * fluid_heat_capa * u_x[uOff_x[1]+d];
  }
}

void g1_thermoporoelas_tt_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_THM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_THM]); 
  const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
  const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
  PetscReal perm_over_visc = perm / fluid_visc;
  // const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_THM]); 
  PetscInt  d;
  // g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
  for (d = 0; d < dim; ++d) {
    // g1[d] = fluid_density * fluid_heat_capa * u[uOff[2]+d];
    g1[d] = - fluid_density * fluid_heat_capa * perm_over_visc * u_x[uOff_x[2]+d]; 
  }

}

void g1_thermoporoelas_tp_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  const PetscReal fluid_heat_capa = PetscRealPart(constants[FLUID_HEAT_CAPA_THM]); 
  const PetscReal fluid_density = PetscRealPart(constants[FLUID_DENSITY_THM]); 
  const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
  const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
  PetscReal perm_over_visc = perm / fluid_visc;
  // const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_THM]); 
  PetscInt  d;
  // g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
  for (d = 0; d < dim; ++d) {
    // g1[d] = fluid_density * fluid_heat_capa * u[uOff[2]+d];
    g1[d] = - fluid_density * fluid_heat_capa * perm_over_visc * u_x[uOff_x[1]+d] ;
  }

}


void f1_thermoporoelas_t_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_THM]);
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond * u_x[uOff_x[1]+d];
}

void f1_thermoporoelas_u_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  PetscReal shear_m = a[SHEAR_M_THM]; 
  PetscReal lambda = a[LAMBDA_THM];
  PetscReal bulk_m = a[BULK_M_THM];  
  PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM]; 
  PetscReal biot_coeff = a[BIOT_COEFFICIENT_THM];
  PetscReal ref_temp = a[REF_TEMP_THM]; // reference temperature for thermal stress analysis
  PetscReal ref_pres = a[REF_PRES_THM]; // reference temperature for thermal stress analysis
  PetscInt  c, d;
  // std::cout << "shear_m= " << shear_m << std::endl;
  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);
      // f1[c*dim+c] += lambda*u_x[d*dim+d] ;
    }
    f1[c*dim+c] += lambda * u[uOff[3]] ;
    // f1[c*dim+c] -= biot_coeff * u[uOff[2]] ; /* pore pressure */
    f1[c*dim+c] -= biot_coeff * (u[uOff[2]] - ref_pres) ; /* pore pressure */
    f1[c*dim+c] -= 3.0 * ther_expan * bulk_m * (u[uOff[1]] - ref_temp) ; /* Thermal stress */
  }

}

void g3_thermoporoelas_uu_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_THM]); 
  // const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]); 
  PetscReal shear_m  = a[SHEAR_M_THM]; 
  // PetscReal lambda = a[LAMBDA_THM]; 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
      // g3[((c*Nc + d)*dim + c)*dim + d] += lambda;
    }
  }
}

void g2_thermoporoelas_ut_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  PetscReal bulk_m = a[BULK_M_THM];  
  PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM];  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 3.0 * ther_expan * bulk_m;
  }
}

void g2_thermoporoelas_up_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
  PetscReal biot_coeff = a[BIOT_COEFFICIENT_THM];  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      // g2[c*dim+c] += biot_coeff;
      g2[c*dim+c] -= biot_coeff;
  }
}

void g2_thermoporoelas_ustrain_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal bulk_m = PetscRealPart(constants[BULK_M_HM]);  
  // const PetscReal biot_coeff = PetscRealPart(constants[BIOT_COEFFICIENT_HM]);  
  const PetscReal lambda = a[LAMBDA_THM];
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      // g2[c*dim+c] += biot_coeff;
      g2[c*dim+c] = lambda;
  }
}

//f1_temp
void f1_thermoporoelas_t_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscReal ther_cond = a[THER_COND_THM];
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = ther_cond * u_x[uOff_x[1] + d];

}

//g3_temp
void g3_thermoporoelas_tt_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscReal ther_cond = a[THER_COND_THM];
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = ther_cond;

}

void f0_thermoporoelas_t_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal heat_capa = a[HEAT_CAPA_THM];
  PetscReal density = a[DENSITY_THM];
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_THM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_THM]; 
  PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
  PetscReal perm = a[PERMEABILITY_THM];
  PetscReal perm_over_visc = perm / fluid_visc;
  // PetscReal flow_vel = a[FLOW_VEL_THM]; 
  PetscReal heat_source = a[HEAT_SOURCE_THM];
  PetscInt  d;
  PetscReal flow_vel_delta_temp(0.0); /* inner product of flow velocity and delta temperature  */
  
  f0[0] = density * heat_capa * u_t[uOff[1]+0];
  // for (d = 0; d < dim; ++d){
  //   // flow_vel_delta_temp += u[uOff[2]+d] * u_x[uOff_x[1]+d];
  //   flow_vel_delta_temp += - perm_over_visc * u_x[uOff_x[2]+d] * u_x[uOff_x[1]+d];
  // }
  // // f0[0] += fluid_density * fluid_heat_capa * flow_vel * u_x[uOff_x[1]+0]; //convection (일단 x방향만 고려)
  // f0[0] += fluid_density * fluid_heat_capa * flow_vel_delta_temp; 
  f0[0] -= heat_source; 

}

void g0_thermoporoelas_tt_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  PetscReal heat_capa = a[HEAT_CAPA_THM];
  PetscReal density = a[DENSITY_THM];
  g0[0] = u_tShift * heat_capa * density;
}

void g0_thermoporoelas_tv_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_THM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_THM]; 
  PetscInt  d;
  
  for (d = 0; d < dim; ++d) {
    g0[d] = fluid_density * fluid_heat_capa * u_x[uOff_x[1]+d];
  }
}

void g1_thermoporoelas_tt_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_THM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_THM]; 
  PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
  PetscReal perm = a[PERMEABILITY_THM];
  PetscReal perm_over_visc = perm / fluid_visc;
  // const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_THM]); 
  PetscInt  d;
  // g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
  for (d = 0; d < dim; ++d) {
    // g1[d] = fluid_density * fluid_heat_capa * u[uOff[2]+d];
    g1[d] = - fluid_density * fluid_heat_capa * perm_over_visc * u_x[uOff_x[2]+d] ;
  }
}


void f1_thermoporoelas_t_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscReal ther_cond = a[THER_COND_THM];
  PetscInt d;
  for (d = 0; d < dim; ++d)
    f1[d] = ther_cond * u_x[uOff_x[1]+d];
}

void g1_thermoporoelas_tp_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{// convection (일단 x방향만 고려)
 
  PetscReal fluid_heat_capa = a[FLUID_HEAT_CAPA_THM]; 
  PetscReal fluid_density = a[FLUID_DENSITY_THM]; 
  PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
  PetscReal perm = a[PERMEABILITY_THM];
  PetscReal perm_over_visc = perm / fluid_visc;
  // const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_THM]); 
  PetscInt  d;
  // g1[0] = fluid_density * fluid_heat_capa * flow_vel ; // convection (일단 x방향만 고려)
//   g1[1] = density * heat_capa * flow_vel;
//   g1[2] = density * heat_capa * flow_vel;
  for (d = 0; d < dim; ++d) {
    // g1[d] = fluid_density * fluid_heat_capa * u[uOff[2]+d];
    g1[d] = - fluid_density * fluid_heat_capa * perm_over_visc * u_x[uOff_x[1]+d] ;
  }

}

//end of mat_fields

void displacement_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;
const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = disp[i] ;}
}

void temperature_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar temperature[])
{
  temperature[0] = u[uOff[1]];
}

void porepressure_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar pore_pressure[])
{
  pore_pressure[0] = u[uOff[2]];
}

void heatflux_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar heatflux[])
{
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_THM]);
  PetscInt d;

  for (d = 0; d < dim; ++d)
    heatflux[d] = -ther_cond * u_x[uOff_x[1]+d];
}

void heatflux_THM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
              const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
              const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
              PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar heatflux[])
{
  PetscReal ther_cond = a[THER_COND_THM];
  PetscInt d;

  for (d = 0; d < dim; ++d)
    heatflux[d] = -ther_cond * u_x[uOff_x[1]+d];
}

void cauchystrain_planestrain_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_planestrain_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_THM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]); 
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) + thermal_stress - biot_pp; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_strain_planestrain_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void axisymmetric_2d_stress_planestrain_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_THM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]);  
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);


    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ;
    PetscScalar stress_zz;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp; // + 2844975.0;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp; // + 2844975.0;
    // const PetscScalar stress_zz = poisson_r * (stress_xx +  2844975.0 + stress_yy + 2844975.0);
    if (dim == 2) {stress_zz = poisson_r * (stress_xx + stress_yy + 2 * biot_pp) - biot_pp; }
    else if (dim == 3) {stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress - biot_pp;}
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

void cauchystrain_3D_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_3D_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);
    const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]);  
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ; 

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress - biot_pp;
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

void cauchystress_planestrain_THM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_THM]; 
    PetscReal lambda = a[LAMBDA_THM];
    PetscReal poisson_r = a[POISSON_R_THM];
    PetscReal bulk_m = a[BULK_M_THM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM]; 
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy) + thermal_stress - biot_pp; 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_stress_planestrain_THM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_THM]; 
    PetscReal lambda = a[LAMBDA_THM];
    PetscReal poisson_r = a[POISSON_R_THM];
    PetscReal bulk_m = a[BULK_M_THM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM];  
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];
    
    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ;
    PetscScalar stress_zz;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp; 
    if (dim == 2) {stress_zz = poisson_r * (stress_xx + stress_yy + 2 * biot_pp) - biot_pp; }
    else if (dim == 3) {stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress - biot_pp;}
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

void cauchystress_3D_THM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    PetscReal shear_m = a[SHEAR_M_THM]; 
    PetscReal lambda = a[LAMBDA_THM];
    PetscReal bulk_m = a[BULK_M_THM];  
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM];  
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1] + disp_x[2*Nc+2];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;
    const PetscReal thermal_stress = - 3.0 * ther_expan * bulk_m * u[uOff[1]];
    const PetscScalar biot_pp = biot_alpha * u[uOff[2]] ; 

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace + thermal_stress - biot_pp;
    const PetscScalar stress_zz = twoshear_m*disp_x[2*Nc+2] + lambda_strainTrace + thermal_stress - biot_pp;
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

void f1_thermoporoelas_u_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]);

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

void g3_thermoporoelas_uu_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_THM]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_THM]); 
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

void g2_thermoporoelas_ut_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  const PetscInt  Nc     = dim;
  // const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
  // const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]);  
  PetscInt        c;

  for (c = 0; c < Nc; ++c) {
      g2[c*dim+c] -= 1;
  }
}

void displacement_THM_dimensionless(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar displacement[])
{
PetscInt i;

const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_THM]); 
const PetscReal bulk_m = PetscRealPart(constants[BULK_M_THM]);  
const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]); 
const PetscReal coeff = 3.0 * ther_expan * bulk_m / shear_m ;

const PetscScalar* disp = &u[uOff[0]];
for (i = 0; i < dim; i++) {displacement[i] = coeff * disp[i] ;}
}

void darcyvelocity_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar darcyvelocity[])
{
PetscInt i;
const PetscScalar* pp_x = &u_x[uOff_x[2]];
const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
PetscReal perm_over_visc = perm / fluid_visc;
const PetscReal grav = PetscRealPart(constants[GRAVITY_THM]);
const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_THM]);
// for (i = 0; i < dim; i++) {darcyvelocity[i] = darcyv[i] * 1.0e-9;}
for (i = 0; i < dim; i++) {darcyvelocity[i] = -perm_over_visc * pp_x[i] ;}
darcyvelocity[dim - 1] += perm_over_visc * grav * dens;

}

void darcyvelocity_THM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar darcyvelocity[])
{
PetscInt i;
const PetscScalar* pp_x = &u_x[uOff_x[2]];
PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
PetscReal perm = a[PERMEABILITY_THM];
PetscReal perm_over_visc = perm / fluid_visc;
PetscReal grav = a[GRAVITY_THM];
PetscReal dens = a[FLUID_DENSITY_THM];
for (i = 0; i < dim; i++) {darcyvelocity[i] = -perm_over_visc * pp_x[i] ;}
darcyvelocity[dim - 1] += perm_over_visc * grav * dens;

}

void f0_thermoporoelas_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    const PetscReal grav = PetscRealPart(constants[GRAVITY_THM]);
    const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_THM]);

    for (c = 0 ; c < dim; ++c) {
            f0[c] = u[uOff[2] + c];
    }

    f0[dim-1] -= grav * dens * perm_over_visc;


}

void f1_thermoporoelas_fluidflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0; c < dim; ++c) {
            f1[c*dim + c] = - perm_over_visc * u[uOff[3]];
    }
}

void g0_thermoporoelas_fluidflux_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;    

    for (c = 0 ; c < dim; ++c) {
            g0[c * dim + c] = 1.0;
    }

}

void g2_thermoporoelas_fluidflux_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[]) {
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0 ; c < dim; ++c) {
            g2[c * dim + c] = - perm_over_visc;
    }
}


void f0_thermoporoelas_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);
    const PetscReal biot_m = PetscRealPart(constants[BIOT_MODULUS_THM]);
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]); 

    f0[0] = u_t[uOff[2]] / biot_m ;
    f0[0] += biot_alpha * u_t[uOff[3]]  ; 
    f0[0] -= ther_expan * u_t[uOff[1]]  ; 
}

void f1_thermoporoelas_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    const PetscReal grav = PetscRealPart(constants[GRAVITY_THM]);
    const PetscReal dens = PetscRealPart(constants[FLUID_DENSITY_THM]);

    for (c = 0; c < dim; ++c) {
            f1[c] = perm_over_visc * u_x[uOff_x[2]+c];
    }
    f1[dim - 1] -= perm_over_visc * grav * dens;
}

void g1_thermoporoelas_pressure_flux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]){
    PetscInt c;
    for (c = 0; c < dim; ++c) {
        g1[c*dim + c] = 1.0;
    }
    
}

void g0_thermoporoelas_pressure_pressure_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;   
    const PetscReal biot_m = PetscRealPart(constants[BIOT_MODULUS_THM]);

    g0[0] = u_tShift / biot_m;

}

void g3_thermoporoelas_pressure_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
    PetscInt c;
    const PetscReal fluid_visc = PetscRealPart(constants[FLUID_VISCOSITY_THM]);
    const PetscReal perm = PetscRealPart(constants[PERMEABILITY_THM]);
    PetscReal perm_over_visc = perm / fluid_visc;
    
    PetscInt d;
    for (d = 0; d < dim; ++d) {
      g3[d * dim + d] = perm_over_visc;
    }    
}

void g0_thermoporoelas_pressure_strain_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;   
    const PetscReal biot_alpha = PetscRealPart(constants[BIOT_COEFFICIENT_THM]);
    
    g0[0] =  u_tShift * biot_alpha;

}

void g0_thermoporoelas_pressure_temperature_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    
    const PetscReal ther_expan = PetscRealPart(constants[LIN_THER_EXPAN_C_THM]); 
    
    g0[0] = -u_tShift * ther_expan;

}

PetscErrorCode initial_pressure_thermoporoelas(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
    PetscScalar *p = (PetscScalar *)ctx;

    *u = *p;

    return 0;
}

void f0_thermoporoelas_strain(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;

  for (d = 0; d < dim; ++d) f0[0] += u_x[uOff_x[0]+ d * dim + d];
  f0[0] -= u[uOff[3]];
}

void g1_thermoporoelas_strain_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) g1[d * dim + d] = 1.0;
}

void g0_thermoporoelas_strain_strain(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = -1.0;
}


void f0_thermoporoelas_fluidflux_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal grav = a[GRAVITY_THM];
    PetscReal dens = a[FLUID_DENSITY_THM];

    for (c = 0 ; c < dim; ++c) {
            f0[c] = u[uOff[2] + c];
    }

    f0[dim-1] -= grav * dens * perm_over_visc;


}

void f1_thermoporoelas_fluidflux_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0; c < dim; ++c) {
            f1[c*dim + c] = - perm_over_visc * u[uOff[3]];
    }
}


void g2_thermoporoelas_fluidflux_pressure_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[]) {
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;

    for (c = 0 ; c < dim; ++c) {
            g2[c * dim + c] = - perm_over_visc;
    }
}

void f0_thermoporoelas_pressure_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];
    PetscReal biot_m = a[BIOT_MODULUS_THM];
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM]; 
    PetscReal phi = a[POROSITY];
    PetscReal vol_fluid_ther_expan = a[VOL_FLUID_THER_EXPAN_C_THM]; 

    f0[0] = u_t[uOff[2]] / biot_m ; 

    f0[0] += biot_alpha * u_t[uOff[3]]  ; 

    f0[0] -= (biot_alpha * 3.0 * ther_expan + phi * (vol_fluid_ther_expan - 3.0 * ther_expan))* u_t[uOff[1]] ; 
}

void f1_thermoporoelas_pressure_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    PetscReal grav = a[GRAVITY_THM];
    PetscReal dens = a[FLUID_DENSITY_THM];

    for (c = 0; c < dim; ++c) {
            f1[c] = perm_over_visc * u_x[uOff_x[2]+c];
    }
    f1[dim - 1] -= perm_over_visc * grav * dens;
}

void g0_thermoporoelas_pressure_pressure_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;   
    PetscReal biot_m = a[BIOT_MODULUS_THM];

    g0[0] = u_tShift / biot_m;

}

void g3_thermoporoelas_pressure_pressure_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
    PetscInt c;
    PetscReal fluid_visc = a[FLUID_VISCOSITY_THM];
    PetscReal perm = a[PERMEABILITY_THM];
    PetscReal perm_over_visc = perm / fluid_visc;
    
    PetscInt d;
    for (d = 0; d < dim; ++d) {
      g3[d * dim + d] = perm_over_visc;
    }    
}

void g0_thermoporoelas_pressure_strain_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscInt c;   
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];
    
    g0[0] =  u_tShift * biot_alpha;

}

void g0_thermoporoelas_pressure_temperature_transient_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
    PetscReal biot_alpha = a[BIOT_COEFFICIENT_THM];
    PetscReal ther_expan = a[LIN_THER_EXPAN_C_THM]; 
    PetscReal phi = a[POROSITY];
    PetscReal vol_fluid_ther_expan = a[VOL_FLUID_THER_EXPAN_C_THM]; 

    g0[0] = - (biot_alpha * 3.0 * ther_expan + phi * (vol_fluid_ther_expan - 3.0 * ther_expan)) * u_tShift; 

}

// 아래 두 함수에 void *ctx 추가함. 
void property_user_THM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar property[])
{
  PetscScalar u1 = u[uOff[0]];
  PetscScalar u2= u[uOff[0]+1];
  PetscScalar u3;
  if (dim == 3) u3 = u[uOff[0]+2];
  PetscScalar T = u[uOff[1]];  // temperature in degC
  PetscScalar pp = u[uOff[3]]; // pore pressure in Pa
  PetscScalar x1= x[0];
  PetscScalar x2= x[1];
  PetscScalar x3;
  if (dim == 3) x3 = x[2];


  std::string *p = (std::string *)ctx;

  const std::string pro_eq = *p;
  try {
      mu::Parser parser;
      // if (temperature > 0) parser.SetExpr(pro_eq);
      // else parser.SetExpr("1.3306 * 0.000001");
      parser.SetExpr(pro_eq);
      parser.DefineVar("T", &T);
      parser.DefineVar("u1", &u1);
      parser.DefineVar("u2", &u2);
      if (dim ==3) parser.DefineVar("u3", &u3);
      parser.DefineVar("P", &pp);
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

