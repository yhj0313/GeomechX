#include "PF_Elasticity.hh"

const char *bdvalue_call_E[6] = {"-bdvalue_elas_left", "-bdvalue_elas_right", "-bdvalue_elas_bottom", "-bdvalue_elas_top", "-bdvalue_elas_front", "-bdvalue_elas_back"};
const char *bdvalues_call_E[6] = {"-bdvalues_elas_left", "-bdvalues_elas_right", "-bdvalues_elas_bottom", "-bdvalues_elas_top", "-bdvalues_elas_front", "-bdvalues_elas_back"};

PetscScalar vBDVALUE = 0.0;
PetscScalar *BDVALUE = &vBDVALUE;

void Getbdvalue_M(PetscScalar &boundaryvalue, PetscInt bdlocation)
{
  PetscScalar bdv = 0.0;
  PetscOptionsGetScalar(NULL, NULL, bdvalue_call_E[bdlocation], &bdv, NULL);
  boundaryvalue = bdv;
}

void Getbdvalues(PetscScalar (&boundaryvalues) [3], PetscInt bdlocation)
{
  PetscScalar bdvs[3] = {0.0, 0.0, 0.0};
  PetscInt dim = 3;

  PetscOptionsGetScalarArray(NULL, NULL, bdvalues_call_E[bdlocation], bdvs, &dim, NULL);
  boundaryvalues[0] = bdvs[0];
  boundaryvalues[1] = bdvs[1];
  boundaryvalues[2] = bdvs[2];

}

namespace Pw_Functions{
  void f1_elas_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  const PetscInt  Nc = dim;
  const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_E]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_E]); 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      f1[c*dim+d] += shear_m*(u_x[c*dim+d] + u_x[d*dim+c]);
      f1[c*dim+c] += lambda*u_x[d*dim+d];
    }
  }

}

void g3_elas_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  const PetscInt  Nc     = dim;
  const PetscReal shear_m  = PetscRealPart(constants[SHEAR_M_E]); 
  const PetscReal lambda = PetscRealPart(constants[LAMBDA_E]); 
  PetscInt        c, d;

  for (c = 0; c < Nc; ++c) {
    for (d = 0; d < dim; ++d) {
      g3[((c*Nc + c)*dim + d)*dim + d] += shear_m;
      g3[((c*Nc + d)*dim + d)*dim + c] += shear_m;
      g3[((c*Nc + d)*dim + c)*dim + d] += lambda;
    }
  }
}

PetscErrorCode zero_disp(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  PetscInt d;
  for (d = 0; d < dim; ++d) u[d] = 0.0;
  return 0;
}

PetscErrorCode disp_bd(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  PetscInt d;
  for (d = 0; d < dim; ++d) u[d] = 0.1;
  return 0;
}

PetscErrorCode boundary_displacement(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{

  BoundaryCondition *p = (BoundaryCondition *)ctx;

  for (int d = 0; d < dim; ++d) u[d] = p->GetBdValue();
  return 0;
}


PetscErrorCode boundary_normalstress(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  PetscScalar *p = (PetscScalar *)ctx;

  *u = *p;
  return 0;
}

void f0_normal_stress_bd_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;

  // {PetscPrintf(PETSC_COMM_WORLD, "aux[0] is %f and x,y is %f, %f \n", a[0], x[0], x[1]);}

  for (d = 0; d < dim; ++d) { 
    f0[d] = - a[0] * n[d]; 
  }
}

void f0_normal_stress_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)
  N = *BDVALUE;
  // {PetscPrintf(PETSC_COMM_WORLD, "N is %f \n", N);}

  for (d = 0; d < dim; ++d) { 
    f0[d] = - N * n[d]; 
  }
}

void f0_normal_stress_bd_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 0);
  for (d = 0; d < dim; ++d) { 
    f0[d] = - N * n[d]; 
    // f0[d] = - N * x[d]; 
  }
  PetscPrintf(PETSC_COMM_WORLD, "x is %f, %f, %f \n normal vector is %f, %f, %f \n", x[0], x[1], x[2], n[0], n[1], n[2]);
}

void f0_normal_stress_bd_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 1);
  for (d = 0; d < dim; ++d) { 
    f0[d] = - N * n[d]; 
  }
  PetscPrintf(PETSC_COMM_WORLD, "x is %f, %f, %f \n normal vector is %f, %f, %f \n", x[0], x[1], x[2], n[0], n[1], n[2]);

}
void f0_normal_stress_bd_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 2);
  for (d = 0; d < dim; ++d) { 
    f0[d] = -N * n[d]; 
  }
}
void f0_normal_stress_bd_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 3);
  for (d = 0; d < dim; ++d) { 
    f0[d] = -N * n[d]; 
  }
}
void f0_normal_stress_bd_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 4);
  for (d = 0; d < dim; ++d) { 
    f0[d] = -N * n[d]; 
  }
}
void f0_normal_stress_bd_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt        d;
  PetscReal N ; // normal stress (sigma_n) in y direction (front and back boundaries in 3D, top and bottom in 2D)

  Getbdvalue_M(N, 5);
  for (d = 0; d < dim; ++d) { 
    f0[d] = -N * n[d]; 
  }
}

void f0_traction_bd_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 0);
  f0[0] = - Traction[0];
  f0[1] = - Traction[1]; 
  if (dim == 3) {f0[2] = - Traction[2]; }
}

void f0_traction_bd_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 1);
  f0[0] = -Traction[0];
  f0[1] = -Traction[1]; 
  if (dim == 3) {f0[2] = -Traction[2]; }
}

void f0_traction_bd_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 2);
  f0[0] = -Traction[0];
  f0[1] = -Traction[1]; 
  if (dim == 3) {f0[2] = -Traction[2]; }
}

void f0_traction_bd_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 3);
  f0[0] = -Traction[0];
  f0[1] = -Traction[1]; 
  if (dim == 3) {f0[2] = -Traction[2]; }
}

void f0_traction_bd_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 4);
  f0[0] = -Traction[0];
  f0[1] = -Traction[1]; 
  if (dim == 3) {f0[2] = -Traction[2]; }
}

void f0_traction_bd_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal Traction[3];

  Getbdvalues(Traction, 5);
  f0[0] = -Traction[0];
  f0[1] = -Traction[1]; 
  if (dim == 3) {f0[2] = -Traction[2]; }
}


void cauchystrain_3D(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_3D(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_E]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_E]);

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

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_zz; // stress_zz
    stress[3] = stress_xy; // stress_xy
    stress[4] = stress_yz; // stress_yz
    stress[5] = stress_xz; // stress_xz
} // stress

void cauchystrain_planestrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void cauchystress_planestrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_E]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_E]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_E]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace;
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace;
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy); 
    const PetscScalar stress_xy = shear_m * (disp_x[0*Nc+1] + disp_x[1*Nc+0]);

    stress[0] = stress_xx; // stress_xx
    stress[1] = stress_yy; // stress_yy
    stress[2] = stress_xy; // stress_xy
    stress[3] = stress_zz; // stress_zz
} // stress

void axisymmetric_2d_stress_planestrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
    const PetscInt Nc = dim;
    const PetscReal shear_m = PetscRealPart(constants[SHEAR_M_E]); 
    const PetscReal lambda = PetscRealPart(constants[LAMBDA_E]);
    const PetscReal poisson_r = PetscRealPart(constants[POISSON_R_E]);

    const PetscScalar* disp_x = &u_x[uOff_x[0]];
    const PetscReal strainTrace = disp_x[0*Nc+0] + disp_x[1*Nc+1];
    const PetscReal lambda_strainTrace = lambda * strainTrace;
    const PetscReal twoshear_m = 2.0*shear_m;

    const PetscScalar stress_xx = twoshear_m*disp_x[0*Nc+0] + lambda_strainTrace; 
    const PetscScalar stress_yy = twoshear_m*disp_x[1*Nc+1] + lambda_strainTrace; 
    const PetscScalar stress_zz = poisson_r * (stress_xx + stress_yy); 
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

void axisymmetric_2d_strain_planestrain(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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

void f1_TI(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
   /* Transversely isotropic material with the (x, z) plane of isotropy, 
   E_x = E_z = E, E_y = E', mu_xz = mu, mu_yx = mu_zy = mu', mu_xy = mu_yz = mu' * E/E'
  G_xy = G_yz = G', G_xz = G */  

  const PetscReal D11 = PetscRealPart(constants[0]) ;    
  const PetscReal D12 = PetscRealPart(constants[1]);
  const PetscReal D13 = PetscRealPart(constants[2]); 
  const PetscReal D16 = PetscRealPart(constants[3]); 
  // const PetscReal D16 = 0.0 ;
  const PetscReal D22 = PetscRealPart(constants[4]); 
  const PetscReal D23 = PetscRealPart(constants[5]); 
  const PetscReal D26 = PetscRealPart(constants[6]); 
  // const PetscReal D26 = 0.0; 
  const PetscReal D33 = PetscRealPart(constants[7]); 
  const PetscReal D36 = PetscRealPart(constants[8]);
  const PetscReal D44 = PetscRealPart(constants[9]); 
  const PetscReal D45 = PetscRealPart(constants[10]); 
  const PetscReal D55 = PetscRealPart(constants[11]); 
  const PetscReal D66 = PetscRealPart(constants[12]);  


  f1[0] = D11 * u_x[0] + D12 * u_x[4] + D13 * u_x[8] + D16 * (u_x[1] + u_x[3]);
  f1[4] = D12 * u_x[0] + D22 * u_x[4] + D23 * u_x[8] + D26 * (u_x[1] + u_x[3]); 
  f1[8] = D13 * u_x[0] + D23 * u_x[4] + D33 * u_x[8] + D36 * (u_x[1] + u_x[3]);

  f1[1] = D16 * u_x[0] + D26 * u_x[4] + D36 * u_x[8] + D66 * (u_x[1] + u_x[3]); 

  f1[2] = D45 * (u_x[5] + u_x[7]) + D55 * (u_x[2] + u_x[6]);
  f1[3] = D16 * u_x[0] + D26 * u_x[4] + D36 * u_x[8] + D66 * (u_x[1] + u_x[3]); 

  f1[5] = D44 * (u_x[5] + u_x[7]) + D45 * (u_x[2] + u_x[6]);
  f1[6] = D45 * (u_x[5] + u_x[7]) + D55 * (u_x[2] + u_x[6]);
  f1[7] = D44 * (u_x[5] + u_x[7]) + D45 * (u_x[2] + u_x[6]);

}                

void g3_TI(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  /* Transversely isotropic material with the (x, z) plane of isotropy, 
   E_x = E_z = E, E_y = E', mu_xz = mu, mu_yx = mu_zy = mu', mu_xy = mu_yz = mu' * E/E'
  G_xy = G_yz = G', G_xz = G */  

  const PetscReal D11 = PetscRealPart(constants[0]) ;    
  const PetscReal D12 = PetscRealPart(constants[1]);
  const PetscReal D13 = PetscRealPart(constants[2]); 
  const PetscReal D16 = PetscRealPart(constants[3]); 
  // const PetscReal D16 = 0.0;
  const PetscReal D22 = PetscRealPart(constants[4]); 
  const PetscReal D23 = PetscRealPart(constants[5]); 
  const PetscReal D26 = PetscRealPart(constants[6]);
  // const PetscReal D26 = 0.0;  
  const PetscReal D33 = PetscRealPart(constants[7]); 
  const PetscReal D36 = PetscRealPart(constants[8]);
  const PetscReal D44 = PetscRealPart(constants[9]); 
  const PetscReal D45 = PetscRealPart(constants[10]); 
  const PetscReal D55 = PetscRealPart(constants[11]); 
  const PetscReal D66 = PetscRealPart(constants[12]);  

  g3[0] = D11;
  g3[4] = D12; 
  g3[8] = D13;
  g3[1] = D16;
  g3[3] = D16;

  g3[36] = D12;
  g3[40] = D22; 
  g3[44] = D23;
  g3[37] = D26;
  g3[39] = D26;

  g3[72] = D13; 
  g3[76] = D23; 
  g3[80] = D33;
  g3[73] = D36;
  g3[75] = D36;

  g3[9] = D16;
  g3[13] = D26;
  g3[17] = D36;
  g3[10] = D66;
  g3[12] = D66;

  g3[23] = D45;
  g3[25] = D45;
  g3[20] = D55;
  g3[24] = D55;
  // g3[26] = D55;

  g3[27] = D16;
  g3[31] = D26;
  g3[35] = D36;
  g3[28] = D66;
  g3[30] = D66;

  g3[50] = D44;
  g3[52] = D44;
  g3[47] = D45;
  g3[51] = D45;

  g3[59] = D45;
  g3[61] = D45;
  g3[56] = D55;
  g3[60] = D55;

  g3[68] = D44; 
  g3[70] = D44; 
  g3[65] = D45; 
  g3[69] = D45; 

}
   
void axisymmetric_3d_stress_TI(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
  const PetscReal D11 = PetscRealPart(constants[0]) ;    
  const PetscReal D12 = PetscRealPart(constants[1]);
  const PetscReal D13 = PetscRealPart(constants[2]); 
  const PetscReal D16 = PetscRealPart(constants[3]); 
  const PetscReal D22 = PetscRealPart(constants[4]); 
  const PetscReal D23 = PetscRealPart(constants[5]); 
  const PetscReal D26 = PetscRealPart(constants[6]); 
  const PetscReal D33 = PetscRealPart(constants[7]); 
  const PetscReal D36 = PetscRealPart(constants[8]);
  const PetscReal D44 = PetscRealPart(constants[9]); 
  const PetscReal D45 = PetscRealPart(constants[10]); 
  const PetscReal D55 = PetscRealPart(constants[11]); 
  const PetscReal D66 = PetscRealPart(constants[12]);  

  const PetscScalar stress_xx = D11 * u_x[0] + D12 * u_x[4] + D13 * u_x[8] + D16 * (u_x[1] + u_x[3]); 
  const PetscScalar stress_yy = D12 * u_x[0] + D22 * u_x[4] + D23 * u_x[8] + D26 * (u_x[1] + u_x[3]); 
  const PetscScalar stress_zz = D13 * u_x[0] + D23 * u_x[4] + D33 * u_x[8] + D36 * (u_x[1] + u_x[3]);
  const PetscScalar stress_xz = D45 * (u_x[5] + u_x[7]) + D55 * (u_x[2] + u_x[6]);

//(PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sy + 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar stress_rr = PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) * stress_zz + PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) * stress_xx + 2 * PetscSinReal(PetscAtan2Real(x[0],x[2])) * PetscCosReal(PetscAtan2Real(x[0],x[2])) * stress_xz;
  // (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sy - 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar stress_thetatheta = PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) * stress_zz + PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) * stress_xx - 2*PetscSinReal(PetscAtan2Real(x[0],x[2])) * PetscCosReal(PetscAtan2Real(x[0],x[2])) * stress_xz;
  // - (PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x))) * solid.sx +(PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x)))  * solid.sy +((PetscCosReal(PetscAtan2Real(y,x)))^2.0 - (PetscSinReal(PetscAtan2Real(y,x)))^2.0 ) * solid.sxy
  const PetscScalar stress_rtheta = - (PetscCosReal(PetscAtan2Real(x[0],x[2]))) * (PetscSinReal(PetscAtan2Real(x[0],x[2]))) * stress_zz +(PetscCosReal(PetscAtan2Real(x[0],x[2]))) * (PetscSinReal(PetscAtan2Real(x[0],x[2])))  * stress_xx +(PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) - PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) ) * stress_xz; 


  stress[0] = stress_rr; // stress_rr (radial stress)
  stress[1] = stress_thetatheta; // stress_thetatheta (tangential stress)
  stress[2] = stress_rtheta; // stress_rtheta
  stress[3] = stress_yy; // stress_yy
} // stress

void axisymmetric_3d_strain_TI(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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
  const PetscScalar strain_xz = 0.5*(disp_x[2*Nc+0] + disp_x[0*Nc+2]);

  //(PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sy + 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar strain_rr = PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) * strain_zz + PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) * strain_xx + 2*PetscSinReal(PetscAtan2Real(x[0],x[2])) * PetscCosReal(PetscAtan2Real(x[0],x[2])) * strain_xz;
  // (PetscSinReal(PetscAtan2Real(y,x)))^2.0 * solid.sx + (PetscCosReal(PetscAtan2Real(y,x)))^2.0 * solid.sy - 2*PetscSinReal(PetscAtan2Real(y,x)) * PetscCosReal(PetscAtan2Real(y,x)) * solid.sxy
  const PetscScalar strain_thetatheta = PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) * strain_zz + PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) * strain_xx - 2*PetscSinReal(PetscAtan2Real(x[0],x[2])) * PetscCosReal(PetscAtan2Real(x[0],x[2])) * strain_xz;
  // - (PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x))) * solid.sx +(PetscCosReal(PetscAtan2Real(y,x))) * (PetscSinReal(PetscAtan2Real(y,x)))  * solid.sy +((PetscCosReal(PetscAtan2Real(y,x)))^2.0 - (PetscSinReal(PetscAtan2Real(y,x)))^2.0 ) * solid.sxy
  const PetscScalar strain_rtheta = - (PetscCosReal(PetscAtan2Real(x[0],x[2]))) * (PetscSinReal(PetscAtan2Real(x[0],x[2]))) * strain_zz +(PetscCosReal(PetscAtan2Real(x[0],x[2]))) * (PetscSinReal(PetscAtan2Real(x[0],x[2])))  * strain_xx +(PetscPowReal(PetscCosReal(PetscAtan2Real(x[0],x[2])), 2) - PetscPowReal(PetscSinReal(PetscAtan2Real(x[0],x[2])), 2) ) * strain_xz; 


  strain[0] = strain_rr; // strain_rr (radial strain)
  strain[1] = strain_thetatheta; // strain_thetatheta (tangential strain)
  strain[2] = strain_rtheta; // strain_rtheta 
  strain[3] = strain_yy; // strain_yy

}

void cauchystress_3d_TI(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar stress[])
{
  const PetscReal D11 = PetscRealPart(constants[0]) ;    
  const PetscReal D12 = PetscRealPart(constants[1]);
  const PetscReal D13 = PetscRealPart(constants[2]); 
  const PetscReal D16 = PetscRealPart(constants[3]); 
  const PetscReal D22 = PetscRealPart(constants[4]); 
  const PetscReal D23 = PetscRealPart(constants[5]); 
  const PetscReal D26 = PetscRealPart(constants[6]); 
  const PetscReal D33 = PetscRealPart(constants[7]); 
  const PetscReal D36 = PetscRealPart(constants[8]);
  const PetscReal D44 = PetscRealPart(constants[9]); 
  const PetscReal D45 = PetscRealPart(constants[10]); 
  const PetscReal D55 = PetscRealPart(constants[11]); 
  const PetscReal D66 = PetscRealPart(constants[12]);  

  const PetscScalar stress_xx = D11 * u_x[0] + D12 * u_x[4] + D13 * u_x[8] + D16 * (u_x[1] + u_x[3]); 
  const PetscScalar stress_yy = D12 * u_x[0] + D22 * u_x[4] + D23 * u_x[8] + D26 * (u_x[1] + u_x[3]); 
  const PetscScalar stress_zz = D13 * u_x[0] + D23 * u_x[4] + D33 * u_x[8] + D36 * (u_x[1] + u_x[3]);
  // const PetscScalar stress_zz = D13 * u_x[0] + D23 * u_x[4] + D33 * u_x[8] ;
  const PetscScalar stress_yz = D44 * (u_x[5] + u_x[7]) + D45 * (u_x[2] + u_x[6]);
  const PetscScalar stress_xz = D45 * (u_x[5] + u_x[7]) + D55 * (u_x[2] + u_x[6]);
  // const PetscScalar stress_xy = D16 * u_x[0] + D26 * u_x[4] + D66 * (u_x[1] + u_x[3]); 
    const PetscScalar stress_xy = D16 * u_x[0] + D26 * u_x[4] + D36 * u_x[8] + D66 * (u_x[1] + u_x[3]);   

  stress[0] = stress_xx; // stress_xx
  stress[1] = stress_yy; // stress_yy
  stress[2] = stress_zz; // stress_zz
  stress[3] = stress_xy; // stress_xy
  stress[4] = stress_yz; // stress_yz
  stress[5] = stress_xz; // stress_xz
} // stress


}
