#include "PF_Thermal.hh"

const char *bdvalue_call_T[6] = {"-bdvalue_heat_left", "-bdvalue_heat_right", "-bdvalue_heat_bottom", "-bdvalue_heat_top", "-bdvalue_heat_front", "-bdvalue_heat_back"};

void Getbdvalue(PetscScalar &boundaryvalue, PetscInt bdlocation)
{
  PetscScalar bdv = 0.0;
  PetscOptionsGetScalar(NULL, NULL, bdvalue_call_T[bdlocation], &bdv, NULL);
  boundaryvalue = bdv;
}

namespace Pw_Functions
{
  PetscErrorCode boundary_temperature(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    BoundaryCondition *p = (BoundaryCondition *)ctx;
    // Thermal::BoundaryCondition *p = (Thermal::BoundaryCondition *)ctx;
    // std::cout << p->GetBdValue() << "yeah" << std::endl;
    *u = p->GetBdValue();
    return 0;
  }

  void f0_insulation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscInt d;
    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);

    for (d = 0; d < dim; ++d)
      f0[0] = u_x[d] * n[d]; 
    f0[0] *= ther_cond;
  }

  void f0_insulation_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscInt d;
    
    for (d = 0; d < dim; ++d){
      f0[0] += u_x[d] * n[d]; 
    }
    f0[0] *= a[THER_COND_T];
  }

  void f0_heatflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscInt d;

    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);
    for (d = 0; d < dim; ++d) {
      f0[0] += u_x[d] * n[d]; 
    }      
    f0[0] = f0[0] * ther_cond - *BDVALUE;
  }

  void f0_heatflux_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscInt d;

    for (d = 0; d < dim; ++d)
      f0[0] = u_x[d] * n[d]; 
    f0[0] = f0[0] * a[THER_COND_T] - *BDVALUE;
  }



  void f0_heatflux_0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 0); // value on the left boundary

    f0[0] = -q0;
  }

  void f0_heatflux_1(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 1); // value on the right boundary

    f0[0] = -q0;
  }

  void f0_heatflux_2(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 2); // value on the bottom boundary    

    f0[0] = -q0;
  }

  void f0_heatflux_3(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 3); // value on the top boundary        

    f0[0] = -q0;
  }

  void f0_heatflux_4(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 4); // value on the front boundary    

    f0[0] = -q0;
  }

  void f0_heatflux_5(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 5); // value on the back boundary        

    f0[0] = -q0;
  }

  void f0_heatflux_0_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 0); // value on the left boundary

    f0[0] = -q0;
  }

  void f0_heatflux_1_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 1); // value on the right boundary

    f0[0] = -q0;
  }

  void f0_heatflux_2_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 2); // value on the bottom boundary    

    f0[0] = -q0;
  }

  void f0_heatflux_3_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 3); // value on the top boundary        

    f0[0] = -q0;
  }

  void f0_heatflux_4_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 4); // value on the front boundary    

    f0[0] = -q0;
  }

  void f0_heatflux_5_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                   const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                   PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal q0;
    
    Getbdvalue(q0, 5); // value on the back boundary        

    f0[0] = -q0;
  }

  void f0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_T]);
    const PetscReal density = PetscRealPart(constants[DENSITY_T]);
    const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_T]);
    f0[0] = density * heat_capa * flow_vel * u_x[0]; // convection (일단 x방향만 고려)
  
  }

  void f1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
  {
    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);
    PetscInt d;
    for (d = 0; d < dim; ++d)
      f1[d] = ther_cond * u_x[d];
  }

  void g1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
  { // convection (일단 x방향만 고려)

    const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_T]);
    const PetscReal density = PetscRealPart(constants[DENSITY_T]);
    const PetscReal flow_vel = PetscRealPart(constants[FLOW_VEL_T]);
    g1[0] = density * heat_capa * flow_vel;
    g1[1] = 0;
    g1[1] = 0;
  }

  void g3_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
  {
    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);
    PetscInt d;
    for (d = 0; d < dim; ++d)
      g3[d * dim + d] = ther_cond;
  }

  
  void f0_temp_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {

    f0[0] = a[DENSITY_T] * a[HEAT_CAPA_T] * a[FLOW_VEL_T] * u_x[0] ;  // convection (일단 x방향만 고려)
  
  }

  void f1_temp_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
  {
    PetscInt d;
    
    for (d = 0; d < dim; ++d)
      f1[d] = a[THER_COND_T] * u_x[d];
    
  }

  void g1_temp_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
  { // convection (일단 x방향만 고려)

    g1[0] = a[DENSITY_T] * a[HEAT_CAPA_T] * a[FLOW_VEL_T];
    g1[1] = 0;
    g1[1] = 0;
  }

  void g3_temp_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
  {

    PetscInt d;
    for (d = 0; d < dim; ++d)
      g3[d * dim + d] = a[THER_COND_T]; 

  }

  void g1_insulation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
  {
    PetscInt d;
    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);
  
    for (d = 0; d < dim; ++d) {
      g1[d] = ther_cond * n[d]; 

    }
      

  }

  void g1_insulation_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
  {
    PetscInt d;
    for (d = 0; d < dim; ++d) {
      g1[d] = a[THER_COND_T] * n[d]; 
    }
  }

  void heatflux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar heatflux[])
  {
    const PetscReal ther_cond = PetscRealPart(constants[THER_COND_T]);
    PetscInt d;

    for (d = 0; d < dim; ++d)
      heatflux[d] = -ther_cond * u_x[d];
  }

  PetscErrorCode bc_auxiliary(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    BoundaryCondition *p = (BoundaryCondition *)ctx;
    *u = p->GetBdValue();
    return 0;
  }

  PetscErrorCode initial_temperature(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    PetscScalar *p = (PetscScalar *)ctx;

    *u = *p;

    return 0;
  }

  PetscErrorCode initial_temperature_user(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    PetscScalar x1= x[0];
    PetscScalar x2= x[1];
    PetscScalar x3;
    if (dim == 3) x3 = x[2];

    std::string *p = (std::string *)ctx;
    const std::string pro_eq = *p;

    try {
      mu::Parser parser;
      parser.SetExpr(pro_eq);
      parser.DefineVar("x", &x1);
      parser.DefineVar("y", &x2);
      if (dim ==3) parser.DefineVar("z", &x3);
      // parser.Compile();
      u[0] = parser.Eval();
      std::cout << u[0] << std::endl;
  } catch (mu::Parser::exception_type& e) {
      std::cout << e.GetMsg() << std::endl;
  }

    return 0;
  }

  PetscErrorCode initial_temperature_t(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    *u = 0.0;
    return 0;
  }

  void f0_test(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_T]);
    const PetscReal density = PetscRealPart(constants[DENSITY_T]);
    f0[0] = density * heat_capa * u_t[0];
  }

  void g0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
  {
    const PetscReal heat_capa = PetscRealPart(constants[HEAT_CAPA_T]);
    const PetscReal density = PetscRealPart(constants[DENSITY_T]);
    g0[0] = u_tShift * heat_capa * density;
  }

  void f0_test_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
  {
    PetscReal heat_source = a[HEAT_SOURCE_T];

    f0[0] = a[DENSITY_T] * a[HEAT_CAPA_T] * u_t[0];
    f0[0] -= heat_source; 
  }

  void g0_temp_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
               const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
               PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
  {
    g0[0] = u_tShift * a[HEAT_CAPA_T] * a[DENSITY_T];
  }

  /* TM */

void f0_insulation_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);

  for (d = 0; d < dim; ++d)
    f0[0] = u_x[d] * n[d]; 
  f0[0] *= ther_cond;
}

void g1_insulation_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  const PetscReal ther_cond = PetscRealPart(constants[THER_COND_TM]);

  for (d = 0; d < dim; ++d)
    g1[d] = ther_cond * n[d]; 

}


void f0_heatflux_0_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 0); // value on the left boundary

  f0[0] = -q0;
}

void f0_heatflux_1_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;

  Getbdvalue(q0, 1); // value on the right boundary

  f0[0] = -q0;
}

void f0_heatflux_2_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 2); // value on the bottom boundary    

  f0[0] = -q0;
}

void f0_heatflux_3_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 3); // value on the top boundary        

  f0[0] = -q0;
}

void f0_heatflux_4_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 4); // value on the front boundary    

  f0[0] = -q0;
}

void f0_heatflux_5_TM(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 5); // value on the back boundary        

  f0[0] = -q0;
}


void f0_insulation_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;

  const PetscReal ther_cond = a[THER_COND_TM];

  for (d = 0; d < dim; ++d)
    f0[0] = u_x[d] * n[d]; 
  f0[0] *= ther_cond;
}

void g1_insulation_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;

  const PetscReal ther_cond = a[THER_COND_TM];

  for (d = 0; d < dim; ++d)
    g1[d] = ther_cond * n[d]; 

}


void f0_heatflux_0_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 0); // value on the left boundary

  f0[0] = -q0;
}

void f0_heatflux_1_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;

  Getbdvalue(q0, 1); // value on the right boundary

  f0[0] = -q0;
}

void f0_heatflux_2_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 2); // value on the bottom boundary    

  f0[0] = -q0;
}

void f0_heatflux_3_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 3); // value on the top boundary        

  f0[0] = -q0;
}

void f0_heatflux_4_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 4); // value on the front boundary    

  f0[0] = -q0;
}

void f0_heatflux_5_TM_mat_fields(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscReal q0;
  
  Getbdvalue(q0, 5); // value on the back boundary        

  f0[0] = -q0;
}

void property_user_T(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar property[])
{

  PetscScalar T = u[uOff[0]];  // temperature in degC
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
