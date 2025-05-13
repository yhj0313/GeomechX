#include "PF_Problem.hh"
// #include "muParser.h"

namespace Pw_Functions_common
{
  PetscErrorCode material_property(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
   
    Property *p = (Property *)ctx;
    *u = p->GetValue();
  
    return 0;
  }

  PetscErrorCode material_property_nochange(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
  {
    *u = *u; 
    return 0;
  }

  void pwr_heatdecay(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                            const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                            PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar decay_to_model[])
  {
    PetscScalar x01(0.7805), x02(101.65), x03(863.15);
    PetscScalar y01(297.95), y02(32.19), y03(1.01);
    PetscScalar t11(2.94), t12(40.56), t13(32479.33);
    PetscScalar t21(1.1), t22(121.29), t23(9417.81);
    PetscScalar t31(42.75), t32(622.19), t33(622.7);
    PetscScalar A11(3218.38), A12(146.76), A13(12.05);
    PetscScalar A21(10394.94), A22(110.4), A23(20.78); 
    PetscScalar A31(2036.43), A32(197.22), A33(56.44); 

    PetscScalar decay_canister;
    PetscScalar t_year = t / 3.1536e7 ;
    PetscScalar *p = (PetscScalar*)ctx;
    const PetscScalar num_canisters_per_volume = *p;  // number per unit volume (m3)

    if (t_year <= 100.0) decay_canister = 1.72 * (y01 + A11 * exp(-(t_year + 40.0 - x01) / t11) + A21 * exp(-(t_year + 40.0- x01) / t21) + A31 * exp(-(t_year + 40.0- x01) / t31));
    else if (t_year <= 1000.0) decay_canister = 1.72 * (y02 + A12 * exp(-(t_year + 40.0 - x02) / t12) + A22 * exp(-(t_year + 40.0- x02) / t22) + A32 * exp(-(t_year + 40.0- x02) / t32));
    else if (t_year <= 100000.0) decay_canister = 1.72 * (y03 + A13 * exp(-(t_year + 40.0 - x03) / t13) + A23 * exp(-(t_year + 40.0- x03) / t23) + A33 * exp(-(t_year + 40.0- x03) / t33));
    else decay_canister = 0.0; 

    decay_to_model[0] = num_canisters_per_volume * decay_canister;

  }

  void nochange_prop_user(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                          const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                          PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], void *ctx, PetscScalar property[])
{


}

}
