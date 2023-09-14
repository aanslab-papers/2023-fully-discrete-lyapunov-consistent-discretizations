#include <ssdc.h>
#include <assert.h>
#include "model-math.h"

#define Sign PetscSignReal
#define Sqrt PetscSqrtReal
#define Pi   PETSC_PI
#define Sin  PetscSinReal
#define Cos  PetscCosReal
#define Exp  PetscExpReal
#define Sqr  PetscSqr
#define Pow  PetscPowReal
#define Log  PetscLogReal

#define NAME "M1VIRUS"
enum {dof = 5};

#define HAVE_FIELDS 1
static const SSDC_Field fields[] = {
  { "S", 1,  0 },
  { "N", 1,  1 },
  { "T", 1,  2 },
  { "V", 1,  3 },
  { "Z", 1,  4 },
  { "",  0,  0 }
};

#define DIM 1
enum {dim = DIM};

#ifndef HAVE_DIFFUSION
#define HAVE_DIFFUSION 1
#endif

#if HAVE_DIFFUSION == 1
#define OPT_DEG "1"
#define OPT_NEL "16"
#if DIM == 3
#define OPT_BOX "-67,+67,-108,+76,-25,+100"
#else
#define OPT_BOX "0,2"
#endif
#endif

static PetscInt    CASE  = 1;

static ssdc_real_t A     = 0.02;
static ssdc_real_t B     = 0.01;
static ssdc_real_t d     = 0.02;
static ssdc_real_t beta1 = 0.15;
static ssdc_real_t beta2 = 0.5;
static ssdc_real_t beta3 = 0.05;
static ssdc_real_t beta4 = 0.6;
static ssdc_real_t r1    = 0.8;
static ssdc_real_t r2    = 0.9;
static ssdc_real_t r3    = 0.5;
static ssdc_real_t r4    = 0.8;
static ssdc_real_t eps1  = 0.008;
static ssdc_real_t eps2  = 0.008;
static ssdc_real_t eps3  = 0.005;
static ssdc_real_t eps4  = 0.02;
static ssdc_real_t DS    = 0.01;
static ssdc_real_t DN    = 0.01;
static ssdc_real_t DT    = 0.01;
static ssdc_real_t DV    = 0.03;
static ssdc_real_t DZ    = 0.03;
static ssdc_real_t Cdiff[dof][dof];

static ssdc_real_t A1 = 0;
static ssdc_real_t A2 = 0;
static ssdc_real_t A3 = 0;

static ssdc_real_t S_eq = 0;
static ssdc_real_t N_eq = 0;
static ssdc_real_t T_eq = 0;
static ssdc_real_t V_eq = 0;
static ssdc_real_t Z_eq = 0;

#define ParamInit(name) \
  do {name = SSDCOptionReal("-"#name, name);} while(0)

// Should include X to compute space dependant parameters Lambda and Mu
PetscErrorCode ModelInitialize(void)
{
  CASE = SSDCOptionInt("-case", CASE);
  assert(CASE>=0);
  assert(CASE<=5);

  ParamInit( A     );
  ParamInit( B     );
  ParamInit( d     );
  ParamInit( beta1 );
  ParamInit( beta2 );
  ParamInit( beta3 );
  ParamInit( beta4 );
  ParamInit( r1    );
  ParamInit( r2    );
  ParamInit( r3    );
  ParamInit( r4    );
  ParamInit( eps1  );
  ParamInit( eps2  );
  ParamInit( eps3  );
  ParamInit( eps4  );

  ParamInit( DS    );
  ParamInit( DN    );
  ParamInit( DT    );
  ParamInit( DV    );
  ParamInit( DZ    );
  Cdiff[0][0] = DS;
  Cdiff[1][1] = DN;
  Cdiff[2][2] = DT;
  Cdiff[3][3] = DV;
  Cdiff[4][4] = DZ;

  switch (CASE) {
  case 0:
    beta1 = 0.03;
    beta2 = 0.03;
    beta3 = 0.1;
    beta4 = 0.03;
    eps1  = 0.04;
    eps2  = 0.01;
    eps3  = 0.008;
    eps4  = 0.01;
    r2    = 0.8;
    break;
  case 1:
    beta1 = 0.1;
    beta2 = 0.03;
    beta3 = 0.1;
    beta4 = 0.03;
    eps1  = 0.008;
    eps2  = 0.01;
    eps3  = 0.006;
    eps4  = 0.01;
    r2    = 0.8;
    break;
  case 2:
    beta1 = 0.03;
    beta2 = 0.1;
    beta3 = 0.1;
    beta4 = 0.2;
    eps1  = 0.04;
    eps2  = 0.008;
    eps3  = 0.008;
    eps4  = 0.1;
    r2    = 0.8;
    break;
  case 3:
    beta1 = 0.15;
    beta2 = 0.35;
    beta3 = 0.1;
    beta4 = 0.1;
    eps1  = 0.008;
    eps2  = 0.008;
    eps3  = 0.008;
    eps4  = 0.01;
    r2    = 0.9;
    break;
  case 4:
    beta1 = 0.04;
    beta2 = 0.5;
    beta3 = 0.1;
    beta4 = 0.6;
    eps1  = 0.05;
    eps2  = 0.008;
    eps3  = 0.005;
    eps4  = 0.02;
    r2    = 0.9;
    break;
  case 5:
    beta1 = 0.15;
    beta2 = 0.5;
    beta3 = 0.05;
    beta4 = 0.6;
    eps1  = 0.008;
    eps2  = 0.008;
    eps3  = 0.005;
    eps4  = 0.02;
    r2    = 0.9;
    break;
  }

  A1 = A*r1*beta1/(d*(d+eps1));
  A2 = A*r2*beta2/(d*(d+eps2));
  A3 = r4*beta4*(d+eps3)/(r3*beta3*(d+eps4));

  switch (CASE) {
  case 0:
    S_eq = A/d;
    N_eq = 0;
    T_eq = 0;
    V_eq = B/(d+eps3);
    Z_eq = 0;
    break;
  case 1:
    S_eq = (d+eps1)/(r1*beta1);
    N_eq = d/beta1*(A1-1);
    T_eq = 0;
    V_eq = B/(d+eps3);
    Z_eq = 0;
    break;
  case 2: {
    ssdc_real_t a1 = beta3*(r3*beta3*d+beta2*(d+eps3));
    ssdc_real_t a2 = a1/beta3*(d+eps2)-beta2*beta3*(B+r2*r3*A);
    ssdc_real_t a3 = -B*beta2*(d+eps2);
    ssdc_real_t delta = a2*a2-4*a1*a3;
    ssdc_real_t V2 = (-a2+sqrt(delta))/(2*a1);
    S_eq = (beta3*V2+d+eps2)/(r2*beta2);
    N_eq = 0;
    T_eq = -d/beta2 + (A*r2)/(beta3*V2+d+eps2);
    V_eq = V2;
    Z_eq = 0;
    } break;
  case 3: {
    ssdc_real_t V3 = (d+eps2)/beta3*(A2/A1-1);
    S_eq = (d+eps1)/(r1*beta1);
    N_eq = (A*r1*r3*beta1*beta3-r3*beta3*d*(d+eps1)-beta2*(d+eps1)*(d+eps3))/(r3*beta1*beta3*(d+eps1))
      /**/ + (B*beta2)/(r3*beta1*(d+eps2)*(A2/A1-1));
    T_eq = (-B+(d+eps3)*V3)/(r3*beta3*V3);
    V_eq = V3;
    Z_eq = 0;
    } break;
  case 4: {
    ssdc_real_t q = (r3*beta3*beta4*(beta2*(d+eps4)+r4*beta4*d)*(d+eps4)*(A3-1));
    S_eq = (A*r4*beta4)/(beta2*(d+eps4)+r4*beta4*d);
    N_eq = 0;
    T_eq = (d+eps4)/(r4*beta4);
    V_eq = (B*r4*beta4)/(r3*beta3*(d+eps4)*(A3-1));
    Z_eq = (A*r2*beta2*r3*beta3*r4*beta4*(d+eps4)*(A3-1))/q
      /**/ - ((beta2*(d+eps4)+r4*beta4*d)*(r3*beta3*(d+eps2)*(d+eps4)*(A3-1)+B*beta3*r4*beta4))/q;
    } break;
  case 5:
    S_eq = (d+eps1)/(r1*beta1);
    N_eq = (A*r1*beta1*r4*beta4-(beta2*(d+eps4)+r4*beta4*d)*(d+eps1))/(beta1*r4*beta4*(d+eps1));
    T_eq = (d+eps4)/(r4*beta4);
    V_eq = B*r4*beta4/(r3*beta3*(d+eps4)*(A3-1));
    Z_eq = (r3*(r2*beta2*(d+eps1)-r1*beta1*(d+eps2))*(d+eps4)*(A3-1)-B*r1*beta1*r4*beta4)/(r1*beta1*r3*beta4*(d+eps4)*(A3-1));
    break;
  }
  return 0;
}

#define EQ 1
static void equilibrium()
{
  printf("S_eq: %g\n", S_eq);
  printf("N_eq: %g\n", N_eq);
  printf("T_eq: %g\n", T_eq);
  printf("V_eq: %g\n", V_eq);
  printf("Z_eq: %g\n", Z_eq);

  printf("A1: %g\n", A1);
  printf("A2: %g\n", A2);
  printf("A3: %g\n", A3);
}

static void reaction(const void *ctx,
                     const ssdc_real_t U[],
                     /* */ ssdc_real_t Ut[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t N = U[1];
  const ssdc_real_t T = U[2];
  const ssdc_real_t V = U[3];
  const ssdc_real_t Z = U[4];

  Ut[0] = A - d*S - beta1*S*N - beta2*S*T;
  Ut[1] =     r1*beta1*S*N - (d+eps1)*N;
  Ut[2] =     r2*beta2*S*T - (d+eps2)*T - beta3*T*V - beta4*T*Z;
  Ut[3] = B + r3*beta3*T*V - (d+eps3)*V;
  Ut[4] =     r4*beta4*T*Z - (d+eps4)*Z;
}

static ssdc_real_t phi(ssdc_real_t n) { return n - Log(n) - 1; }

static void lyapunov(const void *ctx,
                     const ssdc_real_t U[],
                     const ssdc_real_t x[],
                     /**/  ssdc_real_t *L)
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t N = U[1];
  const ssdc_real_t T = U[2];
  const ssdc_real_t V = U[3];
  const ssdc_real_t Z = U[4];

  switch (CASE) {
  case 0:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N                * 1/r1;
    *L += T                * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z                * 1/(r2*r4);
    break;
  case 1:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N_eq*phi(N/N_eq) * 1/r1;
    *L += T                * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z                * 1/(r2*r4);
    break;
  case 2:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N                * 1/r1;
    *L += T_eq*phi(T/T_eq) * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z                * 1/(r2*r4);
    break;
  case 3:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N_eq*phi(N/N_eq) * 1/r1;
    *L += T_eq*phi(T/T_eq) * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z                * 1/(r2*r4);
    break;
  case 4:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N                * 1/r1;
    *L += T_eq*phi(T/T_eq) * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z_eq*phi(Z/Z_eq) * 1/(r2*r4);
    break;
  case 5:
    *L  = 0;
    *L += S_eq*phi(S/S_eq) * 1;
    *L += N_eq*phi(N/N_eq) * 1/r1;
    *L += T_eq*phi(T/T_eq) * 1/r2;
    *L += V_eq*phi(V/V_eq) * 1/(r2*r3);
    *L += Z_eq*phi(Z/Z_eq) * 1/(r2*r4);
    break;
  }
}

static void entropy(const void *ctx,
                    const ssdc_real_t U[],
                    /* */ ssdc_real_t W[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t N = U[1];
  const ssdc_real_t T = U[2];
  const ssdc_real_t V = U[3];
  const ssdc_real_t Z = U[4];

  switch (CASE) {
  case 0:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = 1            * 1/r1;
    W[2] = 1            * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = 1            * 1/(r2*r4);
    break;
  case 1:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = (1 - N_eq/N) * 1/r1;
    W[2] = 1            * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = 1            * 1/(r2*r4);
    break;
  case 2:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = 1            * 1/r1;
    W[2] = (1 - T_eq/T) * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = 1            * 1/(r2*r4);
    break;
  case 3:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = (1 - N_eq/N) * 1/r1;
    W[2] = (1 - T_eq/T) * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = 1            * 1/(r2*r4);
    break;
  case 4:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = 1            * 1/r1;
    W[2] = (1 - T_eq/T) * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = (1 - Z_eq/Z) * 1/(r2*r4);
    break;
  case 5:
    W[0] = (1 - S_eq/S) * 1;
    W[1] = (1 - N_eq/N) * 1/r1;
    W[2] = (1 - T_eq/T) * 1/r2;
    W[3] = (1 - V_eq/V) * 1/(r2*r3);
    W[4] = (1 - Z_eq/Z) * 1/(r2*r4);
    break;
  }
}

PETSC_UNUSED
static void jacobian(const ssdc_real_t U[],
                     /* */ ssdc_real_t dUdW[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t N = U[1];
  const ssdc_real_t T = U[2];
  const ssdc_real_t V = U[3];
  const ssdc_real_t Z = U[4];

  switch (CASE) {
  case 0:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = 0           * r1;
    dUdW[2] = 0           * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = 0           * (r2*r4);
    break;
  case 1:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = Sqr(N)/N_eq * r1;
    dUdW[2] = 0           * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = 0           * (r2*r4);
    break;
  case 2:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = 0           * r1;
    dUdW[2] = Sqr(T)/T_eq * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = 0           * (r2*r4);
    break;
  case 3:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = Sqr(N)/N_eq * r1;
    dUdW[2] = Sqr(T)/T_eq * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = 0           * (r2*r4);
    break;
  case 4:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = 0           * r1;
    dUdW[2] = Sqr(T)/T_eq * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = Sqr(Z)/Z_eq * (r2*r4);
    break;
  case 5:
    dUdW[0] = Sqr(S)/S_eq * 1;
    dUdW[1] = Sqr(N)/N_eq * r1;
    dUdW[2] = Sqr(T)/T_eq * r2;
    dUdW[3] = Sqr(V)/V_eq * (r2*r3);
    dUdW[4] = Sqr(Z)/Z_eq * (r2*r4);
    break;
  }
}

static ssdc_real_t T_final = 400;

static void initial(const void *ctx,
                    const ssdc_real_t t,
                    const ssdc_real_t xyz[],
                    /* */ ssdc_real_t U[])
{
  ssdc_real_t S_ini = 0.30;
  ssdc_real_t N_ini = 0.20;
  ssdc_real_t T_ini = 0.10;
  ssdc_real_t V_ini = 0.10;
  ssdc_real_t Z_ini = 0.01;


  const ssdc_real_t pi = PETSC_PI;
  const ssdc_real_t x = xyz[0];
  const ssdc_real_t cos_pi_x = PetscCosReal(pi*x);
  const ssdc_real_t wave = 1 + 0.2 * PetscSqr(cos_pi_x);

  S_ini *= wave;
  N_ini *= wave;
  T_ini *= wave;
  V_ini *= wave;
  Z_ini *= wave;


  U[0] = S_ini;
  U[1] = N_ini;
  U[2] = T_ini;
  U[3] = V_ini;
  U[4] = Z_ini;
}

#define HAVE_CONVERGENCE 1
static void convergence(ssdc_real_t U[])
{
  U[0] = S_eq;
  U[1] = N_eq;
  U[2] = T_eq;
  U[3] = V_eq;
  U[4] = Z_eq;
}

#include "model-main.c.in"
