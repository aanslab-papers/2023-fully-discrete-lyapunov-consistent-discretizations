#include <ssdc.h>
#include <assert.h>
#include <stdlib.h>
#include "model-math.h"

#define NAME "SCALED-SI"
enum {dof = 2};

#ifndef HAVE_DIFFUSION
#define HAVE_DIFFUSION 1
#endif

#if HAVE_DIFFUSION == 1
#define OPT_DEG "2"
#define OPT_NEL "60"
#define OPT_BOX "0,20,0,20"
#endif


static ssdc_real_t R0    = 0.0;
static ssdc_real_t Rd    = 2.0;
static ssdc_real_t nu    = 0.15;

static ssdc_real_t d1   = 0.01;
static ssdc_real_t d2   = 0.25;
static ssdc_real_t Cdiff[dof][dof];

static ssdc_real_t S0 = 0;
static ssdc_real_t I0 = 0;
static ssdc_real_t N;

static ssdc_real_t lambda;
static ssdc_real_t S_eq_e;
static ssdc_real_t I_eq_e;

#define ParamInit(name) \
  do {name = SSDCOptionReal("-"#name, name);} while(0)

PetscErrorCode ModelInitialize(void)
{
  const int example = (int)SSDCOptionInt("-example", 1);
  switch(example) {
  case 1:
    R0 = 1.14;
    break;
  case 2:
    R0 = 1.154;
    break;
  case 3:
    R0 = 1.17;
    break;
  case 4:
    R0 = 1.2;
    break;
  case 5:
    R0 = 1.213;
    break;
  default:
    break;
  }

  ParamInit( R0    );
  ParamInit( Rd    );
  ParamInit( nu );

  ParamInit( d1 );
  ParamInit( d2 );

  Cdiff[0][0] = d1;
  Cdiff[1][1] = d2;

  ParamInit( S0 );
  ParamInit( I0 );
  N = S0 + I0;

  #if DISEASE_FREE == 1
  S_eq_e = 1-1/Rd;
  I_eq_e = 1e-6;
  #else
  S_eq_e = (nu*R0*Rd-R0+1-nu)/(nu*R0*R0*Rd);
  I_eq_e = S_eq_e*(R0-1);
  #endif
  lambda = S_eq_e/I_eq_e;

  return 0;
}

#define EQ 1
static void equilibrium()
{
  PetscPrintf(PETSC_COMM_WORLD,"S_eq_e: %.10g\n", S_eq_e);
  PetscPrintf(PETSC_COMM_WORLD,"I_eq_e: %.10g\n", I_eq_e);
}

static void reaction(const void *ctx,
                     const ssdc_real_t U[],
                     /* */ ssdc_real_t Ut[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t I = U[1];
  Ut[0] = nu*Rd*(S+I)*(1-(S+I))-R0*S*I/(S+I)-nu*S;
  Ut[1] = R0*S*I/(S+I)-I;
}

static void lyapunov(const void *ctx,
                     const ssdc_real_t U[],
                     const ssdc_real_t x[],
                     /**/  ssdc_real_t *L)
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t I = U[1];
  *L  = 0;
  *L +=        (S - S_eq_e - S_eq_e*Log(S/S_eq_e));
  *L += lambda*(I - I_eq_e - I_eq_e*Log(I/I_eq_e));
}

static void entropy(const void *ctx,
                    const ssdc_real_t U[],
                    /* */ ssdc_real_t W[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t I = U[1];
  W[0] =        (1 - S_eq_e/S);
  W[1] = lambda*(1 - I_eq_e/I);
}

PETSC_UNUSED
static void jacobian(const ssdc_real_t U[],
                     /* */ ssdc_real_t dUdW[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t I = U[1];
  dUdW[0] =          Sqr(S)/S_eq_e;
  dUdW[1] = 1/lambda*Sqr(I)/I_eq_e;
}

ssdc_real_t mkrandom()
{
  return (ssdc_real_t)rand()/(ssdc_real_t)RAND_MAX;
}

static ssdc_real_t T_final = 2; 
// T_final=3750

static void initial(const void *ctx,
                        const ssdc_real_t t,
                        const ssdc_real_t x[],
                        /* */ ssdc_real_t U[])
{

  ssdc_real_t eps = 1e-4;
  N  = S_eq_e + I_eq_e;
  S0 = S_eq_e + eps * mkrandom();
  I0 = N - S0;

  U[0] = S0;
  U[1] = I0;
#if HAVE_DIFFUSION == 1
  {
    U[0] = S0;
    U[1] = I0;
  }
#endif
}

#define HAVE_CONVERGENCE 1
static void convergence(ssdc_real_t U[])
{
  U[0] = S_eq_e;
  U[1] = I_eq_e;
}

#include "model-main.c.in"
