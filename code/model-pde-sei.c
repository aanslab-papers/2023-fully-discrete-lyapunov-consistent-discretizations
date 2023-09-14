#include <ssdc.h>
#include <assert.h>
#include "model-math.h"

#define NAME "SEI"
enum {dof = 3};

#ifndef HAVE_DIFFUSION
#define HAVE_DIFFUSION 1
#endif

#if HAVE_DIFFUSION == 1
#define OPT_DEG "2"
#define OPT_NEL "60"
#define OPT_BOX "0,50,0,50"
#endif

static ssdc_real_t g(ssdc_real_t y) { return y/(y+1); }
static ssdc_real_t beta    = 0;
static ssdc_real_t Lambda1 = 0;
static ssdc_real_t Lambda2 = 0;
static ssdc_real_t Lambda3 = 0;
static ssdc_real_t mu1     = 0;
static ssdc_real_t mu2     = 0;
static ssdc_real_t mu3     = 0;
static ssdc_real_t d1      = 0;
static ssdc_real_t d2      = 0;
static ssdc_real_t d3      = 0;
static ssdc_real_t Cdiff[dof][dof];

static ssdc_real_t S0 = 50;
static ssdc_real_t E0 = 15;
static ssdc_real_t I0 = 10;
static ssdc_real_t N;

static ssdc_real_t S_eq_e;
static ssdc_real_t E_eq_e;
static ssdc_real_t I_eq_e;
static ssdc_real_t B;

static void reaction(const void *,const ssdc_real_t[],ssdc_real_t[]);

#define ParamInit(name) \
  do {name = SSDCOptionReal("-"#name, name);} while(0)

PetscErrorCode ModelInitialize(void)
{
  const int example = (int)SSDCOptionInt("-example", 1);
  switch(example) {
  case 1:
    Lambda1 = 1.0;
    Lambda2 = 1.2;
    Lambda3 = 0.95;
    mu1     = 0.2;
    mu2     = 0.3;
    mu3     = 0.45;
    beta    = 0.02;
    d1      = 1;
    d2      = 1.5;
    d3      = 1.3;
    break;
  case 2:
    Lambda1 = 1.5;
    Lambda2 = 1.2;
    Lambda3 = 0.95;
    mu1     = 0.2;
    mu2     = 0.3;
    mu3     = 0.45;
    beta    = 0.005;
    d1      = 150;
    d2      = 100;
    d3      = 200;
    break;
  default:
    break;
  }

  ParamInit( beta    );
  ParamInit( Lambda1 );
  ParamInit( Lambda2 );
  ParamInit( Lambda3 );
  ParamInit( mu1     );
  ParamInit( mu2     );
  ParamInit( mu3     );

  ParamInit( d1 );
  ParamInit( d2 );
  ParamInit( d3 );
  Cdiff[0][0] = d1;
  Cdiff[1][1] = d2;
  Cdiff[2][2] = d3;

  ParamInit( S0 );
  ParamInit( E0 );
  ParamInit( I0 );
  N = S0 + E0 + I0;


  S_eq_e =  (1.0/2.0)*(2*Lambda1*beta*mu1 + Lambda1*beta + Lambda2*beta*mu1 + Lambda2*beta + Lambda3*beta*mu1 + Lambda3*beta + Lambda3*mu1*mu2 + Lambda3*mu2 + beta*mu1*mu3 + mu1*mu2*mu3 - sqrt(Pow(Lambda1, 2)*Pow(beta, 2) + 2*Lambda1*Lambda2*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda2*Pow(beta, 2) + 2*Lambda1*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda3*Pow(beta, 2) + 2*Lambda1*Lambda3*beta*mu1*mu2 + 2*Lambda1*Lambda3*beta*mu2 - 2*Lambda1*Pow(beta, 2)*mu1*mu3 - 2*Lambda1*beta*mu1*mu2*mu3 + Pow(Lambda2, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda2, 2)*Pow(beta, 2)*mu1 + Pow(Lambda2, 2)*Pow(beta, 2) + 2*Lambda2*Lambda3*Pow(beta, 2)*Pow(mu1, 2) + 4*Lambda2*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda2*Lambda3*Pow(beta, 2) + 2*Lambda2*Lambda3*beta*Pow(mu1, 2)*mu2 + 4*Lambda2*Lambda3*beta*mu1*mu2 + 2*Lambda2*Lambda3*beta*mu2 + 2*Lambda2*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda2*Pow(beta, 2)*mu1*mu3 + 2*Lambda2*beta*Pow(mu1, 2)*mu2*mu3 + 2*Lambda2*beta*mu1*mu2*mu3 + Pow(Lambda3, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda3, 2)*Pow(beta, 2)*mu1 + Pow(Lambda3, 2)*Pow(beta, 2) + 2*Pow(Lambda3, 2)*beta*Pow(mu1, 2)*mu2 + 4*Pow(Lambda3, 2)*beta*mu1*mu2 + 2*Pow(Lambda3, 2)*beta*mu2 + Pow(Lambda3, 2)*Pow(mu1, 2)*Pow(mu2, 2) + 2*Pow(Lambda3, 2)*mu1*Pow(mu2, 2) + Pow(Lambda3, 2)*Pow(mu2, 2) + 2*Lambda3*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda3*Pow(beta, 2)*mu1*mu3 + 4*Lambda3*beta*Pow(mu1, 2)*mu2*mu3 + 4*Lambda3*beta*mu1*mu2*mu3 + 2*Lambda3*Pow(mu1, 2)*Pow(mu2, 2)*mu3 + 2*Lambda3*mu1*Pow(mu2, 2)*mu3 + Pow(beta, 2)*Pow(mu1, 2)*Pow(mu3, 2) + 2*beta*Pow(mu1, 2)*mu2*Pow(mu3, 2) + Pow(mu1, 2)*Pow(mu2, 2)*Pow(mu3, 2)))/(beta*mu1*(mu1 + 1)) ;

  E_eq_e =  (1.0/2.0)*(Lambda1*beta + Lambda2*beta*mu1 + Lambda2*beta - Lambda3*beta*mu1 - Lambda3*beta - Lambda3*mu1*mu2 - Lambda3*mu2 - beta*mu1*mu3 - mu1*mu2*mu3 + sqrt(Pow(Lambda1, 2)*Pow(beta, 2) + 2*Lambda1*Lambda2*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda2*Pow(beta, 2) + 2*Lambda1*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda3*Pow(beta, 2) + 2*Lambda1*Lambda3*beta*mu1*mu2 + 2*Lambda1*Lambda3*beta*mu2 - 2*Lambda1*Pow(beta, 2)*mu1*mu3 - 2*Lambda1*beta*mu1*mu2*mu3 + Pow(Lambda2, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda2, 2)*Pow(beta, 2)*mu1 + Pow(Lambda2, 2)*Pow(beta, 2) + 2*Lambda2*Lambda3*Pow(beta, 2)*Pow(mu1, 2) + 4*Lambda2*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda2*Lambda3*Pow(beta, 2) + 2*Lambda2*Lambda3*beta*Pow(mu1, 2)*mu2 + 4*Lambda2*Lambda3*beta*mu1*mu2 + 2*Lambda2*Lambda3*beta*mu2 + 2*Lambda2*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda2*Pow(beta, 2)*mu1*mu3 + 2*Lambda2*beta*Pow(mu1, 2)*mu2*mu3 + 2*Lambda2*beta*mu1*mu2*mu3 + Pow(Lambda3, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda3, 2)*Pow(beta, 2)*mu1 + Pow(Lambda3, 2)*Pow(beta, 2) + 2*Pow(Lambda3, 2)*beta*Pow(mu1, 2)*mu2 + 4*Pow(Lambda3, 2)*beta*mu1*mu2 + 2*Pow(Lambda3, 2)*beta*mu2 + Pow(Lambda3, 2)*Pow(mu1, 2)*Pow(mu2, 2) + 2*Pow(Lambda3, 2)*mu1*Pow(mu2, 2) + Pow(Lambda3, 2)*Pow(mu2, 2) + 2*Lambda3*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda3*Pow(beta, 2)*mu1*mu3 + 4*Lambda3*beta*Pow(mu1, 2)*mu2*mu3 + 4*Lambda3*beta*mu1*mu2*mu3 + 2*Lambda3*Pow(mu1, 2)*Pow(mu2, 2)*mu3 + 2*Lambda3*mu1*Pow(mu2, 2)*mu3 + Pow(beta, 2)*Pow(mu1, 2)*Pow(mu3, 2) + 2*beta*Pow(mu1, 2)*mu2*Pow(mu3, 2) + Pow(mu1, 2)*Pow(mu2, 2)*Pow(mu3, 2)))/(beta*(beta*mu1 + beta + mu1*mu2 + mu2)) ;

  I_eq_e =  (1.0/2.0)*(Lambda1*beta + Lambda2*beta*mu1 + Lambda2*beta + Lambda3*beta*mu1 + Lambda3*beta + Lambda3*mu1*mu2 + Lambda3*mu2 - beta*mu1*mu3 - mu1*mu2*mu3 + sqrt(Pow(Lambda1, 2)*Pow(beta, 2) + 2*Lambda1*Lambda2*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda2*Pow(beta, 2) + 2*Lambda1*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda1*Lambda3*Pow(beta, 2) + 2*Lambda1*Lambda3*beta*mu1*mu2 + 2*Lambda1*Lambda3*beta*mu2 - 2*Lambda1*Pow(beta, 2)*mu1*mu3 - 2*Lambda1*beta*mu1*mu2*mu3 + Pow(Lambda2, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda2, 2)*Pow(beta, 2)*mu1 + Pow(Lambda2, 2)*Pow(beta, 2) + 2*Lambda2*Lambda3*Pow(beta, 2)*Pow(mu1, 2) + 4*Lambda2*Lambda3*Pow(beta, 2)*mu1 + 2*Lambda2*Lambda3*Pow(beta, 2) + 2*Lambda2*Lambda3*beta*Pow(mu1, 2)*mu2 + 4*Lambda2*Lambda3*beta*mu1*mu2 + 2*Lambda2*Lambda3*beta*mu2 + 2*Lambda2*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda2*Pow(beta, 2)*mu1*mu3 + 2*Lambda2*beta*Pow(mu1, 2)*mu2*mu3 + 2*Lambda2*beta*mu1*mu2*mu3 + Pow(Lambda3, 2)*Pow(beta, 2)*Pow(mu1, 2) + 2*Pow(Lambda3, 2)*Pow(beta, 2)*mu1 + Pow(Lambda3, 2)*Pow(beta, 2) + 2*Pow(Lambda3, 2)*beta*Pow(mu1, 2)*mu2 + 4*Pow(Lambda3, 2)*beta*mu1*mu2 + 2*Pow(Lambda3, 2)*beta*mu2 + Pow(Lambda3, 2)*Pow(mu1, 2)*Pow(mu2, 2) + 2*Pow(Lambda3, 2)*mu1*Pow(mu2, 2) + Pow(Lambda3, 2)*Pow(mu2, 2) + 2*Lambda3*Pow(beta, 2)*Pow(mu1, 2)*mu3 + 2*Lambda3*Pow(beta, 2)*mu1*mu3 + 4*Lambda3*beta*Pow(mu1, 2)*mu2*mu3 + 4*Lambda3*beta*mu1*mu2*mu3 + 2*Lambda3*Pow(mu1, 2)*Pow(mu2, 2)*mu3 + 2*Lambda3*mu1*Pow(mu2, 2)*mu3 + Pow(beta, 2)*Pow(mu1, 2)*Pow(mu3, 2) + 2*beta*Pow(mu1, 2)*mu2*Pow(mu3, 2) + Pow(mu1, 2)*Pow(mu2, 2)*Pow(mu3, 2)))/(mu3*(beta*mu1 + beta + mu1*mu2 + mu2)) ;

  B = S_eq_e*g(I_eq_e)/(beta*E_eq_e);

  {
    ssdc_real_t U[dof],F[dof];
    ssdc_real_t rtol=0,atol=1000*PETSC_MACHINE_EPSILON;
    U[0] = S_eq_e; U[1] = E_eq_e; U[2] = I_eq_e; reaction(NULL,U,F);
    assert(PetscIsCloseAtTol(F[0],0,rtol,atol));
    assert(PetscIsCloseAtTol(F[1],0,rtol,atol));
    assert(PetscIsCloseAtTol(F[2],0,rtol,atol));
  }
  return 0;
}

#define EQ 1
static void equilibrium()
{
  PetscPrintf(PETSC_COMM_WORLD,"S_eq_e: %g\n", S_eq_e);
  PetscPrintf(PETSC_COMM_WORLD,"E_eq_e: %g\n", E_eq_e);
  PetscPrintf(PETSC_COMM_WORLD,"I_eq_e: %g\n", I_eq_e);
}

static void reaction(const void *ctx,
                     const ssdc_real_t U[],
                     /* */ ssdc_real_t Ut[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t E = U[1];
  const ssdc_real_t I = U[2];
  Ut[0] = Lambda1 - S*g(I) - mu1*S;
  Ut[1] = Lambda2 + S*g(I) - (mu2+beta)*E;
  Ut[2] = Lambda3 + beta*E - mu3*I;
}

static void lyapunov(const void *ctx,
                     const ssdc_real_t U[],
                     const ssdc_real_t x[],
                     /**/  ssdc_real_t *L)
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t E = U[1];
  const ssdc_real_t I = U[2];
  *L  = 0;
  *L +=   (S - S_eq_e - S_eq_e*Log(S/S_eq_e));
  *L +=   (E - E_eq_e - E_eq_e*Log(E/E_eq_e));
  *L += B*(I - I_eq_e - I_eq_e*Log(I/I_eq_e));
}

static void entropy(const void *ctx,
                    const ssdc_real_t U[],
                    /* */ ssdc_real_t W[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t E = U[1];
  const ssdc_real_t I = U[2];
  W[0] =   (1 - S_eq_e/S);
  W[1] =   (1 - E_eq_e/E);
  W[2] = B*(1 - I_eq_e/I);
}

PETSC_UNUSED
static void jacobian(const ssdc_real_t U[],
                     /* */ ssdc_real_t dUdW[])
{
  const ssdc_real_t S = U[0];
  const ssdc_real_t E = U[1];
  const ssdc_real_t I = U[2];
  dUdW[0] = Sqr(S)/S_eq_e;
  dUdW[1] = Sqr(E)/E_eq_e;
  dUdW[2] = Sqr(I)/(B*I_eq_e);
}

static ssdc_real_t T_final = 40;

static void initial(const void *ctx,
                        const ssdc_real_t t,
                        const ssdc_real_t x[],
                        /* */ ssdc_real_t U[])
{
  U[0] = S0;
  U[1] = E0;
  U[2] = I0;
#if HAVE_DIFFUSION == 1
  {
    const ssdc_real_t wx = Sin(Pi*(x[0]-0.1)/0.2);
    const ssdc_real_t wy = Sin(Pi*(x[1]-0.1)/0.2);
    const ssdc_real_t w  = Sqr(wx)*Sqr(wy);
    U[0] = 100 - 0.50*w;
    U[1] =   3 - 0.25*w;
    U[2] =  30 - 0.25*w;
  }
#endif
}

#define HAVE_CONVERGENCE 1
static void convergence(ssdc_real_t U[])
{
  U[0] = S_eq_e;
  U[1] = E_eq_e;
  U[2] = I_eq_e;
}

#include "model-main.c.in"
