#undef  I            /* complex.h */
#undef  gamma        /* math.h    */
#define gamma Gamma  /* math.h    */

#define Pi   PETSC_PI
#define Sin  PetscSinReal
#define Cos  PetscCosReal

#define Exp  PetscExpReal
#define Log  PetscLogReal
#define Pow  PetscPowReal
#define Abs  PetscAbsReal

#define Sqr  PetscSqr
#define Sqrt PetscSqrtReal

