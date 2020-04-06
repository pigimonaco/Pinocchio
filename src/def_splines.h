/* #if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC) */
/* #error Trying to compile with ELL_CLASSIC and SCALE_DEPENDENT together */
/* #endif */

#if defined(MOD_GRAV_FR) && !defined(SCALE_DEPENDENT)
#define SCALE_DEPENDENT
#endif

#ifndef SCALE_DEPENDENT
#define NkBINS 1
#else
#define NkBINS 10
#define LOGKMIN ((double)-3.0)
#define DELTALOGK ((double)0.5)
#endif

#define SP_TIME     0
#define SP_INVTIME  1
#define SP_COMVDIST 2
#define SP_DIAMDIST 3
#define SP_INVGROW  4
#define SP_MASSVAR  5
#define SP_DISPVAR  6
#define SP_RADIUS   7
#define SP_DVARDR   8

#define SP_GROW1    (9)
#define SP_GROW2    (9+NkBINS)
#define SP_GROW31   (9+2*NkBINS)
#define SP_GROW32   (9+3*NkBINS)
#define SP_FOMEGA1  (9+4*NkBINS)
#define SP_FOMEGA2  (9+5*NkBINS)
#define SP_FOMEGA31 (9+6*NkBINS)
#define SP_FOMEGA32 (9+7*NkBINS)

#define SP_PK       (9+8*NkBINS)
#define SP_EOS      (9+8*NkBINS+1)
#define SP_INTEOS   (9+8*NkBINS+2)

#define NSPLINES (9+8*NkBINS+3)
