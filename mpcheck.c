/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

/* To adapt to the tested precision */
#define FPPREC 53       /* 24       53         64      113 */

#if (FPPREC <= 24)
#define fptype float
#define EMIN -125
#define EMAX 128
#elif (FPPREC <= 53)
#define fptype double
#define EMIN -1021
#define EMAX 1024
#elif (FPPREC <= 64)
#if defined(__ia64)
#define fptype extended
#define _FPWIDETYPES /* for HP-UX */
#else
#define fptype long double
#endif
#define EMIN -16381
#define EMAX 16384
#elif (FPPREC <= 113)
#if defined(__ia64)
#define _FPWIDETYPES /* for HP-UX */
#define fptype long double
#else
#define fptype quad
#endif
#define EMIN -16381
#define EMAX 16384
#else
#error "not yet implemented"
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"
#include "mpfr.h"

/* tgamma is not defined in some versions of math.h */
#ifndef tgamma
fptype tgamma _PROTO((fptype));
#endif
#ifndef log2
fptype log2 _PROTO((fptype));
#endif
#ifndef exp2
fptype exp2 _PROTO((fptype));
#endif

/* stolen from mpfr-impl.h */
#define MPFR_EXP(x) ((x)->_mpfr_exp)

#ifdef __i386
#include <fpu_control.h>
#ifndef __setfpucw
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif /* ifndef __setfpucw */
#endif

#ifdef __alpha
#ifndef FP_RND_RN
#define FP_RND_RZ       0
#define FP_RND_RN       1
#define FP_RND_RP       2
#define FP_RND_RM       3
#endif
#define TONEAREST write_rnd(FP_RND_RN)
#define TOZERO    write_rnd(FP_RND_RZ)
#define TOINFP    write_rnd(FP_RND_RP)
#define TOINFM    write_rnd(FP_RND_RM)
#else /* ifdef __alpha */
#include <fenv.h>
#define TONEAREST fesetround(FE_TONEAREST)
#define TOZERO    fesetround(FE_TOWARDZERO)
#define TOINFP    fesetround(FE_UPWARD)
#define TOINFM    fesetround(FE_DOWNWARD)
#endif

#ifndef MAX_RND
#define MAX_RND 4
#endif

/* useful for testing monotonicity and symmetry */
#define NO_MONOTON 0
#define INCREASING 1
#define DECREASING -1


#define NO_SYMM 0
#define ODD 3
#define EVEN 4


void test _PROTO ((char *, mp_exp_t, unsigned long, unsigned long));
void test2 _PROTO ((char *, mp_exp_t, mp_exp_t, unsigned long, unsigned long));
double ulp_err _PROTO ((fptype, fptype, mp_rnd_t, fptype));
double ulp_err2 _PROTO ((fptype, fptype, fptype, mp_rnd_t, fptype));
void testall _PROTO ((unsigned long, unsigned long));
void check_fp _PROTO ((void));

/* sets the machine rounding mode to the value rnd_mode */
fptype (*testfun_libm) (fptype) = NULL;
fptype (*testfun_libm2) (fptype, fptype) = NULL;
int (*testfun_mpfr) () = NULL;
int (*mpfr_set_fp) (mpfr_ptr, fptype, mp_rnd_t) = NULL;
fptype (*mpfr_get_fp) (mpfr_srcptr, mp_rnd_t) = NULL;
void (*fprint_fp) (FILE *, fptype) = NULL;
double MAX_ERR_NEAR = 0.0;
double MAX_ERR_DIR = 0.0;
int verbose = 3;
int test_monotonicity = 1;
int test_range = 1;
int test_symmetry = 1;

#define print_fp(x) fprint_fp(stdout, x)

void 
mpfr_set_machine_rnd_mode (mp_rnd_t rnd_mode)
{
  switch (rnd_mode) {
  case GMP_RNDN: TONEAREST; break;
  case GMP_RNDZ: TOZERO; break;
  case GMP_RNDU: TOINFP; break;
  case GMP_RNDD: TOINFM; break;
  default: fprintf(stderr, "invalid rounding mode\n"); exit(1);
  }
}

void
fprint_d (FILE *stream, fptype x)
{
  fprintf (stream, "%1.20e", x);
}

void
fprint_ld (FILE *stream, fptype x)
{
#if (FPPREC <= 64)
  fprintf (stream, "%1.24Le", x);
#else
  fprintf (stream, "%1.38Le", x);
#endif
}

void
check_fp ()
{
  fptype x, y, c, d, dj;
  int j;

  /* sets the functions mpfr_set_fp and mpfr_get_fp */
#if (FPPREC <= 53)
  mpfr_set_fp = mpfr_set_d;
  mpfr_get_fp = mpfr_get_d;
  fprint_fp = fprint_d;
#elif (FPPREC <= 64)
  mpfr_set_fp = mpfr_set_ld;
  mpfr_get_fp = mpfr_get_ld;
  fprint_fp = fprint_ld;
#elif (FPPREC <= 113)
  mpfr_set_fp = mpfr_set_ld;
  mpfr_get_fp = mpfr_get_ld;
  fprint_fp = fprint_ld;
#endif

#if (FPPREC <= 53)
#if defined(__i386)
  /* sets the precision to double */
  __setfpucw((_FPU_DEFAULT & (~_FPU_EXTENDED)) | _FPU_DOUBLE);
#endif
  x = DBL_MIN;
  if (2.0 * (x / 2.0) != x)
    fprintf (stderr, "Warning: no denormalized numbers\n");

  c = 1.46484375e-3;
  dj = 1.0;
  for (j=0; j<54; j++) dj *= 0.5;
  d = 0.75 + dj;
  d /= 1 << 9;
  if (c != d)
    {
      fprintf (stderr, "default seems to use extended precision\n");
      exit (1);
    }
#endif

  /* checks that setting machine rounding mode works */
  mpfr_set_machine_rnd_mode (GMP_RNDD);
  x = 2.0; /* don't write x = sqrt (2.0) and y = sqrt (2.0) otherwise the
              compiler may optimize too much */
  x = sqrt (x);
  mpfr_set_machine_rnd_mode (GMP_RNDU);
  y = 2.0;
  y = sqrt (y);
  if (x == y)
    {
      fprintf (stderr, "setting rounding modes does not work, "
               "you may have to add an option to the C compiler\n");
      fprintf (stderr, "   on Alpha, try:\n");
      fprintf (stderr, "   '-fprm d -ieee_with_inexact' with cc\n");
      fprintf (stderr, "   '-mfp-rounding-mode=d -mieee-with-inexact' with gcc\n");
      exit (1);
    }

  /* checks that the precision of fptype is FPPREC */
  mpfr_set_machine_rnd_mode (GMP_RNDN);
  x = 1.0;
  for (j=0; x + 1.0 != x; j++, x = 2.0 * x);
  if (j != FPPREC)
    {
      fprintf (stderr, "Precision of fptype is not %u but %u\n", FPPREC, j);
      exit (1);
    }
}

/* computes the error in ulps between f(x) and y [computed by libm].
   z is the value computed by mpfr. It assumes y <> z.
 */
double ulp_err (fptype x, fptype y, mp_rnd_t rnd, fptype z)
{
  mpfr_t xx, yy;
  mp_exp_t expy;
  double u;
  
  mpfr_init2 (xx, FPPREC);
  mpfr_init2 (yy, 2 * FPPREC);
  mpfr_set_fp (xx, x, GMP_RNDN); /* exact */
  testfun_mpfr (yy, xx, rnd);
  
  /* checks that the mpfr value z is correct */
  if (mpfr_can_round (yy, 2 * FPPREC, rnd, rnd, FPPREC))
    {
      mpfr_set (xx, yy, rnd);
      if (mpfr_get_fp (xx, GMP_RNDN) != z)
        {
          fprintf (stderr, "Error in mpfr !!!\n");
          fprintf (stderr, "x=");
          fprint_fp (stderr, x);
          fprintf (stderr, " rnd=%s\n", mpfr_print_rnd_mode (rnd));
          fprintf (stderr, "mpfr got   ");
          fprint_fp (stderr, z);
          fprintf (stderr, "\ninstead of ");
          fprint_fp (stderr, mpfr_get_fp (xx, GMP_RNDN));
          fprintf (stderr, "\n");
          exit (1);
        }
    }

  expy = MPFR_EXP (yy);
  mpfr_set_fp (xx, y, GMP_RNDN); /* exact */
  mpfr_sub (xx, xx, yy, GMP_RNDN);
  MPFR_EXP(xx) += FPPREC - expy;
  u = mpfr_get_d (xx, GMP_RNDN);
  mpfr_clear (xx);
  mpfr_clear (yy);
  return u;
}

/* computes the error in ulps between f(x,t) and y [computed by libm],
   where z is the value computed by mpfr.
 */
double ulp_err2 (fptype x, fptype t, fptype y, mp_rnd_t rnd, fptype z)
{
  mpfr_t xx, yy, tt;
  mp_exp_t expy;
  double u;
  
  mpfr_init2 (xx, FPPREC);
  mpfr_init2 (tt, FPPREC);
  mpfr_init2 (yy, 2 * FPPREC);
  mpfr_set_fp (xx, x, GMP_RNDN); /* exact */
  mpfr_set_fp (tt, t, GMP_RNDN); /* exact */
  testfun_mpfr (yy, xx, tt, rnd);

  /* checks that the mpfr value z is correct */
  if (mpfr_can_round (yy, 2 * FPPREC, rnd, rnd, FPPREC))
    {
      mpfr_set (xx, yy, rnd);
      if (mpfr_get_fp (xx, GMP_RNDN) != z)
        {
          fprintf (stderr, "Error in mpfr !!!\n");
          fprintf (stderr, "x=");
          fprint_fp (stderr, x);
          fprintf (stderr, " t=");
          fprint_fp (stderr, t);
          fprintf (stderr, " rnd=%s\n", mpfr_print_rnd_mode (rnd));
          fprintf (stderr, "mpfr got   ");
          fprint_fp (stderr, z);
          fprintf (stderr, "\ninstead of ");
          fprint_fp (stderr, mpfr_get_fp (xx, GMP_RNDN));
          fprintf (stderr, "\n");
          exit (1);
        }
    }

  expy = MPFR_EXP (yy);
  mpfr_set_fp (xx, y, GMP_RNDN); /* exact */
  mpfr_sub (xx, xx, yy, GMP_RNDN);
  MPFR_EXP(xx) += FPPREC - expy;
  u = mpfr_get_d (xx, GMP_RNDN);
  mpfr_clear (xx);
  mpfr_clear (tt);
  mpfr_clear (yy);
  return u;
}

#if (FPPREC==113 && defined(__ia64))
#define SUFFIXL
#endif

fptype my_exp (fptype x)
{
#ifdef SUFFIXL
  return expl (x);
#else
  return exp (x);
#endif
}

fptype my_exp2 (fptype x)
{
#ifdef SUFFIXL
  return exp2l (x);
#else
  return exp2 (x);
#endif
}

fptype my_expm1 (fptype x)
{
#ifdef SUFFIXL
  return expm1l (x);
#else
  return expm1 (x);
#endif
}

fptype my_log (fptype x)
{
#ifdef SUFFIXL
  return logl (x);
#else
  return log (x);
#endif
}

fptype my_log2 (fptype x)
{
#ifdef SUFFIXL
  return log2l (x);
#else
  return log2 (x);
#endif
}

fptype my_log10 (fptype x)
{
#ifdef SUFFIXL
  return log10l (x);
#else
  return log10 (x);
#endif
}

fptype my_log1p (fptype x)
{
#ifdef SUFFIXL
  return log1pl (x);
#else
  return log1p (x);
#endif
}

fptype my_sin (fptype x)
{
#ifdef SUFFIXL
  return sinl (x);
#else
  return sin (x);
#endif
}

fptype my_cos (fptype x)
{
#ifdef SUFFIXL
  return cosl (x);
#else
  return cos (x);
#endif
}

fptype my_tan (fptype x)
{
#ifdef SUFFIXL
  return tanl (x);
#else
  return tan (x);
#endif
}

fptype my_asin (fptype x)
{
#ifdef SUFFIXL
  return asinl (x);
#else
  return asin (x);
#endif
}

fptype my_acos (fptype x)
{
#ifdef SUFFIXL
  return acosl (x);
#else
  return acos (x);
#endif
}

fptype my_atan (fptype x)
{
#ifdef SUFFIXL
  return atanl (x);
#else
  return atan (x);
#endif
}

fptype my_sinh (fptype x)
{
#ifdef SUFFIXL
  return sinhl (x);
#else
  return sinh (x);
#endif
}

fptype my_cosh (fptype x)
{
#ifdef SUFFIXL
  return coshl (x);
#else
  return cosh (x);
#endif
}

fptype my_tanh (fptype x)
{
#ifdef SUFFIXL
  return tanhl (x);
#else
  return tanh (x);
#endif
}

fptype my_asinh (fptype x)
{
#ifdef SUFFIXL
  return asinhl (x);
#else
  return asinh (x);
#endif
}

fptype my_acosh (fptype x)
{
#ifdef SUFFIXL
  return acoshl (x);
#else
  return acosh (x);
#endif
}

fptype my_atanh (fptype x)
{
#ifdef SUFFIXL
  return atanhl (x);
#else
  return atanh (x);
#endif
}

fptype my_tgamma (fptype x)
{
#ifdef SUFFIXL
  return tgammal (x);
#else
  return tgamma (x);
#endif
}

fptype my_sqrt (fptype x)
{
#ifdef SUFFIXL
  return sqrtl (x);
#else
  return sqrt (x);
#endif
}

fptype my_cbrt (fptype x)
{
#ifdef SUFFIXL
  return cbrtl (x);
#else
  return cbrt (x);
#endif
}

fptype my_pow (fptype x, fptype y)
{
#ifdef SUFFIXL
  return powl (x, y);
#else
  return pow (x, y);
#endif
}

fptype my_hypot (fptype x, fptype y)
{
#ifdef SUFFIXL
  return hypotl (x, y);
#else
  return hypot (x, y);
#endif
}

fptype my_add (fptype x, fptype y)
{
  return x + y;
}

fptype my_sub (fptype x, fptype y)
{
  return x - y;
}

fptype my_mul (fptype x, fptype y)
{
  return x * y;
}

fptype my_div (fptype x, fptype y)
{
  return x / y;
}

void
test (char *foo, mp_exp_t e, unsigned long N, unsigned long seed)
{
   unsigned long i, wrong, wrong_range, wrong_monoton, wrong_symm, tot;
   mpfr_t x, y, z;
   fptype xd, xmax, xmax_dir, yd, r, range_min, range_max;
   fptype xdplus, xdminus, rplus, rminus, xdopp, ropp;
   double u, umax, umax_dir, max_err_near, max_err_dir;
   int monoton, symm;
   gmp_randstate_t state;
   mp_rnd_t rnd;

   if (verbose >= 2)
     printf ("Testing function %s for exponent %ld.\n", foo, e);
  
  if (strcmp (foo, "exp") == 0)
    {
      testfun_libm = my_exp;
      testfun_mpfr = mpfr_exp;
      range_min = 0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "exp2") == 0)
    {
      testfun_libm = my_exp2;
      testfun_mpfr = mpfr_exp2;
      range_min = 0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "expm1") == 0)
    {
      testfun_libm = my_expm1;
      testfun_mpfr = mpfr_expm1;
      range_min = -1.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    } 
  else if (strcmp (foo, "log") == 0)
    {
      testfun_libm = my_log;
      testfun_mpfr = mpfr_log;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "log2") == 0)
    {
      testfun_libm = my_log2;
      testfun_mpfr = mpfr_log2;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "log10") == 0)
    {
      testfun_libm = my_log10;
      testfun_mpfr = mpfr_log10;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "log1p") == 0)
    {
      testfun_libm = my_log1p;
      testfun_mpfr = mpfr_log1p;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "sin") == 0)
    {
      testfun_libm = my_sin;
      testfun_mpfr = mpfr_sin;
      range_min = -1.0;
      range_max = 1.0;
      monoton = NO_MONOTON;
      symm = ODD;
    }
  else if (strcmp (foo, "cos") == 0)
    {
      testfun_libm = my_cos;
      testfun_mpfr = mpfr_cos;
      range_min = -1.0;
      range_max = 1.0;
      monoton = NO_MONOTON;
      symm = EVEN;
    }
  else if (strcmp (foo, "tan") == 0)
    {
      testfun_libm = my_tan;
      testfun_mpfr = mpfr_tan;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = NO_MONOTON;
      symm = ODD;
    }
  else if (strcmp (foo, "asin") == 0)
    {
      testfun_libm = my_asin;
      testfun_mpfr = mpfr_asin;
      range_min = -1.570796326794897; /* XXX mettre le flottant = RNDD(-Pi/2) */
      range_max = 1.570796326794897; /* XXX mettre le flottant = RNDU(Pi/2) */
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "acos") == 0)
    {
      testfun_libm = my_acos;
      testfun_mpfr = mpfr_acos;
      range_min = 0.0; 
      range_max = 3.14159265358980; /* XXX mettre le flottant = RNDU(Pi) */
      monoton = DECREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "atan") == 0)
    {
      testfun_libm = my_atan;
      testfun_mpfr = mpfr_atan;
      range_min = -1.570796326794897; /* XXX mettre le flottant = RNDD(-Pi/2) */
      range_max = 1.570796326794897; /* XXX mettre le flottant = RNDU(Pi/2) */
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "sinh") == 0)
    {
      testfun_libm = my_sinh;
      testfun_mpfr = mpfr_sinh;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "cosh") == 0)
    {
      testfun_libm = my_cosh;
      testfun_mpfr = mpfr_cosh;
      range_min = 1.0;
      range_max = 1.0/0.0;
      monoton = NO_MONOTON;
      symm = EVEN;
    }
  else if (strcmp (foo, "tanh") == 0)
    {
      testfun_libm = my_tanh;
      testfun_mpfr = mpfr_tanh;
      range_min = -1.0;
      range_max = 1.0;
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "asinh") == 0)
    {
      testfun_libm = my_asinh;
      testfun_mpfr = mpfr_asinh;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "acosh") == 0)
    {
      testfun_libm = my_acosh;
      testfun_mpfr = mpfr_acosh;
      range_min = 0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "atanh") == 0)
    {
      testfun_libm = my_atanh;
      testfun_mpfr = mpfr_atanh;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = ODD;
    }
  else if (strcmp (foo, "gamma") == 0)
    {
      testfun_libm = my_tgamma;
      testfun_mpfr = mpfr_gamma;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = NO_MONOTON; /* XXX a corriger */
      symm = NO_SYMM;       /* XXX a corriger */
    }
  else if (strcmp (foo, "sqrt") == 0)
    {
      testfun_libm = my_sqrt;
      testfun_mpfr = mpfr_sqrt;
      range_min = 0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = NO_SYMM;
    }
  else if (strcmp (foo, "cbrt") == 0)
    {
      testfun_libm = my_cbrt;
      testfun_mpfr = mpfr_cbrt;
      range_min = -1.0/0.0;
      range_max = 1.0/0.0;
      monoton = INCREASING;
      symm = ODD;
    }
  else
    {
      fprintf (stderr, "Unknown function: %s\n", foo);
      exit (1);
    }

  mpfr_init2 (x, FPPREC);
  mpfr_init2 (y, FPPREC);
  mpfr_init2 (z, FPPREC);

  max_err_near = 0.0;
  max_err_dir = 0.0;

  for (rnd=0; rnd<MAX_RND; rnd++)
    {
      if (verbose >= 3)
        printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd));
      mpfr_set_machine_rnd_mode (rnd);

      /* reset the seed to test the same sequence of numbers with each
         rounding mode */
      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

      tot = 0;
      umax = 0.0; /* (signed) maximal error in ulps */
      umax_dir = 0.0; /* maximal error in ulps when wrong directed rounding */
      wrong = 0;         /* number of wrong directed roundings */
      wrong_range = 0;   /* number of wrong results wrt range        */
      wrong_monoton = 0; /* number of wrong results wrt monotonicity */
      wrong_symm = 0;    /* number of wrong results wrt symmetry     */

  for (i=0; i<N; i++)
    {
      mpfr_urandomb (x, state);
      MPFR_EXP(x) = e;
      testfun_mpfr (y, x, rnd);
      /* Conversion from a mpfr to a fp number : x */
      xd = mpfr_get_fp (x, GMP_RNDN);
      /* Conversion from a mpfr to a fp number : y */
      yd = mpfr_get_fp (y, GMP_RNDN);

      r = testfun_libm (xd);
      /* check for correct rounding */
      if (yd != r) {
        u = ulp_err (xd, r, rnd, yd);
        if (rnd != GMP_RNDN)
          {
            mp_rnd_t rnd2;
            
            if (rnd == GMP_RNDZ)
              rnd2 = (mpfr_cmp_ui (y, 0) > 0) ? GMP_RNDD : GMP_RNDU;
            else
              rnd2 = rnd;
            if ((rnd2 == GMP_RNDU && u < 0.0) || (rnd2 == GMP_RNDD && u > 0.0))
              {
                wrong++;
                if (fabs(u) > fabs(umax_dir))
                  {
                    umax_dir = u;
                    xmax_dir = xd;
                  }
              }
          }
        tot ++;
        if (fabs(u) > fabs(umax))
          {
            umax = u;
            xmax = xd;
          }
      }

      /* Testing if the range is verified */
      if (test_range && ((range_min/2.0 != range_min) || (range_max/2.0 != range_max) )) {
        /* one of the extremity is not infinite */
        if ( (r < range_min) || (r > range_max) ) {
          if (wrong_range == 0) {
              printf ("      outside range for x=");
              print_fp (xd);
              printf ("\n           f(x)=");
              print_fp (r);
              printf ("\n      not between ");
              print_fp (range_min);
              printf ("      and ");
              print_fp (range_max);
              printf (" \n");
          }
          wrong_range++;
        }
      }

      /* Testing if the monotonicity is verified */
      if (test_monotonicity && (monoton != NO_MONOTON)) {
        if (mpfr_cmp_ui(x,0) > 0) {
          mpfr_set(z, x, GMP_RNDN);
          mpfr_sub_one_ulp(z, GMP_RNDD);
          xdminus = mpfr_get_fp(z, GMP_RNDD);
          rminus = testfun_libm (xdminus);
          mpfr_set(z, x, GMP_RNDN);
          mpfr_add_one_ulp(z, GMP_RNDU);
          xdplus = mpfr_get_fp(z, GMP_RNDU);
          rplus = testfun_libm (xdplus);
        }
        else if (mpfr_cmp_ui(x,0) < 0) {
          mpfr_set(z, x, GMP_RNDN);
          mpfr_add_one_ulp(z, GMP_RNDD);
          xdminus = mpfr_get_fp(z, GMP_RNDD);
          rminus = testfun_libm (xdminus);
          mpfr_set(z, x, GMP_RNDN);
          mpfr_sub_one_ulp(z, GMP_RNDU);
          xdplus = mpfr_get_fp(z, GMP_RNDU);
          rplus = testfun_libm (xdplus);
        }
       
        if (monoton == INCREASING) {
          if (rminus > r)  {
            if (wrong_monoton == 0 && verbose >= 3) {
              printf ("      monotonicity not respected for x=");
              print_fp (xd);
              printf ("\n            f(x-)=");
              print_fp (rminus);
              printf ("\n      not <= f(x)=");
              print_fp (r);
              printf (" \n");
            }
            wrong_monoton++;
          }
          if (rplus < r)  {
            if (wrong_monoton == 0 && verbose >= 3) {
              printf ("      monotonicity not respected for x=");
              print_fp (xd);
              printf ("\n              f(x)=");
              print_fp (r);
              printf ("\n      not <= f(x+)=");
              print_fp (rplus);
              printf (" \n");
            }
            wrong_monoton++;
          }
        }
        else {/* monoton == DECREASING */
          if (rminus < r)  {
            if (wrong_monoton == 0 && verbose >= 3) {
              printf ("      monotonicity not respected for x=");
              print_fp (xd);
              printf ("\n            f(x-)=");
              print_fp (rminus);
              printf ("\n      not >= f(x)=");
              print_fp (r);
              printf (" \n");
            }
            wrong_monoton++;
          }
          if (rplus > r)  {
            if (wrong_monoton == 0 && verbose >= 3) {
              printf ("      monotonicity not respected for x=");
              print_fp (xd);
              printf ("\n              f(x)=");
              print_fp (r);
              printf ("\n      not >= f(x+)=");
              print_fp (rplus);
              printf (" \n");
            }
            wrong_monoton++;
          }
        }
      }

      /* Testing if the symmetry is verified */
      if (test_symmetry && (rnd==GMP_RNDN) && (symm != NO_SYMM) ) {
        xdopp = -xd;
        ropp = testfun_libm (xdopp);
        if (symm == ODD) {
          if (ropp != -r) {
            if (wrong_symm == 0) {
              printf ("      symmetry not respected for x=");
              print_fp (xd);
              printf ("\n           f(x)= ");
              print_fp (r);
              printf ("\n      and f(-x)=");
              print_fp (ropp);
              printf (" \n");
            }
            wrong_symm++;
          }
        }
        else { /* symm == EVEN */
          if (r != ropp) {
            if (wrong_symm == 0) {
              printf ("      symmetry not respected for x=");
              print_fp (xd);
              printf ("\n           f(x)=");
              print_fp (r);
              printf ("\n      and f(-x)=");
              print_fp (ropp);
              printf (" \n");
            }
            wrong_symm++;
          }
        }
      }

    }

  mpfr_set_machine_rnd_mode (GMP_RNDN);

  /* if any difference occured, prints the maximal ulp error */
  if (umax != 0.0 && verbose >= 3)
  {
    printf ("      %f ulp(s) for x=", umax);
    print_fp (xmax);
    printf ("\n");
    if (verbose >= 4)
      {
        printf ("      [mpfr: ");
        mpfr_set_fp (x, xmax, GMP_RNDN);
        testfun_mpfr (y, x, rnd);
        mpfr_out_str (stdout, 10, 0, y, GMP_RNDN);
        printf (" libm: ");
        print_fp (testfun_libm(xmax));
        printf ("]\n");
      }
  }

  if (wrong != 0 && verbose >= 3)
    {
      printf ("      wrong directed rounding for x=");
      print_fp (xmax_dir);
      printf (" [%f]\n", umax_dir);
    }

  umax = fabs(umax);

  if (verbose >= 3)
    {
      if (rnd == GMP_RNDN)
        printf ("   nb errors range/monotonicity/symmetry: %lu/%lu/%lu\n", 
                wrong_range, wrong_monoton, wrong_symm);
      else
        printf ("   nb errors range/monotonicity: %lu/%lu\n", 
                wrong_range, wrong_monoton);
      printf ("   nb errors/max ulp diff/wrong directed: %lu/%f/%lu\n", 
              tot, umax, wrong);
    }

  if (rnd == GMP_RNDN)
    {
      if (umax > max_err_near)
        max_err_near = umax;
    }
  else
    {
      if (umax > max_err_dir)
        max_err_dir = umax;
    }
  }

  printf ("Max. errors for %s [exp. %ld]: %f (nearest), %f (directed)\n", foo, e, max_err_near, max_err_dir);
  if (verbose >= 3)
    printf ("\n");

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);

  if (max_err_near > MAX_ERR_NEAR)
    MAX_ERR_NEAR = max_err_near;

  if (max_err_dir > MAX_ERR_DIR)
    MAX_ERR_DIR = max_err_dir;
  }

void
test2 (char *foo, mp_exp_t e, mp_exp_t f, unsigned long N, unsigned long seed)
{
  unsigned long i, wrong, tot;
  mpfr_t x, y, z, t;
  double u, umax, umax_dir, max_err_near, max_err_dir;
  fptype xd, td, r, yd, xmax, xmax_dir, tmax, tmax_dir;
  mp_rnd_t rnd;
  gmp_randstate_t state;

  if (verbose >= 2)
    printf ("Testing function %s for exponents %ld and %ld.\n", foo, e, f);

  if (strcmp (foo, "pow") == 0)
    {
      testfun_libm2 = my_pow;
      testfun_mpfr = mpfr_pow;
    }
  else if (strcmp (foo, "hypot") == 0)
    {
      testfun_libm2 = my_hypot;
      testfun_mpfr = mpfr_hypot;
    }
  else if (strcmp (foo, "add") == 0)
    {
      testfun_libm2 = my_add;
      testfun_mpfr = mpfr_add;
    }
  else if (strcmp (foo, "sub") == 0)
    {
      testfun_libm2 = my_sub;
      testfun_mpfr = mpfr_sub;
    }
  else if (strcmp (foo, "mul") == 0)
    {
      testfun_libm2 = my_mul;
      testfun_mpfr = mpfr_mul;
    }
  else if (strcmp (foo, "div") == 0)
    {
      testfun_libm2 = my_div;
      testfun_mpfr = mpfr_div;
    }
  else
    {
      fprintf (stderr, "Unknown function: %s\n", foo);
      exit (1);
    }

  mpfr_init2 (x, FPPREC);
  mpfr_init2 (y, FPPREC);
  mpfr_init2 (z, FPPREC);
  mpfr_init2 (t, FPPREC);

  max_err_near = 0.0;
  max_err_dir = 0.0;

  for (rnd=0; rnd<MAX_RND; rnd++)
    {
      if (verbose >= 3)
        printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd));
      mpfr_set_machine_rnd_mode (rnd);

      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

  umax = 0.0;
  umax_dir = 0.0;
  tot = 0;
  wrong = 0;
  for (i=0; i<N; i++)
    {
      mpfr_urandomb (x, state);
      MPFR_EXP(x) = e;
      mpfr_urandomb (t, state);
      MPFR_EXP(t) = f;
      testfun_mpfr (y, x, t, rnd);
      xd = mpfr_get_fp (x, GMP_RNDN);
      td = mpfr_get_fp (t, GMP_RNDN);
      r = testfun_libm2 (xd, td);

      /* check for correct directed rounding */
      yd = mpfr_get_fp (y, GMP_RNDN);
      if (yd != r) {
        u = ulp_err2 (xd, td, r, rnd, yd);
        if (rnd != GMP_RNDN)
          {
            mp_rnd_t rnd2;
            
            if (rnd == GMP_RNDZ)
              rnd2 = (mpfr_cmp_ui (y, 0) > 0) ? GMP_RNDD : GMP_RNDU;
            else
              rnd2 = rnd;
            if ((rnd2 == GMP_RNDU && u < 0.0) || (rnd2 == GMP_RNDD && u > 0.0))
              {
                wrong++;
                if (fabs(u) > fabs(umax_dir))
                  {
                    umax_dir = u;
                    xmax_dir = xd;
                    tmax_dir = td;
                  }
              }
          }
        tot ++;
        if (fabs(u) > fabs(umax))
          {
            umax = u;
            xmax = xd;
            tmax = td;
          }
      }
    }

  mpfr_set_machine_rnd_mode (GMP_RNDN);

  if (umax != 0.0 && verbose >= 3)
    {
      printf ("      %f ulp(s) for x=", umax);
      print_fp (xmax);
      printf (" t=");
      print_fp (tmax);
      printf ("\n");
    }

  if (umax_dir != 0.0 && verbose >= 3)
    {
      printf ("      wrong directed rounding for x=");
      print_fp (xmax_dir);
      printf (" t=");
      print_fp (tmax_dir);
      printf (" [%f]\n", umax_dir);
    }

  umax = fabs (umax);

  if (verbose >= 3)
    printf ("   nb errors/max ulp diff/wrong directed: %lu/%f/%lu\n", 
            tot, umax, wrong);

  if (rnd == GMP_RNDN)
    {
      if (umax > max_err_near)
        max_err_near = umax;
    }
  else
    {
      if (umax > max_err_dir)
        max_err_dir = umax;
    }
    }

  printf ("Max. errors for %s [exp. %ld]: %f (nearest), %f (directed)\n", foo, e, max_err_near, max_err_dir);
  if (verbose >= 3)
    printf ("\n");

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);

  if (max_err_near > MAX_ERR_NEAR)
    MAX_ERR_NEAR = max_err_near;

  if (max_err_dir > MAX_ERR_DIR)
    MAX_ERR_DIR = max_err_dir;
}

void
testall (unsigned long N, unsigned long seed)
{
  test ("exp",    0, N, seed);
  test ("exp",    9, N, seed);
  test ("exp2",   0, N, seed);
  test ("exp2",   9, N, seed);
  test ("expm1",  0, N, seed);
  test ("expm1", -9, N, seed);
  test ("log",    0, N, seed);
  test ("log", EMAX, N, seed);
  test ("log2",   0, N, seed);
  test ("log2", EMAX, N, seed);
  test ("log10",  0, N, seed);
  test ("log10", EMAX, N, seed);
  test ("log1p",  0, N, seed);
  test ("log1p", EMAX, N, seed);
  test ("sin",    0, N, seed);
  test ("sin",   10, N, seed); /* mpfr-2.0.1 is too slow for 1024 */
  test ("cos",    0, N, seed);
  test ("cos",   10, N, seed);
  test ("tan",    0, N, seed);
  test ("tan",   10, N, seed);
  test ("asin",   0, N, seed);
  test ("asin", -10, N, seed); /* mpfr-2.0.1 is too slow for -1021 */
  test ("acos",   0, N, seed);
  test ("acos", -10, N, seed); /* mpfr-2.0.1 is too slow for -1021 */
  test ("atan",   0, N, seed);
  test ("atan",  53, N, seed); /* why 53 ??? mpfr too slow ? */
  test ("sinh",   0, N, seed);
  test ("sinh",   9, N, seed);
  test ("cosh",   0, N, seed);
  test ("cosh",   9, N, seed);
  test ("tanh",   0, N, seed);
  test ("tanh",   4, N, seed);
  test ("asinh",  0, N, seed);
  test ("asinh", EMAX, N, seed);
  test ("acosh",  1, N, seed);
  test ("acosh", EMAX, N, seed);
  test ("atanh",  0, N, seed);
  test ("atanh", -10, N, seed);
  test ("gamma",  0, N, seed);
#if (FPPREC <= 53)
  test ("gamma",  7, N, seed);
#else
  test ("gamma", 10, N, seed);
#endif
  test ("sqrt",  0, N, seed);
  test ("sqrt",  EMAX, N, seed);
  test ("sqrt",  EMIN, N, seed);
  test ("cbrt",  0, N, seed);
  test ("cbrt",  EMAX, N, seed);
  test ("cbrt",  EMIN, N, seed);
  test2 ("pow", 0, 0, N, seed);
#if (FPPREC <= 53)
  test2 ("pow", 8, 7, N, seed);
#else
  test2 ("pow", 16, 10, N, seed);
#endif
  test2 ("hypot", 0, 0, N, seed);
  test2 ("hypot", EMAX-1, EMAX-1, N, seed);
  test2 ("hypot", EMIN, EMIN, N, seed);
  test2 ("add", 0, 0, N, seed);
  test2 ("add", EMAX-1, EMAX-1, N, seed);
  test2 ("sub", EMAX, EMAX, N, seed);
  test2 ("sub", 0, 0, N, seed);
  test2 ("mul", 0, 0, N, seed);
  test2 ("mul", EMAX/2, EMAX/2, N, seed);
  test2 ("div", 0, 0, N, seed);
  test2 ("div", EMAX, EMAX, N, seed);

  printf ("Maximal errors for all functions: %f (nearest), %f (directed)\n",
          MAX_ERR_NEAR, MAX_ERR_DIR);
}

int
main (int argc, char *argv[])
{
  mp_exp_t exponent, exp2 = 0;
  unsigned long N, nargs, seed = 1;

  check_fp ();

  /* sets the minimum and maximum exponents */
  mpfr_set_emax (EMAX);
  mpfr_set_emin (EMIN);

  while (argc > 1 && argv[1][0] == '-')
    {
      
      if (strcmp (argv[1], "-seed") == 0)
        {
          seed = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (strncmp (argv[1], "-verb", 4) == 0)
        {
          verbose = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (strncmp (argv[1], "-mono", 4) == 0)
        {
          test_monotonicity = 0;
          argc -= 1;
          argv += 1;
        }
      else if (strncmp (argv[1], "-rang", 4) == 0)
        {
          test_range = 0;
          argc -= 1;
          argv += 1;
        }
      else if (strncmp (argv[1], "-symm", 4) == 0)
        {
          test_symmetry = 0;
          argc -= 1;
          argv += 1;
        }
    }

  if (argc == 1 || argc == 3)
    {
      fprintf (stderr, "Usage: mpcheck [options] N\n");
      fprintf (stderr, "Usage: mpcheck [options] <function> <exponent> [N]\n");
      fprintf (stderr, "Usage: mpcheck [options] <function> <exp1> <exp2> [N]\n");
      fprintf (stderr, "where options are:\n");
      fprintf (stderr, "-seed s: set random seed to s [default 1]\n");
      fprintf (stderr, "-verb k: set verbose level to k [default 3]\n");
      fprintf (stderr, "-range : do not test output range\n");
      fprintf (stderr, "-mono  : do not test monotonicity\n");
      exit (1);
    }

  fprintf (stderr, "*******************************************************************\n");
  fprintf (stderr, "*                                                                 *\n");
  fprintf (stderr, "* MpCheck version 1.0 (c) INRIA 2002 (projects Arenaire & Spaces) *\n");
  fprintf (stderr, "*                                                                 *\n");
  fprintf (stderr, "*******************************************************************\n");

  fprintf (stderr, "[precision=%u, seed=%u]\n", FPPREC, seed);

  if (argc == 2)
    {
      N = atoi(argv[1]);
      testall (N, seed);
    }
  else
    {
      nargs = 1;
      if (strcmp (argv[1], "pow") == 0 ||
          strcmp (argv[1], "hypot") == 0 ||
          strcmp (argv[1], "add") == 0 ||
          strcmp (argv[1], "sub") == 0 ||
          strcmp (argv[1], "mul") == 0 ||
          strcmp (argv[1], "div") == 0)
        {
          nargs = 2;
          exp2  = atoi (argv[3]);
        }
  
      exponent = atoi (argv[2]);
      N = (argc > 2 + nargs) ? atoi(argv[2 + nargs]) : 1;

      if (nargs == 1)
        test (argv[1], exponent, N, seed);
      else
        test2 (argv[1], exponent, exp2, N, seed);
    }

  return 0;
}
