/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

/* To adapt to the tested precision */
#define FPPREC 53       /* 24       53         64       */

#if (FPPREC == 24)
#define fptype float
#define EMIN -125
#define EMAX 128
#elif (FPPREC == 53)
#define fptype double
#define EMIN -1021
#define EMAX 1024
#elif (FPPREC == 64)
#define fptype extended
#define _FPWIDETYPES /* for HP-UX */
#define EMIN -16381
#define EMAX 16384
#elif (FPPREC == 113)
#define _FPWIDETYPES /* for HP-UX */
#define fptype quad
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

void test _PROTO ((char *, mp_exp_t, unsigned long, unsigned long));
void test2 _PROTO ((char *, mp_exp_t, mp_exp_t, unsigned long, unsigned long));
fptype ulp_err _PROTO ((fptype, fptype, mp_rnd_t));
fptype ulp_err2 _PROTO ((fptype, fptype, fptype, mp_rnd_t));
void testall _PROTO ((unsigned long, unsigned long));
void check_fp _PROTO ((void));

/* sets the machine rounding mode to the value rnd_mode */
fptype (*testfun_libm) (fptype) = NULL;
fptype (*testfun_libm2) (fptype, fptype) = NULL;
int (*testfun_mpfr) () = NULL;
double MAX_ERR_NEAR = 0.0;
double MAX_ERR_DIR = 0.0;

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
check_fp ()
{
  double x, y, c, d, dj;
  int j;

#if (FPPREC == 53) && defined(__i386)
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
}

fptype ulp_err (fptype x, fptype y, mp_rnd_t rnd)
{
  mpfr_t xx, yy;
  mp_exp_t expy;
  fptype u;
  
  mpfr_init2 (xx, 2 * FPPREC);
  mpfr_init2 (yy, 2 * FPPREC);
  /*  if (FPPREC > 53) abort (); */
  mpfr_set_d (xx, (double) x, GMP_RNDN); /* exact for FPPREC <= 53 */
  testfun_mpfr (yy, xx, rnd);
  expy = MPFR_EXP (yy);
  mpfr_set_d (xx, (double) y, GMP_RNDN); /* exact */
  mpfr_sub (xx, xx, yy, GMP_RNDN);
  MPFR_EXP(xx) += FPPREC - expy;
  u = (fptype) mpfr_get_d1 (xx);
  mpfr_clear (xx);
  mpfr_clear (yy);
  return u;
}

fptype ulp_err2 (fptype x, fptype t, fptype y, mp_rnd_t rnd)
{
  mpfr_t xx, yy, tt;
  mp_exp_t expy;
  fptype u;
  
  mpfr_init2 (xx, 2 * FPPREC);
  mpfr_init2 (tt, 2 * FPPREC);
  mpfr_init2 (yy, 2 * FPPREC);
  /*  if (FPPREC > 53) abort (); */
  mpfr_set_d (xx, (double) x, GMP_RNDN); /* exact */
  mpfr_set_d (tt, (double) t, GMP_RNDN); /* exact */
  testfun_mpfr (yy, xx, tt, rnd);
  expy = MPFR_EXP (yy);
  mpfr_set_d (xx, (double) y, GMP_RNDN); /* exact */
  mpfr_sub (xx, xx, yy, GMP_RNDN);
  MPFR_EXP(xx) += FPPREC - expy;
  u = (fptype) mpfr_get_d1 (xx);
  mpfr_clear (xx);
  mpfr_clear (tt);
  mpfr_clear (yy);
  return u;
}

fptype my_exp (fptype x)
{
  return exp (x);
}

fptype my_exp2 (fptype x)
{
  return exp2 (x);
}

fptype my_expm1 (fptype x)
{
  return expm1 (x);
}

fptype my_log (fptype x)
{
  return log (x);
}

fptype my_log10 (fptype x)
{
  return log10 (x);
}

fptype my_log1p (fptype x)
{
  return log1p (x);
}

fptype my_sin (fptype x)
{
  return sin (x);
}

fptype my_tgamma (fptype x)
{
  return tgamma (x);
}

fptype my_pow (fptype x, fptype y)
{
  return pow (x, y);
}

void
test (char *foo, mp_exp_t e, unsigned long N, unsigned long seed)
{
   unsigned long i, wrong, tot;
   mpfr_t x, y, z;
   fptype xd, yd, r;
   double u, umax, xmax, max_err_near, max_err_dir;
   gmp_randstate_t state;
   mp_rnd_t rnd;

   printf ("Testing function %s for exponent %ld.\n", foo, e);
   fflush (stdout);
  
  if (strcmp (foo, "exp") == 0)
    {
      testfun_libm = my_exp;
      testfun_mpfr = mpfr_exp;
    }
  else if (strcmp (foo, "exp2") == 0)
    {
      testfun_libm = my_exp2;
      testfun_mpfr = mpfr_exp2;
    }
  else if (strcmp (foo, "expm1") == 0)
    {
      testfun_libm = my_expm1;
      testfun_mpfr = mpfr_expm1;
    } 
  else if (strcmp (foo, "log") == 0)
    {
      testfun_libm = my_log;
      testfun_mpfr = mpfr_log;
    }
  else if (strcmp (foo, "log10") == 0)
    {
      testfun_libm = my_log10;
      testfun_mpfr = mpfr_log10;
    }
  else if (strcmp (foo, "log1p") == 0)
    {
      testfun_libm = my_log1p;
      testfun_mpfr = mpfr_log1p;
    }
  else if (strcmp (foo, "sin") == 0)
    {
      testfun_libm = my_sin;
      testfun_mpfr = mpfr_sin;
    }
  else if (strcmp (foo, "cos") == 0)
    {
      testfun_libm = (void*) cos;
      testfun_mpfr = mpfr_cos;
    }
  else if (strcmp (foo, "tan") == 0)
    {
      testfun_libm = (void*) tan;
      testfun_mpfr = mpfr_tan;
    }
  else if (strcmp (foo, "asin") == 0)
    {
      testfun_libm = (void*) asin;
      testfun_mpfr = mpfr_asin;
    }
  else if (strcmp (foo, "acos") == 0)
    {
      testfun_libm = (void*) acos;
      testfun_mpfr = mpfr_acos;
    }
  else if (strcmp (foo, "atan") == 0)
    {
      testfun_libm = (void*) atan;
      testfun_mpfr = mpfr_atan;
    }
  else if (strcmp (foo, "sinh") == 0)
    {
      testfun_libm = (void*) sinh;
      testfun_mpfr = mpfr_sinh;
    }
  else if (strcmp (foo, "cosh") == 0)
    {
      testfun_libm = (void*) cosh;
      testfun_mpfr = mpfr_cosh;
    }
  else if (strcmp (foo, "tanh") == 0)
    {
      testfun_libm = (void*) tanh;
      testfun_mpfr = mpfr_tanh;
    }
  else if (strcmp (foo, "asinh") == 0)
    {
      testfun_libm = (void*) asinh;
      testfun_mpfr = mpfr_asinh;
    }
  else if (strcmp (foo, "acosh") == 0)
    {
      testfun_libm = (void*) acosh;
      testfun_mpfr = mpfr_acosh;
    }
  else if (strcmp (foo, "atanh") == 0)
    {
      testfun_libm = (void*) atanh;
      testfun_mpfr = mpfr_atanh;
    }
#if 0
  else if (strcmp (foo, "gamma") == 0)
    {
      testfun_libm = my_tgamma;
      testfun_mpfr = mpfr_gamma;
    }
#endif
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
      printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd));
      fflush (stdout);
      mpfr_set_machine_rnd_mode (rnd);

      /* reset the seed to test the same sequence of numbers with each
         rounding mode */
      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

      tot = 0;
  umax = 0.0;
  xmax = 0.0;
  wrong = 0;
  for (i=0; i<N; i++)
    {
      mpfr_urandomb (x, state);
      MPFR_EXP(x) = e;
      testfun_mpfr (y, x, rnd);
      /* Conversion from a mpfr to a fp number : x */
      xd = (fptype) mpfr_get_d1 (x);
      /* Conversion from a mpfr to a fp number : y */
      yd = (fptype) mpfr_get_d1 (y);

      r = testfun_libm (xd);
      /* check for correct rounding */
      if (yd != r) {
        u = ulp_err (xd, r, rnd);
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
                if (wrong == 1)
                  {
                    printf ("      wrong directed rounding for x=%1.20e"
                            " [%f]\n", xd, u);
                    /* printf ("      got %1.20e instead of %1.20e\n", r, yd); */
                    fflush (stdout);
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
    }

  if (umax != 0.0)
    printf ("      %f ulp(s) for x=%1.20e\n", umax, xmax);

  umax = fabs(umax);

  printf ("   nb errors/max ulp diff/wrong directed: %lu/%f/%lu\n", 
          tot, umax, wrong);
  fflush (stdout);

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

  printf ("Maximal errors for %s: %f (nearest), %f (directed)\n", foo,
          max_err_near, max_err_dir);
  fflush (stdout);

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
  double xd, td, r, u, umax, xmax, tmax, max_err_near, max_err_dir;
  mp_rnd_t rnd;
  gmp_randstate_t state;

  printf ("Testing function %s for exponents %ld and %ld.\n", foo, e, f);
  fflush (stdout);

  if (strcmp (foo, "pow") == 0)
    {
      testfun_libm2 = (void*) pow;
      testfun_mpfr = mpfr_pow;
    }
  else
    {
      fprintf (stderr, "Unknown function: %s\n", foo);
      exit (1);
    }

  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);
  mpfr_init2 (z, 53);
  mpfr_init2 (t, 53);

  max_err_near = 0.0;
  max_err_dir = 0.0;

  for (rnd=0; rnd<MAX_RND; rnd++)
    {
      printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd));
      fflush (stdout);
      mpfr_set_machine_rnd_mode (rnd);

      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

  umax = 0.0;
  xmax = 0.0;
  tmax = 0.0;
  tot = 0;
  wrong = 0;
  for (i=0; i<N; i++)
    {
      mpfr_urandomb (x, state);
      MPFR_EXP(x) = e;
      mpfr_urandomb (t, state);
      MPFR_EXP(t) = f;
      testfun_mpfr (y, x, t, rnd);
      xd = mpfr_get_d1 (x);
      td = mpfr_get_d1 (t);
      r = testfun_libm2 (xd, td);

      /* check for correct directed rounding */
      if (mpfr_get_d1 (y) != r) {
        u = ulp_err2 (xd, td, r, rnd);
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
                if (wrong == 1)
                  {
                    printf ("      wrong directed rounding for x=%1.20e "
                            " t=%1.20e [%f]\n", xd, td, u);
                    fflush (stdout);
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

  if (umax != 0.0)
    printf ("      %f ulp(s) for x=%1.20e t=%1.20e\n", umax, xmax, tmax);

  umax = fabs (umax);

  printf ("   nb errors/max ulp diff/wrong directed: %lu/%f/%lu\n", 
          tot, umax, wrong);
  fflush (stdout);

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

  printf ("Maximal errors for %s: %f (nearest), %f (directed)\n", foo,
          max_err_near, max_err_dir);
  fflush (stdout);

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
  test ("log", 1024, N, seed);
  test ("log10",  0, N, seed);
  test ("log10", 1024, N, seed);
  test ("log1p",  0, N, seed);
  test ("log1p", 1024, N, seed);
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
  test ("atan",  53, N, seed);
  test ("sinh",   0, N, seed);
  test ("sinh",   9, N, seed);
  test ("cosh",   0, N, seed);
  test ("cosh",   9, N, seed);
  test ("tanh",   0, N, seed);
  test ("tanh",   4, N, seed);
  test ("asinh",  0, N, seed);
  test ("asinh", 1024, N, seed);
  test ("acosh",  1, N, seed);
  test ("acosh", 1024, N, seed);
  test ("atanh",  0, N, seed);
  test ("atanh", -10, N, seed);
  test ("gamma",  0, N, seed);
  test ("gamma",  7, N, seed);
  test2 ("pow", 0, 0, N, seed);
  test2 ("pow", 8, 7, N, seed);

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

  if (strcmp (argv[1], "-seed") == 0)
    {
      seed = atoi (argv[2]);
      argc -= 2;
      argv += 2;
    }

  printf ("[precision=%u, seed=%u]\n", FPPREC, seed);

  if (argc == 1 || argc == 3)
    {
      fprintf (stderr, "Usage: test [-seed s] N\n");
      fprintf (stderr, "Usage: test [-seed s] <function> <exponent> [N]\n");
      fprintf (stderr, "Usage: test [-seed s] <function> <exp1> <exp2> [N]\n");
      exit (1);
    }

  if (argc == 2)
    {
      N = atoi(argv[1]);
      testall (N, seed);
    }
  else
    {
      nargs = 1;
      if (strcmp (argv[1], "pow") == 0)
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
