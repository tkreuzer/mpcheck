/*
 * MPCHECK - Check native Math functions
 * Copyright (C) 2002, 2004, 2005 INRIA
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef fptype
# error "NOT FP-TYPE selected!"
#endif

/* Include how to select a native rounding mode */
#include "rnd_mode.c"

#if defined (__i386__) && defined (__GNUC__)
# define _FPU_EXTENDED 0x0300
# define _FPU_DOUBLE   0x0200
# define _FPU_DEFAULT  0x137f
# define __setfpucw(cw) do { int toto = (cw); __asm__ ("fldcw %0" : : "m" (toto)); } while (0)
# define _fpu_ieee ((_FPU_DEFAULT & (~_FPU_EXTENDED)) | _FPU_DOUBLE)
# define set_double() __setfpucw(_fpu_ieee)
# define set_float()  __setfpucw(_FPU_DEFAULT & (~_FPU_EXTENDED))
#else
# define set_double()
# define set_float()
#endif

static mp_prec_t prec;
static mp_exp_t  emin, emax;
#ifndef RND
static int RND = ALL_RND;
#endif

static void
setup_native (void)
{
  volatile fptype x, y;
  int j;

  /* Get the precison of fptype */
  set_rnd_mode (GMP_RNDN);
  x = 1.0;
  for (j = 0; (y = x + 1.0) != x; j++, x = 2.0 * x);
  prec = j;

  /* Special checking for double numbers */
  if (prec == 53)
    {
      fptype c, d, dj;
      set_double ();
      x = DBL_MIN;
      if (2.0 * (x / 2.0) != x)
	fprintf (stderr, "WARNING: no denormalized numbers\n");
      
      c = 1.46484375e-3;
      dj = 1.0;
      for (j=0; j<54; j++) dj *= 0.5;
      d = 0.75 + dj;
      d /= 1 << 9;
      if (c != d)
	{
	  fprintf (stderr, "Default seems to use extended precision\n");
	  exit (1);
	}
    }
  else if (prec == 24)
    set_float ();

  /* Checks that the exponents EMIN and EMAX are correct */
  x = 0.5; /* EXP(x)=0 */
  for (j=0; x + x != x; j++, x = 2.0 * x);
  emax = j-1;
  for (x=1.0, j=0; j<prec; j++, x = x / 2.0);
  x = 0.5 + x;
  /* x = 0.5 + 2^(-FPPREC), EXP(x)=0 */
  for (j=0; ; j--)
    {
      y = x;
      x = x / 2.0;
      if (2.0 * x != y)
        break;
    }
  emin = j;

#ifndef RND
  /* Checks which rounding mode works */
  /* Can't use NAME(exp) seems expf or expl may not be available */
  RND = 0;
  x = 5.0;
  if (set_rnd_mode (GMP_RNDN))
    {
      y = exp (x);
      if (y > 140 && y < 160)
	RND |= 1;
      else
	fprintf (stderr, "WARNING: Rounding To Nearest doesn't work.\n");
    }
  if (set_rnd_mode (GMP_RNDZ))
    {
      y = exp (x);
      if (y > 140 && y < 160)
        RND |= 2;
      else
        fprintf (stderr, "WARNING: Rounding To Zero doesn't work.\n");
    }
  if (set_rnd_mode (GMP_RNDU))
    {
      y = exp (x);
      if (y > 140 && y < 160)
        RND |= 4;
      else
        fprintf (stderr, "WARNING: Rounding Up doesn't work.\n");
    }
  if (set_rnd_mode (GMP_RNDD))
    {
      y = exp (x);
      if (y > 140 && y < 160)
        RND |= 8;
      else
        fprintf (stderr, "WARNING: Rounding Down doesn't work.\n");
    }
#endif

}

static void *new_fp (mp_prec_t p)
{
  void *fp;
  if (p != prec)
    {
      fprintf (stderr, "ERROR: Requested prec: %lu. FPTYPE prec: %lu\n",
	       p, prec);
      abort ();
    }
  fp = malloc (sizeof (fptype));
  if (fp == NULL)
    {
      fprintf (stderr, "ERROR: Can't allocate memory");
      abort ();
    }
  return fp;
}
static void del_fp (void *fp)
{
  free (fp);
}

/* Define the wrapper for the libm */
#ifdef NAME

static void my_add (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = *(fptype*)a + *(fptype*) b;
}
static void my_sub (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = *(fptype*)a - *(fptype*) b;
}
static void my_mul (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = *(fptype*)a * *(fptype*) b;
}
static void my_div (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = *(fptype*)a / *(fptype*) b;
}

#if HAVE_SQRT
static void my_sqrt (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(sqrt) (*(fptype*)a);
}
#endif
#if HAVE_EXP
static void my_exp (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(exp) (*(fptype*)a);
}
#endif
#if HAVE_LOG
static void my_log (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(log) (*(fptype*)a);
}
#endif
#if HAVE_SIN
static void my_sin (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(sin) (*(fptype*)a);
}
#endif
#if HAVE_COS
static void my_cos (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(cos) (*(fptype*)a);
}
#endif
#if HAVE_TAN
static void my_tan (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(tan) (*(fptype*)a);
}
#endif
#if HAVE_ASIN
static void my_asin (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(asin) (*(fptype*)a);
}
#endif
#if HAVE_ACOS
static void my_acos (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(acos) (*(fptype*)a);
}
#endif
#if HAVE_ATAN
static void my_atan (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(atan) (*(fptype*)a);
}
#endif
#if HAVE_ATAN2
static void my_atan2 (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(atan2) (*(fptype*)a, *(fptype*)b);
}
#endif
#if HAVE_SINH
static void my_sinh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(sinh) (*(fptype*)a);
}
#endif
#if HAVE_COSH
static void my_cosh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(cosh) (*(fptype*)a);
}
#endif
#if HAVE_TANH
static void my_tanh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(tanh) (*(fptype*)a);
}
#endif
#if HAVE_ASINH
static void my_asinh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(asinh) (*(fptype*)a);
}
#endif
#if HAVE_ACOSH
static void my_acosh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(acosh) (*(fptype*)a);
}
#endif
#if HAVE_ATANH
static void my_atanh (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(atanh) (*(fptype*)a);
}
#endif
#if HAVE_TGAMMA
static void my_gamma (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(tgamma) (*(fptype*)a);
}
#endif
#if HAVE_EXP2
static void my_exp2 (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(exp2) (*(fptype*)a);
}
#endif
#if HAVE_LOG2
static void my_log2 (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(log2) (*(fptype*)a);
}
#endif
#if HAVE_EXPM1
static void my_expm1 (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(expm1) (*(fptype*)a);
}
#endif
#if HAVE_LOG1P
static void my_log1p (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(log1p) (*(fptype*)a);
}
#endif
#if HAVE_LOG10
static void my_log10 (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(log10) (*(fptype*)a);
}
#endif
#if HAVE_CBRT
static void my_cbrt (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(cbrt) (*(fptype*)a);
}
#endif
#if HAVE_HYPOT
static void my_hypot (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(hypot) (*(fptype*)a, *(fptype*)b);
}
#endif
#if HAVE_POW
static void my_pow (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(pow) (*(fptype*)a, *(fptype*)b);
}
#endif
#if HAVE_ERF
static void my_erf (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(erf) (*(fptype*)a);
}
#endif
#if HAVE_ERFC
static void my_erfc (void *dest, const void *a, const void *b)
{
  *(fptype*) dest = NAME(erfc) (*(fptype*)a);
}
#endif

static mpcheck_user_func_t tab[] = {
  {"add", my_add, 0, 0},
  {"add", my_add, LONG_MAX, LONG_MAX},
  {"sub", my_sub, LONG_MAX, LONG_MAX},
  {"sub", my_sub, 0, 0},
  {"mul", my_mul, 0, 0},
  {"mul", my_mul, LONG_MAX-1, LONG_MAX-1},
  {"div", my_div, 0, 0},
  {"div", my_div, LONG_MAX, LONG_MAX},
#if HAVE_SQRT
  {"sqrt", my_sqrt, 0, 0},
  {"sqrt", my_sqrt, LONG_MAX, 0},
  {"sqrt", my_sqrt, LONG_MAX-2, 0},
#endif
#if HAVE_EXP
  {"exp", my_exp, 0, 0},
  {"exp", my_exp, 9, 0}, /* FIXME: Improve by autodection of overflow? */
#endif
#if HAVE_LOG
  {"log", my_log, 0, 0},
  {"log", my_log, LONG_MAX, 0},
#endif
#if HAVE_SIN
  {"sin", my_sin, 0, 0},
  {"sin", my_sin, 10, 0},
#endif
#if HAVE_COS
  {"cos", my_cos, 0, 0},
  {"cos", my_cos, 10, 0},
#endif
#if HAVE_TAN
  {"tan", my_tan, 0, 0},
  {"tan", my_tan, 10, 0},
#endif
#if HAVE_ATAN
  {"atan", my_atan, 0, 0},
  {"atan", my_atan, 53, 0},
#endif
#if HAVE_ATAN2
  {"atan2", my_atan2, 0, 0},
  {"atan2", my_atan2, 53, 0},
#endif
#if HAVE_ASIN
  {"asin", my_asin, 0, 0},
  {"asin", my_asin, -10, 0},
#endif
#if HAVE_ACOS
  {"acos", my_acos, 0, 0},
  {"acos", my_acos, -10, 0},
#endif
#if HAVE_SINH
  {"sinh", my_sinh, 0, 0},
  {"sinh", my_sinh, 9, 0}, /* TODO: Improve overflow detection */
#endif
#if HAVE_COSH
  {"cosh", my_cosh, 0, 0},
  {"cosh", my_cosh, 9, 0}, /* TODO: Improve overflow detection */
#endif
#if HAVE_TANH
  {"tanh", my_tanh, 0, 0},
  {"tanh", my_tanh, 4, 0}, 
#endif
#if HAVE_ASINH
  {"asinh", my_asinh, 0, 0},
  {"asinh", my_asinh, 9, 0}, /* TODO */
#endif
#if HAVE_ACOSH
  {"acosh", my_acosh, 1, 0},
  {"acosh", my_acosh, 9, 0},
#endif
#if HAVE_ATANH
  {"atanh", my_atanh, 0, 0},
  {"atanh", my_atanh,-10, 0},
#endif
#if HAVE_CBRT
  {"cbrt", my_cbrt, 0, 0},
  {"cbrt", my_cbrt, LONG_MAX, 0},
  {"cbrt", my_cbrt, -1010, 0},
#endif
#if HAVE_HYPOT
  {"hypot", my_hypot, 0, 0},
  {"hypot", my_hypot, LONG_MAX, LONG_MAX},
  {"hypot", my_hypot, -1010, -1010},
#endif
#if HAVE_TGAMMA
  {"gamma", my_gamma, 0, 0},
  {"gamma", my_gamma, 10, 0},
#endif
#if HAVE_EXP2
  {"exp2", my_exp2, 0, 0},
  {"exp2", my_exp2, 9, 0},
#endif
#if HAVE_LOG2
  {"log2", my_log2,  0, 0},
  {"log2", my_log2, LONG_MAX,0},
#endif
#if HAVE_EXPM1
  {"expm1", my_expm1, 0, 0},
  {"expm1", my_expm1, -9, 0},
#endif
#if HAVE_LOG10
  {"log10", my_log10, 0, 0},
  {"log10", my_log10, LONG_MAX},
#endif
#if HAVE_LOG1P
  {"log1p", my_log1p, 0, 0},
  {"log1p", my_log1p, LONG_MAX, 0},
#endif
#if HAVE_POW
  {"pow", my_pow, 0, 0},
  {"pow", my_pow, 5, 4},
  /* {"pow", my_pow, 16, 10},  old snapshot of MPFR 2.2.0 are even too buggy*/
#endif
#if HAVE_ERF
  {"erf", my_erf, 0, 0},
  {"erf", my_erf, 9, 0},
#endif
#if HAVE_ERFC
  {"erfc", my_erfc, 0, 0},
  {"erfc", my_erfc, 2, 0},
#endif
  {NULL, NULL, 0, 0}
};
#endif /* NAME */

int main (int argc, const char *const argv[])
{
  int i;

  /* print command line */
  fprintf (stderr, "%s", argv[0]);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  setup_native ();
#ifdef LIB_INIT
  LIB_INIT ();
#endif
  mpcheck_init (argc, argv, prec, emin, emax,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		RND, ALL_TEST, 0, DEFAULT_N, 2); 
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
#ifdef LIB_EXIT
  LIB_EXIT ();
#endif
  return 0;
}
