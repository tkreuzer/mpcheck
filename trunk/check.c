/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2002, 2004, 2005, 2010 INRIA
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

#include "mpcheck.h"

#define MAX_REPORTS 3 /* maximal number of errors reported by category */

#define OUT 10 /* output base */

mpfr_t    mpcheck_max_err_dir, mpcheck_max_err_near;
mp_rnd_t  mpcheck_rnd_mode;
char function_to_check[128] = "all"; /* default=all */
unsigned long tot_wrong_errno = 0, tot_wrong_inexact = 0,
  tot_wrong_dir = 0, tot_wrong_range = 0, tot_wrong_monoton = 0,
  tot_wrong_symm = 0, tot_wrong = 0, tot_wrong_special = 0,
  tot_wrong_overflow = 0, tot_wrong_underflow = 0,
  tot_wrong_basic = 0;
unsigned long suppressed_wrong_inexact = 0, suppressed_wrong_errno = 0,
  suppressed_wrong_underflow = 0, suppressed_wrong_range = 0,
  suppressed_wrong_overflow = 0;

/* output x in decimal and hexadecimal */
static void
out_value (FILE *fp, mpfr_t x)
{
  mpfr_out_str (fp, 10, 0, x, GMP_RNDN);
  fprintf (fp, " [");
  mpfr_out_str (fp, 16, 0, x, GMP_RNDN);
  fprintf (fp, "]");
}

/* Since we can't store the range directly (we need to compute
   it for any prec), we compute it here */
static void
mpcheck_set_range (mpfr_t dest, mpcheck_range_e range)
{
  switch ((int) range)
    {
    case RANGE_INF:
      mpfr_set_inf (dest, 1);
      break;
    case RANGE_ZERO:
      mpfr_set_ui (dest, 0, GMP_RNDU);
      break;
    case RANGE_ONE:
      mpfr_set_ui (dest, 1, GMP_RNDU);
      break;
    case RANGE_PI:
      mpfr_const_pi (dest, GMP_RNDU);
      break;
    case RANGE_PI2:
      mpfr_const_pi (dest, GMP_RNDU);
      mpfr_div_2ui (dest, dest, 1, GMP_RNDU);
      break;
    case RANGE_TWO:
      mpfr_set_ui (dest, 2, GMP_RNDU);
      break;
    case -RANGE_INF:
      mpfr_set_inf (dest, -1);
      break;
    case -RANGE_ZERO:
      mpfr_set_ui (dest, 0, GMP_RNDD);
      break;
    case -RANGE_ONE:
      mpfr_set_si (dest, -1, GMP_RNDD);
      break;
    case -RANGE_PI:
      mpfr_const_pi (dest, GMP_RNDU);
      mpfr_neg (dest, dest, GMP_RNDD);
      break;
    case -RANGE_PI2:
      mpfr_const_pi (dest, GMP_RNDU);
      mpfr_div_2ui (dest, dest, 1, GMP_RNDU);
      mpfr_neg (dest, dest, GMP_RNDD);
      break;
    case -RANGE_TWO:
      mpfr_set_si (dest, -2, GMP_RNDU);
      break;
    default:
      fprintf (stderr, "MPCHECK ERROR: Range undefined\n");
      abort ();
    }
}

/* library is the result computed by the tested library (precision prec)
   reference is the result computed by MPFR using more precision */
static void
mpcheck_ulp (mpfr_t ulp, mpfr_t library, mpfr_t reference, mp_prec_t prec)
{
  mpfr_exp_t emin, eulp;

  /* enlarge the exponent range to avoid getting zero */
  emin = mpfr_get_emin ();
  mpfr_set_emin (mpfr_get_emin_min ());
  mpfr_sub (ulp, library, reference, GMP_RNDN);
  if (mpfr_cmp_ui (ulp, 0) == 0)
    { mpfr_dump (reference); mpfr_dump (library); }
  assert (mpfr_cmp_ui (ulp, 0) != 0);
  /* ulp(reference) = 2^(EXP(reference) - prec), unless it is smaller than
     1/2*2^emin = 2^(emin-1) */
  eulp = mpfr_get_exp (reference) - prec;
  if (eulp < emin -1)
    eulp = emin - 1;
  mpfr_div_2si (ulp, ulp, eulp, GMP_RNDN);
  mpfr_set_emin (emin);
}

static void
usage (void)
{
  fprintf (stderr,
	   "Usage: mpcheck [options]\n"
	   "where options are:\n"
	   "--seed=s    : set random seed to s [default 1]\n"
	   "--verbose=k : set verbose level to k [default 2]\n"
	   "--range=bool: test output range\n"
	   "--mono=bool : test monotonicity\n"
	   "--symm=bool : test symmetry\n"
	   "--RNDN=bool : test Round To Nearest\n"
	   "--RNDZ=bool : test Round To Zero\n"
	   "--RNDD=bool : test Round To -Inf\n"
	   "--RNDU=bool : test Round To +Inf\n"
	   "--prec=val  : Set prec\n"
	   "--num=n     : N [default %u]\n", DEFAULT_N);
  exit (1);
}

/* Global variables to avoid passing too much args to functions */
static mp_prec_t prec;
static mp_exp_t emin, emax;
static void *(*new)(mp_prec_t);
static void (*del)(void *);
static void (*getfp)(void *, mpfr_srcptr);
static void (*setfp)(mpfr_ptr, const void *);
static int (*setrnd) (mp_rnd_t);
static unsigned long N;
static unsigned long seed;
static mpcheck_test_e test;
static int rnd_mode;
/* verbose levels:
   1: only print one summary line for each function + exponent range
   2 (default): print two lines per function + exp. (Testing ... + Max. errors)
   3: for each (function,exponent), and for each rounding mode,
      print the maximal error in ulps, the corresponding input value(s),
      and a summary "nb errors/wrong directed/max ulp diff"
   4: print also the results return by mpfr and the other library
      (only for the worst cases)
   5: print ERROR: ... for each difference between mpfr and the other library
*/
static int verbose;
static void (*fp_feclearexcept)(void);
static int (*fp_fetestexcept)(int flag);
    
void mpcheck_set_exception_functions(void (*fp_feclearexcept1)(void),
                                     int (*fp_fetestexcept1)(int flag))
{
    fp_feclearexcept = fp_feclearexcept1;
    fp_fetestexcept = fp_fetestexcept1;
}

void
mpcheck_init (int argc, const char *const argv[],
	      mp_prec_t prec2, mp_exp_t emin2, mp_exp_t emax2,
	      void *(*new2)(mp_prec_t), void (*del2)(void *),
	      void (*getfp2)(void *, mpfr_srcptr),
	      void (*setfp2)(mpfr_ptr, const void *),
	      int (*setrnd2)(mp_rnd_t),
	      int rnd_mode2, mpcheck_test_e test2, unsigned long seed2,
	      unsigned long N2, int verbose2)
{
  char Buffer[100];
  char *val;
  long value;
  int i;

  /* Setup default values */
  new = new2;
  del = del2;
  getfp = getfp2;
  setfp = setfp2;
  rnd_mode = rnd_mode2;
  test = test2;
  seed = seed2;
  prec = prec2;
  emin = emin2;
  emax = emax2;
  setrnd = setrnd2;
  N = N2;
  verbose = verbose2;

  /* Read input of programs and update values if needed */
  /* Options begin with '--' */
  for (i = 1 ; i < argc ; i++)
    {
      if (argv[i][0] == '-' && argv[i][1] == '-')
	{
	  /* Options are --NAME or --NAME=VAL */
	  strncpy (Buffer, argv[i]+2, 99); Buffer[99] = 0;
	  val = strchr (Buffer, '=');
	  if (val != NULL)
	    *val++ = 0;
	  value = (val == NULL || *val == 0) ? 0 : atol (val);
	  /* Deal with them */
	  if (strcmp (Buffer, "mono") == 0)
	    if (value == 0)
	      test &= ~MPCHECK_TEST_MONOTON;
	    else
	      test |= MPCHECK_TEST_MONOTON;
	  else if (strcmp (Buffer, "range") == 0)
	    if (value == 0)
	      test &= ~MPCHECK_TEST_RANGE;
	    else
	      test |= MPCHECK_TEST_RANGE;
	  else if (strcmp (Buffer, "symm") == 0)
	    if (value == 0)
	      test &= ~MPCHECK_TEST_SYMM;
	    else
	      test |= MPCHECK_TEST_SYMM;
	  else if (strcmp (Buffer, "RNDN") == 0)
	    if (value == 0)
	      rnd_mode &= ~1;
	    else
	      rnd_mode |= 1;
	  else if (strcmp (Buffer, "RNDZ") == 0)
	    if (value == 0)
	      rnd_mode &= ~2;
	    else
	      rnd_mode |= 2;
	  else if (strcmp (Buffer, "RNDU") == 0)
	    if (value == 0)
	      rnd_mode &= ~4;
	    else
	      rnd_mode |= 4;
	  else if (strcmp (Buffer, "RNDD") == 0)
	    if (value == 0)
	      rnd_mode &= ~8;
	    else
	      rnd_mode |= 8;
	  else if (strcmp (Buffer, "seed") == 0)
	    seed = value;
	  else if (strcmp (Buffer, "verbose") == 0)
	    verbose = value;
	  else if (strcmp (Buffer, "num") == 0)
	    N = value;
	  else if (strcmp (Buffer, "prec") == 0)
	    prec = value;
	  else if (strcmp (Buffer, "func") == 0)
	    strcpy (function_to_check, val);
	  else
	    {
	      usage ();
	      exit(0);
	    }
	}
      else
	{
	  usage ();
	  exit(0);
	}
    }

  /* Init */
  mpfr_set_emin (emin - prec + 1);
  mpfr_set_emax (emax);
  mpfr_set_default_prec (prec);
  mpfr_inits (mpcheck_max_err_dir, mpcheck_max_err_near, NULL);
  mpfr_set_ui (mpcheck_max_err_dir, 0, GMP_RNDN);
  mpfr_set_ui (mpcheck_max_err_near, 0, GMP_RNDN);
}

void
mpcheck_clear (FILE *out)
{
  fprintf (out, "Max. errors : ");
  mpfr_out_str (out, 10, 3,  mpcheck_max_err_near, MPFR_RNDA);
  fprintf (out, " (nearest), ");
  mpfr_out_str (out, 10, 3, mpcheck_max_err_dir, MPFR_RNDA);
  fprintf (out, " (directed) [seed=%lu]\n", seed);
  if (tot_wrong > 0)
    fprintf (out, "Incorrect roundings: %lu (basic %lu)\n",
             tot_wrong, tot_wrong_basic);
  if (tot_wrong_dir > 0)
    fprintf (out, "Wrong side of directed rounding: %lu\n", tot_wrong_dir);
  if (tot_wrong_range > 0)
    fprintf (out, "Wrong range: %lu (suppressed %lu)\n", tot_wrong_range,
             suppressed_wrong_range);
  if (tot_wrong_monoton > 0)
    fprintf (out, "Wrong monotonicity: %lu\n", tot_wrong_monoton);
  if (tot_wrong_symm > 0)
    fprintf (out, "Wrong symmetry: %lu\n", tot_wrong_symm);
  if (tot_wrong_errno > 0)
    fprintf (out, "Wrong errno: %lu (suppressed %lu)\n",
             tot_wrong_errno, suppressed_wrong_errno);
  if (tot_wrong_inexact > 0)
    fprintf (out, "Wrong inexact: %lu (suppressed %lu)\n",
             tot_wrong_inexact, suppressed_wrong_inexact);
  if (tot_wrong_special > 0)
    fprintf (out, "Wrong special values: %lu\n",
             tot_wrong_special);
  if (tot_wrong_overflow > 0)
    fprintf (out, "Wrong overflow: %lu (suppressed %lu)\n",
             tot_wrong_overflow, suppressed_wrong_overflow);
  if (tot_wrong_underflow > 0)
    fprintf (out, "Wrong underflow: %lu (suppressed %lu)\n",
             tot_wrong_underflow, suppressed_wrong_underflow);

  mpfr_clears (mpcheck_max_err_dir, mpcheck_max_err_near, NULL);
}

static void
set_special (mpfr_t op, int i)
{
  switch (i)
    {
    case 0:
      mpfr_set_zero (op, 1);
      return;
    case 1:
      mpfr_set_zero (op, -1);
      return;
    case 2:
      mpfr_set_inf (op, 1);
      return;
    case 3:
      mpfr_set_inf (op, -1);
      return;
    case 4:
      mpfr_set_nan (op);
      return;
    default:
      abort ();
    }
}

/* Put in x the largest number such that foo(x) < v for a monotone
   function foo. Return zero if no overflow is found.
   nb = 0: v is MAX_FLOAT
   nb = 1: v is -MAX_FLOAT
   nb = 2: v is MIN_FLOAT (subnormal)
   nb = 3: v is -MIN_FLOAT (subnormal)
   nb = 4: v is the smallest normal number
   nb = 5: v is minus the smallest normal number
*/
static int
find_overflow (mpfr_t x, mpcheck_func_t *ref, int nb)
{
  mpfr_t y, xmax, z, t, u, v;
  mpfr_prec_t prec = mpfr_get_prec (x);
  int ret = 0;
  mp_rnd_t rnd;

  assert (ref->NumArg == 1);
  assert (ref->monoton == INCREASING || ref->monoton == DECREASING);

  rnd = (ref->monoton == INCREASING) ? MPFR_RNDU : MPFR_RNDD;

  mpfr_init2 (y, prec);
  mpfr_init2 (z, prec);
  mpfr_init2 (t, prec);
  mpfr_init2 (u, prec);
  mpfr_init2 (v, prec);
  mpfr_init2 (xmax, prec);

  /* set xmax to (1-2^(-prec))*2^EMAX */
  mpfr_set_ui (xmax, 1, MPFR_RNDN);
  mpfr_nextbelow (xmax);
  mpfr_mul_2exp (xmax, xmax, mpfr_get_emax (), MPFR_RNDN);
  
  /* first find the smallest x such that foo(x) is not NaN */

  /* initialize x to -xmax */
  mpfr_neg (x, xmax, MPFR_RNDN);
  (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (t, x, rnd);

  if (mpfr_nan_p (t))
    {
      /* initialize y to 0 or 1, we assume foo(0) or foo(1) is not NaN */
      mpfr_set_ui (y, 0, MPFR_RNDN);
      (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (z, y, rnd);
      if (mpfr_nan_p (z)) /* needed for acosh */
        {
          mpfr_set_ui (y, 1, MPFR_RNDN);
          (*((int(*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (z, y, rnd);
          assert (mpfr_nan_p (z) == 0);
        }

      while (1) /* invariant: foo(x)=NaN, foo(y)<>NaN */
        {
          mpfr_div_2exp (z, x, 1, MPFR_RNDN);
          mpfr_div_2exp (t, y, 1, MPFR_RNDN);
          mpfr_add (z, z, t, MPFR_RNDN);
          if (mpfr_cmp (z, x) == 0 || mpfr_cmp (z, y) == 0)
            {
              mpfr_set (x, y, MPFR_RNDN);
              break;
            }
          (*((int(*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (t, z, rnd);
          if (mpfr_nan_p (t))
            mpfr_set (x, z, MPFR_RNDN);
          else
            mpfr_set (y, z, MPFR_RNDN);
        }
    }

  /* same on the right side */
  mpfr_set (y, xmax, MPFR_RNDN);
  (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (t, y, rnd);

  if (mpfr_nan_p (t))
    {
      mpfr_set_ui (u, 0, MPFR_RNDN); /* we assume foo(0) is not NaN */
      while (1) /* invariant: foo(u)<>NaN, foo(y)=NaN */
        {
          mpfr_div_2exp (z, u, 1, MPFR_RNDN);
          mpfr_div_2exp (t, y, 1, MPFR_RNDN);
          mpfr_add (z, z, t, MPFR_RNDN);
          if (mpfr_cmp (z, u) == 0 || mpfr_cmp (z, y) == 0)
            {
              mpfr_set (y, u, MPFR_RNDN);
              break;
            }
          (*((int(*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (t, z, rnd);
          if (mpfr_nan_p (t))
            mpfr_set (y, z, MPFR_RNDN);
          else
            mpfr_set (u, z, MPFR_RNDN);
        }
    }

  switch (nb)
    {
    case 0:
      mpfr_set (v, xmax, MPFR_RNDN);
      break;
    case 1:
      mpfr_neg (v, xmax, MPFR_RNDN);
      break;
    case 2:
      mpfr_set_ui (v, 0, MPFR_RNDN);
      mpfr_nextabove (v);
      mpfr_subnormalize (v, 0, MPFR_RNDU);
      break;
    case 3:
      mpfr_set_ui (v, 0, MPFR_RNDN);
      mpfr_nextbelow (v);
      mpfr_subnormalize (v, 0, MPFR_RNDD);
      break;
    case 4:
      mpfr_set_ui (v, 0, MPFR_RNDN);
      mpfr_nextabove (v);
      mpfr_subnormalize (v, 0, MPFR_RNDU);
      /* if x is the smallest subnormal number and y the smallest norma
         number, then y = 2^(prec-1) * x */
      mpfr_mul_2exp (v, v, mpfr_get_prec (v) - 1, MPFR_RNDU);
      break;
    case 5:
      mpfr_set_ui (v, 0, MPFR_RNDN);
      mpfr_nextbelow (v);
      mpfr_subnormalize (v, 0, MPFR_RNDD);
      mpfr_mul_2exp (v, v, mpfr_get_prec (v) - 1, MPFR_RNDD);
      break;
    default:
      abort();
    }

  /* check foo(x) <= v if foo is increasing, and foo(x) >= v if decreasing */
  (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (z, x, rnd);
  if ((ref->monoton == INCREASING && mpfr_cmp (z, v) > 0) ||
      (ref->monoton == DECREASING && mpfr_cmp (z, v) < 0))
    goto failure;

  /* check v <= foo(y) if foo is increasing, and v >= foo(y) if decreasing */
  (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr)) (z, y, rnd);
  if ((ref->monoton == INCREASING && mpfr_cmp (z, v) < 0) ||
      (ref->monoton == DECREASING && mpfr_cmp (z, v) > 0))
    goto failure;

  while (mpfr_cmp (x, y) < 0)
    {
      mpfr_div_2exp (z, x, 1, MPFR_RNDN);
      mpfr_div_2exp (t, y, 1, MPFR_RNDN);
      mpfr_add (z, z, t, MPFR_RNDN);
      if (mpfr_cmp (z, x) == 0 || mpfr_cmp (z, y) == 0)
        break;
      (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
         ref->mpfr)) (t, z, rnd);
      if ((ref->monoton == INCREASING && mpfr_cmp (t, v) <= 0) ||
          (ref->monoton == DECREASING && mpfr_cmp (t, v) >= 0))
        mpfr_set (x, z, MPFR_RNDN);
      else
        mpfr_set (y, z, MPFR_RNDN);
    }
  ret = 1;

 failure:
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (v);
  mpfr_clear (xmax);
  return ret;
}

/* replace op1 by op1 + j * ulp(op1) where j is a small random number */
static void
next (mpfr_t op1)
{
  long j = (lrand48 () % 64) - 32;

  while (j > 0)
    {
      mpfr_nextabove (op1);
      mpfr_subnormalize (op1, 0, MPFR_RNDU);
      j --;
    }
  while (j < 0)
    {
      mpfr_nextbelow (op1);
      mpfr_subnormalize (op1, 0, MPFR_RNDD);
      j ++;
    }
}

static void
gen (mpfr_t op1, mp_exp_t e1, mpfr_t op2, mp_exp_t e2, mpfr_t op3, mp_exp_t e3,
     gmp_randstate_t state, int i, mpcheck_func_t *ref)
{
  int Signed = ref->signed_input == IN_POSNEG;
  int numarg = ref->NumArg;

// #define NO_SPECIAL /* useful to test say Pari/GP */
#ifndef NO_SPECIAL
  if (numarg == 1 && i < 5)
    {
      set_special (op1, i);
      return;
    }
  else if (numarg == 2 && i < 25)
    {
      set_special (op1, i % 5);
      set_special (op2, i / 5);
      return;
    }
  else if (numarg == 3 && i < 125)
    {
      set_special (op1, i % 5);
      set_special (op2, (i / 5) % 5);
      set_special (op3, i / 25);
      return;
    }
  if (numarg == 1 && (5 <= i && i <= 6) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 0)) /* around +Inf */
        {
          if (i == 6)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (7 <= i && i <= 8) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 1)) /* around -Inf */
        {
          if (i == 8)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (9 <= i && i <= 10) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 2)) /* around MIN_FLOAT */
        {
          mpfr_subnormalize (op1, 0, MPFR_RNDN);
          if (i == 10)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (11 <= i && i <= 12) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 3)) /* around -MIN_FLOAT */
        {
          mpfr_subnormalize (op1, 0, MPFR_RNDN);
          if (i == 12)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (13 <= i && i <= 14) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 4)) /* around smallest normal */
        {
          mpfr_subnormalize (op1, 0, MPFR_RNDN);
          if (i == 14)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (15 <= i && i <= 16) &&
      (ref->monoton == INCREASING || ref->monoton == DECREASING))

    {
      if (find_overflow (op1, ref, 5)) /* around -smallest normal */
        {
          mpfr_subnormalize (op1, 0, MPFR_RNDN);
          if (i == 16)
            next (op1);
          return;
        }
    }
  if (numarg == 1 && (17 <= i && i <= 18))
    {
      mpfr_set_inf (op1, 1);
      mpfr_nextbelow (op1);
      if (i == 18)
        next (op1);
      return;
    }
  if (numarg == 1 && (19 <= i && i <= 20))
    {
      mpfr_set_inf (op1, -1);
      mpfr_nextabove (op1);
      if (i == 20)
        next (op1);
      return;
    }
  if (numarg == 1 && (21 <= i && i <= 22))
    {
      mpfr_set_zero (op1, 1);
      mpfr_nextabove (op1);
      if (i == 22)
        next (op1);
      return;
    }
  if (numarg == 1 && (23 <= i && i <= 24))
    {
      mpfr_set_zero (op1, -1);
      mpfr_nextbelow (op1);
      if (i == 24)
        next (op1);
      return;
    }
#endif

  do mpfr_urandomb (op1, state); while (mpfr_zero_p (op1));
  mpfr_set_exp (op1, 0);
  mpfr_mul_2si (op1, op1, e1, GMP_RNDN);
  if (numarg >= 2)
    {
      do mpfr_urandomb (op2, state); while (mpfr_zero_p (op2));
      mpfr_set_exp (op2, 0);
      mpfr_mul_2si (op2, op2, e2, GMP_RNDN);
    }
  if (numarg >= 3)
    {
      do mpfr_urandomb (op3, state); while (mpfr_zero_p (op3));
      mpfr_set_exp (op3, 0);
      mpfr_mul_2si (op3, op3, e3, GMP_RNDN);
    }
  if (Signed)
    {
      if ((rand () & 1) == 0)
        mpfr_neg (op1, op1, GMP_RNDN);
      if (numarg >= 2)
        if ((rand () & 1) == 0)
          mpfr_neg (op2, op2, GMP_RNDN);
      if (numarg >= 3)
        if ((rand () & 1) == 0)
          mpfr_neg (op3, op3, GMP_RNDN);
    }
  /* Functions which take only positive arguments
     may have a strange behaviour for 0 */
  else while (mpfr_zero_p (op1)) {
      mpfr_urandomb (op1, state);
      mpfr_mul_2si (op1, op1, e1, GMP_RNDN);
    }
}

static int
is_largest (mpfr_t x)
{
  if (mpfr_regular_p (x) && mpfr_get_exp (x) == mpfr_get_emax ())
    {
      mpfr_t y;
      int ret;

      mpfr_init2 (y, mpfr_get_prec (x));
      mpfr_set (y, x, MPFR_RNDN);
      if (mpfr_sgn (x) > 0)
        mpfr_nextabove (y);
      else
        mpfr_nextbelow (y);
      ret = mpfr_inf_p (y);
      mpfr_clear (y);
      return ret;
    }
  return 0;
}

static int
is_smallest_subnormal (mpfr_t x)
{
  if (mpfr_regular_p (x))
    {
      int sign = mpfr_signbit (x) ? -1 : 1;
      return mpfr_cmp_si_2exp (x, sign, mpfr_get_emin () - 1) == 0;
    }
  else
    return 0;
}

/* cases where we don't check the inexact flag */
static int
suppress (int err, const char *name, mpfr_t result, mpfr_t op1, mpfr_t op2)
{
#ifdef HAVE_GLIBC
  if (err == 0) /* inexact exception */
    {
      /* cbrt does not guarantee a correct inexact flag, except for NaN */
      if (strcmp (name, "cbrt") == 0 && !mpfr_nan_p (result))
        return 1;
      /* same for gamma, exp10, log2, exp2 */
      if ((strcmp (name, "gamma") == 0 || strcmp (name, "exp10") == 0 ||
           strcmp (name, "log2") == 0 || strcmp (name, "exp2") == 0) &&
          !mpfr_nan_p (result))
        return 1;
      /* glibc does not try for accurate "inexact" for NaN results, unless
         coming from NaN arguments */
      if (mpfr_nan_p (result) && !mpfr_nan_p (op1))
        return 1;
    }
  else if (err == 1) /* errno */
    {
      /* errno is only supposed to be set for <math.h> standard library
         functions, not C operators (or <complex.h> functions) */
      if (strcmp (name, "add") == 0 || strcmp (name, "sub") == 0 ||
          strcmp (name, "mul") == 0 || strcmp (name, "div") == 0)
        return 1;
      if (strcmp (name, "log1p") == 0)
        return 1; /* bug 6792 */
      if (strcmp (name, "fma") == 0)
        return 1; /* bug 6801 */
      /* man sincos says it does not set errno for x=Inf, cf bug 15467 */
      if ((strcmp (name, "sincos1") == 0 || strcmp (name, "sincos2") == 0)
          && (mpfr_inf_p (op1) || mpfr_inf_p (op2)))
        return 1;
    }
  else if (err == 2) /* underflow flags don't agree (with same result) */
    {
      /* if the result is the smallest subnormal, it is not a real bug */
      if (is_smallest_subnormal (result))
        return 1;
      if (strcmp (name, "expm1") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16539 (fixed) and 16353 */
      if (strcmp (name, "log1p") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16339 */
#if 0 /* fixed 27 Mar 2014 */
      if (strcmp (name, "exp") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16348 */
#endif
#if 0 /* fixed 23 Jun 2014 */
      if (strcmp (name, "cosh") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16354 */
#endif
#if 0 /* fixed 14 May 2014 */
      if (strcmp (name, "erf") == 0 && (mpfr_get_prec (result) == 24 ||
                                        mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16516 */
#endif
      if (strcmp (name, "sincos2") == 0)
        return 1; /* underflow might be in sin, not in cos */
      if (strcmp (name, "exp2") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 16560 */
#if 0 /* fixed 25 Jun 2014 */
      if (strcmp (name, "exp10") == 0 && (mpfr_get_prec (result) == 53 ||
                                          mpfr_get_prec (result) == 64))
        return 1; /* bug 16560 */
#endif
    }
  else if (err == 3) /* inconsistent underflow/result */
    {
      if (strcmp (name, "y0") == 0 && mpfr_get_prec (result) == 24)
        return 1; /* bug 14471 */
      if (strcmp (name, "atan") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 15319 */
      if (strcmp (name, "log1p") == 0 && (mpfr_get_prec (result) == 24 ||
                                          mpfr_get_prec (result) == 53))
        return 1; /* bug 16339 */
      if (strcmp (name, "asinh") == 0 && (mpfr_get_prec (result) == 24 ||
                                          mpfr_get_prec (result) == 53 ||
                                          mpfr_get_prec (result) == 64))
        return 1; /* bug 16350 */
      if (strcmp (name, "asin") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 16351 */
      if (strcmp (name, "atanh") == 0 && (mpfr_get_prec (result) == 24 ||
                                          mpfr_get_prec (result) == 53 ||
                                          mpfr_get_prec (result) == 64))
        return 1; /* bug 16352 */
      if (strcmp (name, "expm1") == 0 && (mpfr_get_prec (result) == 24 ||
                                          mpfr_get_prec (result) == 53 ||
                                          mpfr_get_prec (result) == 64))
        return 1; /* bug 16353 */
#if 0 /* fixed 23 Jun 2014 */
      if (strcmp (name, "cosh") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16354 */
#endif
#if 0 /* fixed 14 May 2014 */
      if (strcmp (name, "erf") == 0 && (mpfr_get_prec (result) == 24 ||
                                        mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16516 */
#endif
      if (strcmp (name, "exp") == 0 && (mpfr_get_prec (result) == 64))
        return 1; /* bug 16361 */
      if (strcmp (name, "exp10") == 0 && (mpfr_get_prec (result) == 64))
        return 1; /* bug 16361 */
      if (strcmp (name, "tan") == 0 && (mpfr_get_prec (result) == 24 ||
                                        mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16517 */
      if (strcmp (name, "sinh") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 16519 */
      if (strcmp (name, "tanh") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 16520 */
      if (strcmp (name, "exp2") == 0 && (mpfr_get_prec (result) == 24 ||
                                         mpfr_get_prec (result) == 53 ||
                                         mpfr_get_prec (result) == 64))
        return 1; /* bug 16521 */
      if (strcmp (name, "sincos1") == 0 && (mpfr_get_prec (result) == 53 ||
                                            mpfr_get_prec (result) == 64))
        return 1; /* bug 16526 */
      if (strcmp (name, "sincos2") == 0) /* underflow might be in sin! */
        return 1;
      if (strcmp (name, "sin") == 0 && (mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16538 */
      if (strcmp (name, "j42") == 0 && (mpfr_get_prec (result) == 24 ||
                                        mpfr_get_prec (result) == 53))
        return 1; /* bug 16559 */
      if (strcmp (name, "j1") == 0 && (mpfr_get_prec (result) == 24 ||
                                       mpfr_get_prec (result) == 53 ||
                                       mpfr_get_prec (result) == 64))
        return 1; /* bug 16559 */
      if (strcmp (name, "exp10") == 0 && mpfr_get_prec (result) == 53)
        return 1; /* bug 16560 */
    }
  else if (err == 4) /* inconsistent overflow/result */
    {
#if 0 /* fixed in 2.20 (24 Mar 2014) */
      if (strcmp (name, "exp") == 0)  /* bug 16284 */
        return 1;
#endif
      if (strcmp (name, "exp10") == 0 && mpfr_get_prec (result) == 24)
        return 1; /* bug 16284 */
#if 0 /* fixed 27 Jun 2014 */
      if (strcmp (name, "y17") == 0 && (mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16561 */
      if (strcmp (name, "y42") == 0 && (mpfr_get_prec (result) == 53 ||
                                        mpfr_get_prec (result) == 64))
        return 1; /* bug 16562 */
#endif
#if 0 /* fixed 14 May 2014 */
      if (strcmp (name, "log1p") == 0 && mpfr_get_prec (result) == 64)
        return 1; /* bug 16564 */
#endif
    }
  else
    abort ();
#endif
#ifdef HAVE_LIBBF
    if (err == 1)
        return 1; /* libbf never sets errno */
#endif
  return 0;
}

static int
is_normal (mpfr_t x)
{
  if (mpfr_regular_p (x))
    {
      mpfr_prec_t prec = mpfr_get_prec (x);
      return mpfr_signbit (x)
        ? mpfr_cmp_si_2exp (x, -1, mpfr_get_emin () + prec - 2) <= 0
        : mpfr_cmp_si_2exp (x, +1, mpfr_get_emin () + prec - 2) >= 0;
    }
  else
    return 0;
}

static int
is_smallest_normal (mpfr_t x)
{
  if (mpfr_regular_p (x))
    {
      int sign = mpfr_signbit (x) ? -1 : 1;
      mpfr_prec_t prec = mpfr_get_prec (x);
      return mpfr_cmp_si_2exp (x, sign, mpfr_get_emin () + prec - 2) == 0;
    }
  else
    return 0;
}

static int
is_rndz (mpfr_t x, mpfr_rnd_t rnd)
{
  if (rnd == MPFR_RNDZ)
    return 1;
  if (mpfr_cmp_ui (x, 0) > 0 && rnd == MPFR_RNDD)
    return 1;
  if (mpfr_cmp_ui (x, 0) < 0 && rnd == MPFR_RNDU)
    return 1;
  return 0;
}

static int
is_rnda (mpfr_t x, mpfr_rnd_t rnd)
{
  if ((mpfr_cmp_ui (x, 0) > 0) && (rnd == MPFR_RNDU))
    return 1;
  if ((mpfr_cmp_ui (x, 0) < 0) && (rnd == MPFR_RNDD))
    return 1;
  return 0;
}

/* MPFR considers underflow *after* rounding, thus we have to recover underflow
 *before* rounding.
From Joseph Myers, 21 May 2014:
For functions with exactly defined results and exceptions (in particular, 
fma), glibc follows the architecture-specific definition of whether 
tininess is detected before or after rounding, since it's required to be 
the same for all operations with a binary result.  That is, after rounding 
on alpha, hppa, ia64, mips, sh, x86, before rounding on other supported 
architectures.
For other functions, the valid results from glibc are not intended to 
depend on that architecture-specific definition.  If the returned value is 
plus or minus the least normal, you should take it that underflow 
exceptions are permitted but not required, whether or not the result 
should be exact or "inexact" is raised (except that they are not permitted 
if the rounding is towards zero so the implied infinite precision result 
cannot possibly have underflowed).
*/
static int
fixed_underflow_p (mpfr_t result, int inex)
{
  int sign = mpfr_signbit (result) ? -1 : 1;

  /* check if |result| is the smallest normal number. The smallest subnormal
     is 1/2*2^emin, which is 2^(prec-1) smaller than the smallest normal
     number, which is thus 2^(emin+prec-2) */
  if (is_smallest_normal (result))
    if ((sign > 0 && inex > 0) || (sign < 0 && inex < 0))
        mpfr_set_underflow ();
  return mpfr_underflow_p ();
}

static void
dichotomy (mpfr_t x, int sign_x, mpfr_t y, mpcheck_func_t *ref)
{
  mpfr_t z, t;
  mpfr_exp_t emax;

  emax = mpfr_get_emax ();
  mpfr_set_emax (mpfr_get_emax_max ());
  mpfr_init2 (z, mpfr_get_prec (x));
  mpfr_init2 (t, mpfr_get_prec (x));
  while (1)
    {
      mpfr_add (z, x, y, MPFR_RNDN);
      mpfr_div_2exp (z, z, 1, MPFR_RNDN);
      if (mpfr_cmp (z, x) == 0 || mpfr_cmp (z, y) == 0)
        break;
      (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) ref->mpfr))
        (t, z, MPFR_RNDN);
      if (sign_x * mpfr_sgn (t) < 0)
        mpfr_set (y, z, MPFR_RNDN);
      else
        mpfr_set (x, z, MPFR_RNDN);
    }
  if (lrand48 () & 1)
    mpfr_set (x, y, MPFR_RNDN);
  mpfr_subnormalize (x, 0, MPFR_RNDN);
  mpfr_clear (z);
  mpfr_clear (t);
  mpfr_set_emax (emax);
}

#if !defined(HAVE_PARI)
int
is_accepted (const char *name, int numarg, mpfr_t op1, mpfr_t op2)
{
  return 1;
}

#endif

void
mpcheck (FILE *out, mp_exp_t e1, mp_exp_t e2, mp_exp_t e3,
	const char *const name,
	 void (*func) (void*,const void*, const void*, const void*))
{
  static int print_done = 0;
  int reduction_done;
  gmp_randstate_t state;
  mpfr_t op1, op2, op3, result, result_lib, result_more_prec, rmax, rlibmax;
  mpfr_t u, umax, umax_dir, op1max_dir, op2max_dir, op3max_dir, max_err_near, max_err_dir;
   mpfr_t op1max, op2max, op3max, range_min, range_max;
   mpfr_t xdplus, xdminus, last_x, last_result;
   void *rop1, *rop2, *rop3, *rresult;
   mpcheck_func_t *ref;
   unsigned long i, wrong_dir, wrong_range, wrong_monoton, wrong_symm;
   unsigned long wrong_errno, wrong_inexact, wrong, wrong_special;
   unsigned long wrong_overflow, wrong_underflow;
   mp_rnd_t rnd;
   int saved_errno, inex, rinex, overflow, roverflow, underflow, runderflow;

   for (i = 0 ; mpcheck_tab[i].name != NULL ; i++)
     if (strcmp (mpcheck_tab[i].name, name) == 0)
       break;
   if (mpcheck_tab[i].name == NULL)
     {
      fprintf (out, "Unknown function: %s\n", name);
      return;
     }
   ref = &mpcheck_tab[i];
   if (strcmp (function_to_check, "all") != 0 &&
       strcmp (ref->name, function_to_check) != 0)
     return;

   /* if seed is zero, we generate a random one */
   if (seed == 0)
     seed = getpid ();

   if (print_done == 0)
     print_done = fprintf (out,
			   "[precision=%lu, seed=%lu, emin=%ld, emax=%ld]\n",
			   (unsigned long) prec, seed,
			   (long) emin, (long) emax);

   /* Translate some absolute exponent to relative one's */
   if (e1 == LONG_MAX)     e1 = emax - 1;
   if (e2 == LONG_MAX)     e2 = emax - 1;
   if (e3 == LONG_MAX)     e3 = emax - 1;
   if (e1 == LONG_MAX-1)   e1 = emax / 2;
   if (e2 == LONG_MAX-1)   e2 = emax / 2;
   if (e3 == LONG_MAX-1)   e3 = emax / 2;
   if (e1 == LONG_MAX-2)   e1 = emin;
   if (e2 == LONG_MAX-2)   e2 = emin;
   if (e3 == LONG_MAX-2)   e3 = emin;
   if (e1 == LONG_MAX-3)   e1 = emax / 2 + 1;
   if (e2 == LONG_MAX-3)   e2 = emax / 2 + 1;
   if (e3 == LONG_MAX-3)   e3 = emax / 2 + 1;

   mpfr_inits2 (prec, op1, op2, op3, result, result_lib, u, umax, umax_dir,
                max_err_near, max_err_dir, op1max_dir, op2max_dir, op3max_dir,
                op1max, op2max, op3max, range_min, range_max, rmax, rlibmax,
	       xdplus, xdminus, last_x, last_result, NULL);
  mpfr_init2 (result_more_prec, prec+32);
  gmp_randinit (state, GMP_RAND_ALG_LC, 128);
  rop1 = (*new) (prec);
  rop2 = (*new) (prec);
  rop3 = (*new) (prec);
  rresult = (*new) (prec);

  mpfr_set_ui (max_err_near, 0, GMP_RNDN);
  mpfr_set_ui (max_err_dir, 0, GMP_RNDN);

  mpcheck_set_range (range_min, ref->min);
  mpcheck_set_range (range_max, ref->max);

  /* Check if func(e1) doesn't overflow/underflow */
  reduction_done = 0;
  int loop = 0;
  for (;;)
    {
      if (strcmp (ref->name, "fma") == 0)
        break;
      if (loop++ > 100)
        {
          fprintf (stderr, "Infinite loop in exponent range reduction\n");
          exit (1);
        }
      /* we should ensure that 2^(e1-1) <= op1 < 2^e1 below */
      mpfr_urandomb (op1, state);
      mpfr_set_exp (op1, 0);
      mpfr_mul_2si (op1, op1, e1, GMP_RNDN);
      /* we should ensure that 2^(e2-1) <= op1 < 2^e2 below */
      if (ref->NumArg >= 2)
        {
          mpfr_urandomb (op2, state);
          mpfr_set_exp (op2, 0);
          mpfr_mul_2si (op2, op2, e2, GMP_RNDN);
        }
      if (ref->NumArg >= 3)
        {
          mpfr_urandomb (op3, state);
          mpfr_set_exp (op3, 0);
          mpfr_mul_2si (op3, op3, e3, GMP_RNDN);
        }
      if (ref->NumArg == 1)
	inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
                  ref->mpfr)) (result, op1, GMP_RNDN);
      else if (ref->NumArg == 2)
	inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
                  ref->mpfr)) (result, op1, op2, GMP_RNDN);
      else /* 3 arguments */
	inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
                  ref->mpfr)) (result, op1, op2, op3, GMP_RNDN);
      inex = mpfr_subnormalize (result, inex, GMP_RNDN);

      if (mpfr_number_p (result))
	break;
      if (e1 > 0) e1 --; else e1++;
      if (e2 > 0) e2 --; else e2++;
      if (e3 > 0) e3 --; else e3++;
      reduction_done = 1;
    }
  if (reduction_done)
    {
      if (e1 > 0) e1 --; else e1++;
      if (e2 > 0) e2 --; else e2++;
      if (e3 > 0) e3 --; else e3++;
    }

   if (verbose >= 2)
     {
       if (ref->NumArg == 1)
	 fprintf (out, "Testing function %s for precision %lu, exponent %ld [seed=%lu]\n",
		  name, prec, e1, seed);
       else if (ref->NumArg == 2)
	 fprintf (out, "Testing function %s for precision %lu, exponents %ld and %ld [seed=%lu]\n",
		  name, prec, e1, e2, seed);
       else
	 fprintf (out, "Testing function %s for precision %lu, exponents %ld, %ld and %ld [seed=%lu]\n",
		  name, prec, e1, e2, e3, seed);
     }

  for (rnd = 0 ; rnd < 4; rnd++)
    {
      /* Skip rounding mode if not requested */
      if ((rnd_mode & (1 << rnd)) == 0)
	continue;
      if (!(*setrnd) ((mp_rnd_t) rnd))
	continue;
      mpcheck_rnd_mode = rnd;

      if (verbose >= 3)
	fprintf (out, " rounding mode %s:\n",mpfr_print_rnd_mode (rnd));
      /* reset the seed to test the same sequence of numbers with each
         rounding mode */
      gmp_randseed_ui (state, seed);
      srand48 (seed);

      mpfr_set_ui (umax, 0, GMP_RNDN);  /* (signed) maximal error in ulps */
      mpfr_set_ui (umax_dir, 0, GMP_RNDN); /* maximal error in ulps when 
					      wrong directed rounding */
      mpfr_set_ui (op2max, 0, GMP_RNDN);
      mpfr_set_ui (op2max_dir, 0, GMP_RNDN);

      mpfr_set_ui (op3max, 0, GMP_RNDN);
      mpfr_set_ui (op3max_dir, 0, GMP_RNDN);

      wrong = 0;
      wrong_dir = 0;     /* number of wrong (side of) directed roundings */
      wrong_range = 0;   /* number of wrong results wrt range        */
      wrong_monoton = 0; /* number of wrong results wrt monotonicity */
      wrong_symm = 0;    /* number of wrong results wrt symmetry     */
      wrong_errno = 0;   /* number of wrong results wrt errno     */
      wrong_inexact = 0; /* number of wrong inexact flags */
      wrong_special = 0; /* number of wrong special values */
      wrong_overflow = 0; /* number of wrong overflow exceptions */
      wrong_underflow = 0; /* number of wrong underflow exceptions */

      for (i = 0; i < N ; i++)
	{
          if (ref->NumArg == 1 && mpfr_regular_p (result) &&
              mpfr_sgn (last_result) * mpfr_sgn (result) < 0)
            dichotomy (op1, mpfr_sgn (result), last_x, ref);
          else
            /* Generate random numbers */
            gen (op1, e1, op2, e2, op3, e3, state, i, ref);
          assert (!mpfr_regular_p (op1) || mpfr_get_exp (op1) -
                  mpfr_min_prec (op1) >= mpfr_get_emin () - 1);

          if (mpfr_regular_p (result))
            {
              mpfr_set (last_x, op1, MPFR_RNDN);
              mpfr_set (last_result, result, MPFR_RNDN);
            }

          if (is_accepted (name, ref->NumArg, op1, op2) == 0)
            continue;

	  /* Convert the input to the format used by the library */
	  (*getfp) (rop1, op1);
          if (ref->NumArg >= 2)
            (*getfp) (rop2, op2);
          if (ref->NumArg >= 3)
            (*getfp) (rop3, op3);
	  /* Compute the result with MPFR and the LIBRARY */
          mpfr_clear_overflow ();
          mpfr_clear_underflow ();
          if (ref->NumArg == 1)
	    inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
                      ref->mpfr)) (result, op1, rnd);
	  else if (ref->NumArg == 2)
	    inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
                      ref->mpfr)) (result, op1, op2, rnd);
          else /* 3 arguments */
	    inex = (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
                      ref->mpfr)) (result, op1, op2, op3, rnd);
          inex = mpfr_subnormalize (result, inex, rnd);
          overflow = mpfr_overflow_p ();
          underflow = fixed_underflow_p (result, inex);
          errno = 0;
          if (fp_feclearexcept) {
              fp_feclearexcept();
          } else {
              feclearexcept (FE_ALL_EXCEPT); /* clear all exceptions */
          }
	  (*func) (rresult, rop1, rop2, rop3);
          if (fp_fetestexcept) {
              rinex = fp_fetestexcept (FE_INEXACT);
              roverflow = fp_fetestexcept (FE_OVERFLOW);
              runderflow = fp_fetestexcept (FE_UNDERFLOW);
          } else {
              rinex = fetestexcept (FE_INEXACT);
              roverflow = fetestexcept (FE_OVERFLOW);
              runderflow = fetestexcept (FE_UNDERFLOW);
          }
          saved_errno = errno;
	  /* Convert back the result to MPFR */
	  (*setfp) (result_lib, rresult);

          if (verbose >= 5)
            {
              printf ("op1="); mpfr_dump (op1);
              if (ref->NumArg >= 2)
                { printf ("op2="); mpfr_dump (op2); }
              printf ("result="); mpfr_dump (result);
              printf ("result_lib="); mpfr_dump (result_lib);
              printf ("errno: library %d\n", errno);
              printf ("inexact: mpfr %d, library %d\n", inex, rinex);
              printf ("overflow: mpfr %d, library %d\n", overflow, roverflow);
              printf ("underflow: mpfr %d, library %d\n", underflow, runderflow);
            }

          /* check if the library result gives Inf for rounding to zero
             and an inexact result */
          if (ref->NumArg == 1 && mpfr_inf_p (result_lib) && inex &&
              ((mpfr_sgn (result_lib) > 0 && rnd != MPFR_RNDU
                && rnd != MPFR_RNDN) || (mpfr_sgn (result_lib) < 0 &&
                                    rnd != MPFR_RNDD && rnd != MPFR_RNDN)) &&
              mpfr_regular_p (op1) && mpfr_cmp_ui (op1, 1) != 0 &&
              mpfr_cmp_si (op1, -1) != 0)
            {
              wrong_overflow ++;
              if (suppress (4, name, result, op1, op2))
                suppressed_wrong_overflow ++;
              else if (verbose >= 3 && wrong_overflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong overflow for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      overflow: mpfr %d, library %d\n",
                           overflow, roverflow);
                }
            }

          /* check if the library result gives 0 for rounding away:
             this is not considered as a real bug unless the error in ulps
             is large (checked elsewhere) or the underflow flag is not
             raised */
          if (ref->NumArg == 1 && mpfr_zero_p (result_lib) &&
              ((mpfr_sgn (result_lib) > 0 && rnd == MPFR_RNDU) ||
               (mpfr_sgn (result_lib) < 0 && rnd == MPFR_RNDD)) &&
              mpfr_regular_p (op1))
            {
              wrong_underflow ++;
              if (runderflow != 0)
                suppressed_wrong_underflow ++;
              else if (verbose >= 3 && wrong_underflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong underflow for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      underflow: mpfr %d, library %d\n",
                           underflow, runderflow);
                }
            }

          /* check if the library result gives MAX_FLOAT for rounding away,
             when MPFR gives Inf */
          if (ref->NumArg == 1 &&
              mpfr_inf_p (result) && !mpfr_inf_p (result_lib) &&
              ((mpfr_sgn (result_lib) > 0 && rnd != MPFR_RNDD
                && rnd != MPFR_RNDN) || (mpfr_sgn (result_lib) < 0 &&
                                    rnd != MPFR_RNDU && rnd != MPFR_RNDN)) &&
              mpfr_regular_p (op1) && mpfr_cmp_ui (op1, 1) != 0 &&
              mpfr_cmp_si (op1, -1) != 0)
            {
              wrong_overflow ++;
              if (verbose >= 3 && wrong_overflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong overflow for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      overflow: mpfr %d, library %d\n",
                           overflow, roverflow);
                }
            }

          /* check if the library result and underflow flag are inconsistent:
             (1) underflow flag is set and result exceeds smallest normal
                 (or equals smallest normal and rounding is toward zero)
             (2) underflow flag is not set and result is smaller than
                 smallest normal */
          if (ref->NumArg == 1 &&
              ((runderflow && is_normal (result_lib) && (!is_smallest_normal (result_lib) || is_rndz (result_lib, rnd))) ||
               (!runderflow && inex && (!is_normal (result_lib))
                && !mpfr_inf_p (result_lib))))
            {
              wrong_underflow ++;
              if (suppress (3, name, result, op1, op2))
                suppressed_wrong_underflow ++;
              else if (verbose >= 3 && wrong_underflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong underflow flag for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      underflow: mpfr %d, library %d\n",
                           underflow, runderflow);
                }
            }

          /* check if the library result and overflow flag are inconsistent:
             (1) overflow flag is set and result is below largest float
                 (or equals largest float and rounding is away from zero)
             (2) overflow flag is not set and result is Inf,
                 or largest float and rounding is to zero
                 (in case of inexact result) */
          if (ref->NumArg == 1 &&
              ((roverflow && mpfr_number_p (result_lib) && (!is_largest (result_lib) || is_rnda (result_lib, rnd))) ||
               (!roverflow && inex && (mpfr_inf_p (result_lib) ||
                (is_largest (result_lib) && is_rndz (result_lib, rnd))))))
            {
              wrong_overflow ++;
              if (suppress (4, name, result, op1, op2))
                suppressed_wrong_overflow ++;
              else if (verbose >= 3 && wrong_overflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong overflow flag for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      overflow: mpfr %d, library %d\n",
                           overflow, roverflow);
                }
            }

          /* check cases where the (inexact) result is the largest float,
             and rounding is towards zero */
          if (ref->NumArg == 1 && is_rndz (result_lib, rnd) &&
              mpfr_inf_p (result_lib) && inex)
            {
              wrong_overflow ++;
              if (suppress (4, name, result, op1, op2))
                suppressed_wrong_overflow ++;
              else if (verbose >= 3 && wrong_overflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong overflow flag for x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      overflow: mpfr %d, library %d\n",
                           overflow, roverflow);
                }
            }

          /* we can only compare the inexact flags when the results agree */
          if (((inex == 0 && rinex != 0) || (inex != 0 && inex == 0)) &&
              mpfr_cmp (result, result_lib) == 0)
            {
              wrong_inexact ++;
              if (suppress (0, name, result, op1, op2))
                suppressed_wrong_inexact ++;
              else if (verbose >= 3 && wrong_inexact <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong inexact flag: mpfr gives %d, "
                           "library %d\n", inex, rinex);
                  fprintf (out, "      x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  if (ref->NumArg >= 2)
                    {
                      fprintf (out, "      t=");
                      mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                    }
                  if (ref->NumArg >= 3)
                    {
                      fprintf (out, "      u=");
                      mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                    }
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, " \n");
                }
            }

          /* check overflow */
          if ((overflow == 0 && roverflow != 0) ||
              (overflow != 0 && overflow == 0))
            {
              wrong_overflow ++;
              if (suppress (4, name, result, op1, op2))
                suppressed_wrong_overflow ++;
              else if (verbose >= 3 && wrong_overflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong overflow flag: mpfr gives %d, "
                           "library %d\n", overflow, roverflow);
                  fprintf (out, "      x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  if (ref->NumArg >= 2)
                    {
                      fprintf (out, "      t=");
                      mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                    }
                  if (ref->NumArg >= 3)
                    {
                      fprintf (out, "      u=");
                      mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                    }
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, "\n      overflow: mpfr %d, library %d\n",
                           overflow, roverflow);
                }
            }

          /* check underflow flags agree (if results are identical) */
          if (((underflow == 0 && runderflow != 0) ||
               (underflow != 0 && underflow == 0)) &&
              mpfr_cmp (result, result_lib) == 0 &&
              /* we don't compare for the smallest normal, because of
                 before-rounding/after-rounding issues */
              !is_smallest_normal (result))
            {
              wrong_underflow ++;
              if (suppress (2, name, result, op1, op2))
                suppressed_wrong_underflow ++;
              else if (verbose >= 3 && wrong_underflow <= MAX_REPORTS)
                {
                  fprintf (out, "      wrong underflow flag: mpfr gives %d, "
                           "library %d\n", underflow, runderflow);
                  fprintf (out, "      x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  if (ref->NumArg >= 2)
                    {
                      fprintf (out, "      t=");
                      mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                    }
                  if (ref->NumArg >= 3)
                    {
                      fprintf (out, "      u=");
                      mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                    }
                  fprintf (out, "\n      library gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\n      mpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, " \n");
                }
            }

	  /* Check for correct rounding */
	  if (mpfr_cmp (result, result_lib) != 0)
	    {
              if (verbose >= 5) {
                printf("ERROR:\n");
                printf ("op1="); mpfr_dump (op1);
                if (ref->NumArg >= 2)
                  { printf ("op2="); mpfr_dump (op2); }
                printf ("result    ="); mpfr_dump (result);
                printf ("result_lib="); mpfr_dump (result_lib);
              }
	      /* Error ++ */
	      wrong ++;
              if (!mpfr_regular_p (result) || !mpfr_regular_p (result_lib))
                {
                  wrong_special ++;
                  break; /* ulp difference does not make sense */
                }
	      /* Recompute MPFR result with more prec */
	      if (ref->NumArg == 1)
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
                   ref->mpfr)) (result_more_prec, op1, rnd);
	      else if (ref->NumArg == 2)
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		   ref->mpfr)) (result_more_prec, op1, op2, rnd);
	      else /* three arguments */
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		   ref->mpfr)) (result_more_prec, op1, op2, op3, rnd);
	      /* Compute ULP */
	      mpcheck_ulp (u, result_lib, result_more_prec, prec);
	      if (rnd != GMP_RNDN)
		{
		  mp_rnd_t rnd2;
		  rnd2 = (rnd==GMP_RNDZ)
		    ? ((mpfr_sgn (result) > 0 ? GMP_RNDD : GMP_RNDU))
		    : rnd;
		  if ((rnd2 == GMP_RNDU && mpfr_sgn (u) < 0)
		      || (rnd2 == GMP_RNDD && mpfr_sgn (u) > 0))
		    {
		      wrong_dir ++;
		      if (mpfr_cmpabs (u, umax_dir) > 0)
			{
			  mpfr_set (umax_dir, u, GMP_RNDN);
			  mpfr_set (op1max_dir, op1, GMP_RNDN);
			  mpfr_set (op2max_dir, op2, GMP_RNDN);
			  mpfr_set (op3max_dir, op3, GMP_RNDN);
			}
		    }
		}
	      if (mpfr_cmpabs (u, umax) > 0)
		{
		  mpfr_set (umax, u, GMP_RNDN);
		  mpfr_set (op1max, op1, GMP_RNDN);
		  mpfr_set (op2max, op2, GMP_RNDN);
		  mpfr_set (op3max, op3, GMP_RNDN);
                  mpfr_set (rmax, result, GMP_RNDN);
                  mpfr_set (rlibmax, result_lib, GMP_RNDN);
		}
	    }

	  /* Testing if the range is verified */
	  if ((test & MPCHECK_TEST_RANGE) != 0)
	    {
	      if (mpfr_cmp (result_lib, range_min) < 0
		  || mpfr_cmp (result_lib, range_max) > 0)
		{
		  wrong_range++;
                  /* for glibc it is not a bug, unless the error is large,
                     which is checked elsewhere */
                  if (1)
                    suppressed_wrong_range ++;
		  else if (wrong_range <= MAX_REPORTS && verbose >= 3)
		    {
		      fprintf (out, "      wrong range for x=");
		      mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                      if (ref->NumArg >= 2)
                        {
                          fprintf (out, " t=");
                          mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                        }
                      if (ref->NumArg >= 3)
                        {
                          fprintf (out, " u=");
                          mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                        }
		      fprintf (out, "\n           f(x)=");
		      mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
		      fprintf (out, "\n      not between ");
		      mpfr_out_str (out, OUT, 0, range_min, GMP_RNDN);
		      fprintf (out, "      and ");
		      mpfr_out_str (out, OUT, 0, range_max, GMP_RNDN);
		      fprintf (out, " \n");
		    }
		}
	    }

	  /* Testing if the monotonicity is verified */
	  if (ref->NumArg == 1
	      && (test & MPCHECK_TEST_MONOTON)
	      && (ref->monoton != NO_MONOTON))
	    {
	      mpfr_set (xdminus, op1, GMP_RNDN);
	      mpfr_nextbelow (xdminus);
              if (is_accepted (name, 1, xdminus, NULL) == 0)
                continue;
	      (*getfp) (rop1, xdminus);
	      (*func) (rresult, rop1, rop2, rop3);
	      (*setfp) (xdminus, rresult);
	      
	      mpfr_set (xdplus, op1, GMP_RNDN);
	      mpfr_nextabove (xdplus);
              if (is_accepted (name, 1, xdplus, NULL) == 0)
                continue;
              (*getfp) (rop1, xdplus);
              (*func) (rresult, rop1, rop2, rop3);
              (*setfp) (xdplus, rresult);

	      if (ref->monoton == INCREASING)
		{
		  if (mpfr_cmp (xdminus, result_lib) > 0)
		    {
		      wrong_monoton++;
		      if (wrong_monoton <= MAX_REPORTS && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n            f(x-)=");
			  mpfr_out_str (out, OUT, 0, xdminus, GMP_RNDN);
			  fprintf (out, "\n      not <= f(x)=");
			  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		  if (mpfr_cmp (xdplus, result_lib) < 0)
		    {
		      wrong_monoton++;
		      if (wrong_monoton <= MAX_REPORTS && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n              f(x)=");
			  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      not <= f(x+)=");
			  mpfr_out_str (out, OUT, 0, xdplus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		}
	      else  /* monoton == DECREASING */
		{
		  if (mpfr_cmp (xdminus, result_lib) < 0)
		    {
		      wrong_monoton++;
		      if (wrong_monoton <= MAX_REPORTS && verbose >= 3)
			{
			  fprintf (out,
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n            f(x-)=");
			  mpfr_out_str (out, OUT, 0, xdminus, GMP_RNDN);
			  fprintf (out, "\n      not >= f(x)=");
                          mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		  if (mpfr_cmp (xdplus, result_lib) > 0)
		    {
		      wrong_monoton++;
		      if (wrong_monoton <= MAX_REPORTS && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
                          mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n              f(x)=");
                          mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      not >= f(x+)=");
                          mpfr_out_str (out, OUT, 0, xdplus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		}
	    }
	  
	  /* Testing if the symmetry is verified */
	  if (ref->NumArg == 1 && rnd == GMP_RNDN
              && (test & MPCHECK_TEST_SYMM)!=0 && (ref->symm != NO_SYMM))
	    {
	      mpfr_neg (xdminus, op1, GMP_RNDN);
              if (is_accepted (name, 1, xdminus, NULL) == 0)
                continue;
              (*getfp) (rop1, xdminus);
              (*func) (rresult, rop1, rop2, rop3);
              (*setfp) (xdminus, rresult);

	      if (ref->symm == ODD)
		{
		  if (mpfr_cmp3 (xdminus, result_lib, -1) != 0)
		    {
		      wrong_symm++;
		      if (wrong_symm <= MAX_REPORTS && verbose >= 3)
			{
			  fprintf (out, "      symmetry not respected for x=");
                          mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n           f(x)= ");
			  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      and f(-x)=");
			  mpfr_out_str (out, OUT, 0, xdminus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		}
	      else /* symm == EVEN */
		{ 
		  if (mpfr_cmp (result_lib, xdminus) != 0) 
		    {
		      wrong_symm++;
		      if (wrong_symm <= MAX_REPORTS && verbose >= 3) 
			{
			  fprintf (out, "      symmetry not respected for x=");
			  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
			  fprintf (out, "\n           f(x)=");
                          mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      and f(-x)=");
                          mpfr_out_str (out, OUT, 0, xdminus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		    }
		}
	    }

          /* Testing errno */
          if (saved_errno != 0)
            {
              if (saved_errno == EDOM)
                {
                  /* result and result_lib should be NaN */
                  if (mpfr_nan_p (result) == 0 || mpfr_nan_p (result_lib) == 0)
                    {
                     wrong_errno ++;
                     if (suppress (1, name, result, op1, op2))
                       suppressed_wrong_errno ++;
                     else if (wrong_errno <= MAX_REPORTS && verbose >= 3) 
                        {
                          fprintf (out, "errno=%d: out of domain of function\n",
                                   saved_errno);
                          fprintf (out, "x=");
                          mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                          if (ref->NumArg >= 2)
                            {
                              fprintf (out, " t=");
                              mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                            }
                          if (ref->NumArg >= 3)
                            {
                              fprintf (out, " u=");
                              mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                            }
                        fprintf (out, "\nlibrary gives ");
                        mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                        fprintf (out, "\nmpfr    gives ");
                        mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                        fprintf (out, " \n");
                        }
                    }
                }
              else if (saved_errno == ERANGE) /* 34 */
                {
                  if (mpfr_nan_p (result) || mpfr_nan_p (result_lib) ||
                      mpfr_cmp (result, result_lib))
                    {
                      wrong_errno ++;
                      if (wrong_errno <= MAX_REPORTS && verbose >= 3)
                        {
                        fprintf (out, "errno=%d: result not representable\n",
                                 saved_errno);
                        fprintf (out, "x=");
                        mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                        if (ref->NumArg >= 2)
                          {
                            fprintf (out, " t=");
                            mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                          }
                        if (ref->NumArg >= 3)
                          {
                            fprintf (out, " u=");
                            mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                          }
                        fprintf (out, "\nlibrary gives ");
                        mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                        fprintf (out, "\nmpfr    gives ");
                        mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                        fprintf (out, " \n");
                        }
                    }
                }
              else
                {
                 wrong_errno ++;
                 if (wrong_errno <= MAX_REPORTS && verbose >= 3)
                  {
                  fprintf (out, "errno=%d\n", saved_errno);
                  fprintf (out, "x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  if (ref->NumArg >= 2)
                    {
                    fprintf (out, " t=");
                    mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                    }
                  if (ref->NumArg >= 3)
                    {
                    fprintf (out, " u=");
                    mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                    }
                  fprintf (out, "\nlibrary gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\nmpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, " \n");
                  }
                }
            }
          if (mpfr_nan_p (result) && !mpfr_nan_p (op1) &&
              (ref->NumArg < 2 || !mpfr_nan_p (op2)))
            /* check errno is set */
            {
              if (saved_errno == 0)
                {
                  wrong_errno++;
                  if (suppress (1, name, result, op1, op2))
                    suppressed_wrong_errno ++;
                  else
                {
                 if (wrong_errno <= MAX_REPORTS && verbose >= 3)
                  {
                  fprintf (out, "Result is NaN but errno is not set\n");
                  fprintf (out, "x=");
                  mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                  if (ref->NumArg >= 2)
                    {
                    fprintf (out, " t=");
                    mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                    }
                  if (ref->NumArg >= 3)
                    {
                    fprintf (out, " u=");
                    mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                    }
                  fprintf (out, "\nlibrary gives ");
                  mpfr_out_str (out, OUT, 0, result_lib, GMP_RNDN);
                  fprintf (out, "\nmpfr    gives ");
                  mpfr_out_str (out, OUT, 0, result, GMP_RNDN);
                  fprintf (out, " \n");
                  }
                }
                }
              else if (saved_errno != EDOM)
                {
                  wrong_errno ++;
                  if (wrong_errno <= MAX_REPORTS && verbose >= 3)
                    {
                      fprintf (out, "errno=%d\n", saved_errno);
                      fprintf (out, "x=");
                      mpfr_out_str (out, OUT, 0, op1, GMP_RNDN);
                      if (ref->NumArg >= 2)
                        {
                          fprintf (out, " t=");
                          mpfr_out_str (out, OUT, 0, op2, GMP_RNDN);
                        }
                      if (ref->NumArg >= 3)
                        {
                          fprintf (out, " u=");
                          mpfr_out_str (out, OUT, 0, op3, GMP_RNDN);
                        }
                      fprintf (out, " \n");
                    }
                }
            }
	} /* for i to N */

      tot_wrong += wrong;
      if (strcmp (ref->name, "add") == 0 || strcmp (ref->name, "sub") == 0 ||
          strcmp (ref->name, "mul") == 0 || strcmp (ref->name, "div") == 0 ||
          strcmp (ref->name, "fma") == 0 || strcmp (ref->name, "sqrt") == 0)
        tot_wrong_basic += wrong;
      tot_wrong_dir += wrong_dir;
      tot_wrong_range += wrong_range;
      tot_wrong_monoton += wrong_monoton;
      tot_wrong_symm += wrong_symm;
      tot_wrong_errno += wrong_errno;
      tot_wrong_inexact += wrong_inexact;
      tot_wrong_overflow += wrong_overflow;
      tot_wrong_underflow += wrong_underflow;
      
      /* if any difference occurred, prints the maximal ulp error */
      if (mpfr_cmp_ui (umax, 0) != 0.0 && verbose >= 3)
	{
	  fprintf (out, "      ");
	  mpfr_out_str (out, 10, 3, umax, MPFR_RNDA);
	  fprintf (out, " ulp(s) for x=");
          out_value (out, op1max);
	  if (ref->NumArg >= 2)
	    {  
	      fprintf (out, " t=");
              out_value (out, op2max);
	    }
	  if (ref->NumArg >= 3)
	    {  
	      fprintf (out, " u=");
              out_value (out, op3max);
	    }
	  fprintf (out, "\n");
          if (verbose >= 4)
            {
              fprintf (out, "         mpfr gives ");
              out_value (out, rmax);
              fprintf (out, "\n");
              fprintf (out, "         lib  gives ");
              out_value (out, rlibmax);
              fprintf (out, "\n");
            }
	}
      
      if (wrong_dir != 0 && verbose >= 3)
	{
	  fprintf (out, "      wrong directed rounding for x=");
	  mpfr_out_str (out, OUT, 0, op1max_dir, GMP_RNDN);
	  if (ref->NumArg >= 2)
	    {
	      fprintf (out, " t=");
	      mpfr_out_str (out, OUT, 0, op2max_dir, GMP_RNDN);
	    }
	  if (ref->NumArg >= 3)
	    {
	      fprintf (out, " u=");
	      mpfr_out_str (out, OUT, 0, op3max_dir, GMP_RNDN);
	    }
	  fprintf (out, " [");
	  mpfr_out_str (out, 10, 3, umax_dir, MPFR_RNDA);
	  fprintf (out, "]\n");
	}
      
      mpfr_abs (umax, umax, GMP_RNDN);
      
      if (verbose >= 3)
	{
	  if (rnd == GMP_RNDN)
            {
              if (wrong_range + wrong_monoton + wrong_symm + wrong_errno + wrong_inexact)
                fprintf(out, 
                        "   nb errors range/monotonicity/symmetry/errno/inexact: %lu/%lu/%lu/%lu/%lu\n", 
                        wrong_range, wrong_monoton, wrong_symm, wrong_errno,
                        wrong_inexact);
            }
	  else
            {
              if (wrong_range + wrong_monoton + wrong_errno + wrong_inexact)
                fprintf(out, 
                        "   nb errors range/monotonicity/errno/inexact: %lu/%lu/%lu/%lu\n", 
                        wrong_range, wrong_monoton, wrong_errno, wrong_inexact);
            }
	}
      if (verbose >= 3)
	{
          if (wrong + wrong_dir > 0 || !mpfr_zero_p (umax))
            {
              fprintf (out, 
                       "   nb errors/wrong directed/max ulp diff: %lu/%lu/", 
                       wrong, wrong_dir);
              mpfr_out_str (out, 10, 3, umax, MPFR_RNDA);
              fprintf (out, "\n");
            }
	}

      if (rnd == GMP_RNDN)
	{
	  if (mpfr_cmp (umax, max_err_near) > 0)
	    mpfr_set (max_err_near, umax, GMP_RNDN);
	}
      else
	{
	  if (mpfr_cmp (umax, max_err_dir) > 0)
            mpfr_set (max_err_dir, umax, GMP_RNDN);
	}
      fflush (out);
    } /* For each rnd */

  fprintf (out, "Max. errors for %s [exp. %ld]: ", name, e1);
  mpfr_out_str (out, 10, 3, max_err_near, MPFR_RNDA);
  fprintf (out, " (nearest), ");
  mpfr_out_str (out, 10, 3, max_err_dir, MPFR_RNDA);
  fprintf (out, " (directed)\n");
  if (verbose >= 3)
    fprintf (out, "\n");

  if (mpfr_cmp (max_err_near, mpcheck_max_err_near) > 0)
    mpfr_set (mpcheck_max_err_near, max_err_near, GMP_RNDN);
  
  if (mpfr_cmp (max_err_dir, mpcheck_max_err_dir) > 0)
    mpfr_set (mpcheck_max_err_dir, max_err_dir, GMP_RNDN);
  
  (*del) (rop1);
  (*del) (rop2);
  (*del) (rop3);
  (*del) (rresult);
  gmp_randclear (state);
  mpfr_clears (op1, op2, op3, result, result_lib, u, umax, umax_dir,
               max_err_near, max_err_dir, op1max_dir, op2max_dir, op3max_dir,
	       op1max, op2max, op3max, range_min, range_max, rmax, rlibmax,
               xdplus, xdminus, last_x, last_result, NULL);
  mpfr_clear (result_more_prec);
  mpfr_free_cache ();
  (*setrnd) (GMP_RNDN);
}

void 
mpcheck_check (FILE *out, mpcheck_user_func_t *tab)
{
  int i;
  for (i = 0 ; tab[i].name != NULL ; i++)
    mpcheck (out, tab[i].e1, tab[i].e2, tab[i].e3, tab[i].name, tab[i].func);
}

