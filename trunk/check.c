/*
 * MPCHECK - Check mathematical functions
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

#include "mpcheck.h"

mpfr_t  mpcheck_max_err_dir, mpcheck_max_err_near;
mp_rnd_t mpcheck_rnd_mode;

static void
mpcheck_set_range (mpfr_t dest, mpcheck_range_e range)
{
  switch (range)
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
    case RANGE_PI2:
      mpfr_const_pi (dest, GMP_RNDU);
      mpfr_div_2ui (dest, dest, 1, GMP_RNDU);
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
    case -RANGE_PI2:
      mpfr_const_pi (dest, GMP_RNDU);
      mpfr_div_2ui (dest, dest, 1, GMP_RNDU);
      mpfr_neg (dest, dest, GMP_RNDD);
      break;
    default:
      fprintf (stderr, "MPCHECK ERROR: Range undefined\n");
      abort ();
    }
}

static void
mpcheck_ulp (mpfr_t ulp, mpfr_t library, mpfr_t reference, mp_prec_t prec)
{
  mpfr_sub (ulp, library, reference, GMP_RNDN);
  if (mpfr_cmp_ui (ulp, 0) == 0)
    {
      fprintf (stderr, "MPCHECK ERROR: Ulp couldn't be 0!\n");
      fprintf (stderr, "LIBRARY="); 
      mpfr_out_str (stderr, 10, 0, library, GMP_RNDN);
      fprintf (stderr, "\nREFEREN="); 
      mpfr_out_str (stderr, 10, 0, reference, GMP_RNDN);
      fprintf (stderr, "\n");
      abort ();
    }
  if (mpfr_number_p (ulp))
    mpfr_set_exp (ulp, mpfr_get_exp (ulp) - mpfr_get_exp (reference) + prec);
}

static void 
usage (void)
{
  fprintf (stderr, 
	   "Usage: mpcheck [options]\n"
	   "Usage: mpcheck [options] <function> <exponent>\n"
	   "Usage: mpcheck [options] <function> <exp1> <exp2>\n"
	   "where options are:\n"
	   "--seed=s    : set random seed to s [default 1]\n"
	   "--verbose=k : set verbose level to k [default 3]\n"
	   "--range=bool: test output range\n"
	   "--mono=bool : test monotonicity\n"
	   "--symm=bool : test symmetry\n"
	   "--RNDN=bool : test Round To Nearest\n"
	   "--RNDZ=bool : test Round To Zero\n"
	   "--RNDD=bool : test Round To -Inf\n"
	   "--RNDU=bool : test Round To +Inf\n"
	   "--prec=val  : Set prec\n"
	   "--num=n     : N\n");
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
static int verbose;

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
	  else
	    {
	      fprintf (stderr, "Invalid option: %s\n", argv[i]);
	      usage ();
	    }
	}
    }

  /* Init */
  mpfr_inits2 (prec, mpcheck_max_err_dir, mpcheck_max_err_near, NULL);
  mpfr_set_ui (mpcheck_max_err_dir, 0, GMP_RNDN);
  mpfr_set_ui (mpcheck_max_err_near, 0, GMP_RNDN);
  mpfr_set_emin (emin);
  mpfr_set_emax (emax);  
}

void
mpcheck_clear (FILE *out)
{
  fprintf (out, "Max. errors : ");
  mpfr_out_str (out, 10, 0,  mpcheck_max_err_near, GMP_RNDN);
  fprintf (out, " (nearest), ");
  mpfr_out_str (out, 10, 0, mpcheck_max_err_dir, GMP_RNDN);
  fprintf (out, " (directed)\n");
 
  mpfr_clears (mpcheck_max_err_dir, mpcheck_max_err_near, NULL);
}

void 
mpcheck (FILE *out, mp_exp_t e1, mp_exp_t e2,
	const char *const name, 
	 void (*func) (void*,const void*, const void*))
{
  static int print_done = 0;
  gmp_randstate_t state;
   mpfr_t op1, op2, result, result_lib, result_more_prec;
   mpfr_t u, umax, umax_dir, op1max_dir, op2max_dir, max_err_near, max_err_dir;
   mpfr_t op1max, op2max, range_min, range_max;
   mpfr_t xdplus, xdminus;
   void *rop1, *rop2, *rresult;
   mpcheck_func_t *ref;
   unsigned long i, wrong, wrong_range, wrong_monoton, wrong_symm, tot;
   mp_rnd_t rnd;

   for (i = 0 ; mpcheck_tab[i].name != NULL ; i++)
     if (strcmp (mpcheck_tab[i].name, name) == 0)
       break;
   if (mpcheck_tab[i].name == NULL)     
     {
      fprintf (out, "Unknown function: %s\n", name);
      return;
     }
   ref = &mpcheck_tab[i];
   
   if (print_done == 0)
     print_done = fprintf (out, 
			   "[precision=%lu, seed=%lu, emin=%ld, emax=%ld]\n", 
			   (unsigned long) prec, seed,
			   (long) emin, (long) emax);

   if (e1 == LONG_MAX)     e1 = emax - 1;
   if (e2 == LONG_MAX)     e2 = emax -1; 
   if (e1 == LONG_MAX-1)   e1 = emax /2;
   if (e2 == LONG_MAX-1)   e2 = emax /2;
   if (e1 == LONG_MAX-2)   e1 = emin;
   if (e2 == LONG_MAX-2)   e2 = emin;

  mpfr_inits2 (prec, op1, op2, result, result_lib, u, umax, umax_dir, 
	       max_err_near, max_err_dir, op1max_dir, op2max_dir, 
	       op1max, op2max, range_min, range_max,
	       xdplus, xdminus, NULL);
  mpfr_init2 (result_more_prec, prec+32);
  gmp_randinit (state, GMP_RAND_ALG_LC, 128);
  rop1 = (*new) (prec);
  rop2 = (*new) (prec);
  rresult = (*new) (prec);

  mpfr_set_ui (max_err_near, 0, GMP_RNDN);
  mpfr_set_ui (max_err_dir, 0, GMP_RNDN);

  mpcheck_set_range (range_min, ref->min);
  mpcheck_set_range (range_max, ref->max);

  /* Check if func(e1) doesn't overflow/underflow */
  for (;;)
    {
      /* mpfr_set_ui_2exp (op1, 1, e1, GMP_RNDN);
      mpfr_sub_ui (op1, op1, 1, GMP_RNDN);
      mpfr_set_ui_2exp (op2, 1, e2, GMP_RNDN);
      mpfr_sub_ui (op2, op2, 1, GMP_RNDN); --> Gamma*/

      mpfr_urandomb (op1, state);
      mpfr_mul_2si (op1, op1, e1, GMP_RNDN);
      mpfr_urandomb (op2, state);
      mpfr_mul_2si (op2, op2, e2, GMP_RNDN); 
      
      if (ref->NumArg == 2)
	(*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
	   ref->mpfr)) (result, op1, op2, GMP_RNDN);
      else
	(*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
	   ref->mpfr)) (result, op1, GMP_RNDN);
      if (mpfr_number_p (result))
	break;
      if (e1 > 0) e1 --; else e1++;
      if (e2 > 0) e2 --; else e2++;
    }

   if (verbose >= 2)
     {
       if (ref->NumArg == 1)
	 fprintf (out, "Testing function %s for exponent %ld.\n", 
		  name, e1);
       else
	 fprintf (out, "Testing function %s for exponents %ld and %ld.\n", 
		  name, e1, e2);
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

      mpfr_set_ui (umax, 0, GMP_RNDN);  /* (signed) maximal error in ulps */
      mpfr_set_ui (umax_dir, 0, GMP_RNDN); /* maximal error in ulps when 
					      wrong directed rounding */
      mpfr_set_ui (op2max, 0, GMP_RNDN);
      mpfr_set_ui (op2max_dir, 0, GMP_RNDN);

      tot = 0;
      wrong = 0;         /* number of wrong directed roundings */
      wrong_range = 0;   /* number of wrong results wrt range        */
      wrong_monoton = 0; /* number of wrong results wrt monotonicity */
      wrong_symm = 0;    /* number of wrong results wrt symmetry     */

      for (i = 0; i < N ; i++)
	{
	  /* Generate random numbers */
 	  mpfr_urandomb (op1, state);
	  mpfr_mul_2si (op1, op1, e1, GMP_RNDN);
	  mpfr_urandomb (op2, state);
	  mpfr_mul_2si (op2, op2, e2, GMP_RNDN);
	  /* Convert the input to the format used by the library */
	  (*getfp) (rop1, op1);
	  (*getfp) (rop2, op2);
	  /* Compute the result with MPFR and the LIBRARY */
	  if (ref->NumArg == 2)
	    {
	      (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		 ref->mpfr)) (result, op1, op2, rnd); 
	      (*func) (rresult, rop1, rop2);
	    }
	  else
	    {
	      (*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
                 ref->mpfr)) (result, op1, rnd);
              (*func) (rresult, rop1, rop2);
	    }
	  /* Convert back the result to MPFR */
	  (*setfp) (result_lib, rresult);

	  /* Check for correct rounding */
	  if (mpfr_cmp (result, result_lib) != 0)
	    {
	      /* Recompute MPFR result with more prec */
	      if (ref->NumArg == 2)
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		   ref->mpfr)) (result_more_prec, op1, op2, rnd);
	      else
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t))
                   ref->mpfr)) (result_more_prec, op1, rnd);
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
		      wrong ++;
		      if (mpfr_cmpabs (u, umax_dir) > 0)
			{
			  mpfr_set (umax_dir, u, GMP_RNDN);
			  mpfr_set (op1max_dir, op1, GMP_RNDN);
			  mpfr_set (op2max_dir, op2, GMP_RNDN);
			}
		    }
		}
	      /* Error ++ */
	      tot ++;
	      if (mpfr_cmpabs (u, umax) > 0)
		{
		  mpfr_set (umax, u, GMP_RNDN);
		  mpfr_set (op1max, op1, GMP_RNDN);
		  mpfr_set (op2max, op2, GMP_RNDN);
		}
	    }

	  /* Testing if the range is verified */
	  if ((test&MPCHECK_TEST_RANGE) != 0)
	    {
	      if (mpfr_cmp (result_lib, range_min) < 0
		  || mpfr_cmp (result_lib, range_max) > 0)
		{
		  if (wrong_range == 0 && verbose >= 3)
		    {
		      fprintf (out, "      outside range for x=");
		      mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
		      fprintf (out, "\n           f(x)=");
		      mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
		      fprintf (out, "\n      not between ");
		      mpfr_out_str (out, 10, 0, range_min, GMP_RNDN);
		      fprintf (out, "      and ");
		      mpfr_out_str (out, 10, 0, range_max, GMP_RNDN);
		      fprintf (out, " \n");
		    }
		  wrong_range++;
		}
	    }

	  /* Testing if the monotonicity is verified */
	  if (ref->NumArg == 1
	      && (test & MPCHECK_TEST_MONOTON)
	      && (ref->monoton != NO_MONOTON))
	    {
	      mpfr_set (xdminus, op1, GMP_RNDN);
	      mpfr_nextbelow (xdminus);
	      (*getfp) (rop1, xdminus);
	      (*func) (rresult, rop1, rop2);
	      (*setfp) (xdminus, rresult);
	      
	      mpfr_set (xdplus, op1, GMP_RNDN);
	      mpfr_nextabove (xdplus);
              (*getfp) (rop1, xdplus);
              (*func) (rresult, rop1, rop2);
              (*setfp) (xdplus, rresult);

	      if (ref->monoton == INCREASING)
		{
		  if (mpfr_cmp (xdminus, result_lib) > 0)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n            f(x-)=");
			  mpfr_out_str (out, 10, 0, xdminus, GMP_RNDN);
			  fprintf (out, "\n      not <= f(x)=");
			  mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_monoton++;
		    }
		  if (mpfr_cmp (xdplus, result_lib) < 0)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n              f(x)=");
			  mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      not <= f(x+)=");
			  mpfr_out_str (out, 10, 0, xdplus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_monoton++;
		    }
		}
	      else  /* monoton == DECREASING */
		{
		  if (mpfr_cmp (xdminus, result_lib) < 0)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  fprintf (out,
				   "      monotonicity not respected for x=");
			  mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n            f(x-)=");
			  mpfr_out_str (out, 10, 0, xdminus, GMP_RNDN);
			  fprintf (out, "\n      not >= f(x)=");
                          mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_monoton++;
		    }
		  if (mpfr_cmp (xdplus, result_lib) > 0)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  fprintf (out, 
				   "      monotonicity not respected for x=");
                          mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n              f(x)=");
                          mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      not >= f(x+)=");
                          mpfr_out_str (out, 10, 0, xdplus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_monoton++;
		    }
		}
	    }
	  
	  /* Testing if the symmetry is verified */
	  if (ref->NumArg == 1 && rnd == GMP_RNDN
              && (test & MPCHECK_TEST_SYMM)!=0 && (ref->symm != NO_SYMM))
	    {
	      mpfr_neg (xdminus, op1, GMP_RNDN);
              (*getfp) (rop1, xdminus);
              (*func) (rresult, rop1, rop2);
              (*setfp) (xdminus, rresult);

	      if (ref->symm == ODD)
		{
		  if (mpfr_cmp3 (xdminus, result_lib, -1) != 0)
		    {
		      if (wrong_symm == 0 && verbose >= 3)
			{
			  fprintf (out, "      symmetry not respected for x=");
                          mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n           f(x)= ");
			  mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      and f(-x)=");
			  mpfr_out_str (out, 10, 0, xdminus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_symm++;
		    }
		}
	      else /* symm == EVEN */
		{ 
		  if (mpfr_cmp (result_lib, xdminus) != 0) 
		    {
		      if (wrong_symm == 0 && verbose >= 3) 
			{
			  fprintf (out, "      symmetry not respected for x=");
			  mpfr_out_str (out, 10, 0, op1, GMP_RNDN);
			  fprintf (out, "\n           f(x)=");
                          mpfr_out_str (out, 10, 0, result_lib, GMP_RNDN);
			  fprintf (out, "\n      and f(-x)=");
                          mpfr_out_str (out, 10, 0, xdminus, GMP_RNDN);
			  fprintf (out, " \n");
			}
		      wrong_symm++;
		    }
		}
	    }
	} /* for i to N */
      
      /* if any difference occured, prints the maximal ulp error */
      if (mpfr_cmp_ui (umax, 0) != 0.0 && verbose >= 3)
	{
	  fprintf (out, "      ");
	  mpfr_out_str (out, 10, 0, umax, GMP_RNDN);
	  fprintf (out, " ulp(s) for x=");
	  mpfr_out_str (out, 10, 0, op1max, GMP_RNDN);
	  if (ref->NumArg == 2)
	    {  
	      fprintf (out, " t=");
	      mpfr_out_str (out, 10, 0, op2max, GMP_RNDN);
	    }
	  fprintf (out, "\n");
	}
      
      if (wrong != 0 && verbose >= 3)
	{
	  fprintf (out, "      wrong directed rounding for x=");
	  mpfr_out_str (out, 10, 0, op1max_dir, GMP_RNDN);
	  if (ref->NumArg == 2)
	    {
	      fprintf (out, " t=");
	      mpfr_out_str (out, 10, 0, op2max_dir, GMP_RNDN);
	    }
	  fprintf (out, " [");
	  mpfr_out_str (out, 10, 0, umax_dir, GMP_RNDN);
	  fprintf (out, "]\n");
	}
      
      mpfr_abs (umax, umax, GMP_RNDN);
      
      if (verbose >= 3 && ref->NumArg == 1)
	{
	  if (rnd == GMP_RNDN)
	    fprintf(out, 
		    "   nb errors range/monotonicity/symmetry: %lu/%lu/%lu\n", 
		    wrong_range, wrong_monoton, wrong_symm);
	  else
	    fprintf(out, 
		    "   nb errors range/monotonicity: %lu/%lu\n", 
		    wrong_range, wrong_monoton);
	}
      if (verbose >= 3)
	{
	  fprintf (out, 
		   "   nb errors/wrong directed/max ulp diff: %lu/%lu/", 
		   tot, wrong);
	  mpfr_out_str (out, 10, 0, umax, GMP_RNDN);
	  fprintf (out, "\n");
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
  mpfr_out_str (out, 10, 0, max_err_near, GMP_RNDN);
  fprintf (out, " (nearest), ");
  mpfr_out_str (out, 10, 0, max_err_dir, GMP_RNDN);
  fprintf (out, " (directed)\n");
  if (verbose >= 3)
    fprintf (out, "\n");
  
  if (mpfr_cmp (max_err_near, mpcheck_max_err_near) > 0)
    mpfr_set (mpcheck_max_err_near, max_err_near, GMP_RNDN);
  
  if (mpfr_cmp (max_err_dir, mpcheck_max_err_dir) > 0)
    mpfr_set (mpcheck_max_err_dir, max_err_dir, GMP_RNDN);
  
  (*del) (rop1);
  (*del) (rop2);
  (*del) (rresult);
  gmp_randclear (state);
  mpfr_clears (op1, op2, result, result_lib, u, umax, umax_dir,
               max_err_near, max_err_dir, result_more_prec, 
	       op1max_dir, op2max_dir, op1max, op2max, range_min, 
	       range_max, xdplus, xdminus, NULL);
  (*setrnd) (GMP_RNDN);
}

void 
mpcheck_check (FILE *out, mpcheck_user_func_t *tab)
{
  int i;
  for (i = 0 ; tab[i].name != NULL ; i++)
    mpcheck (out, tab[i].e1, tab[i].e2, tab[i].name, tab[i].func);
}
 