/*
 * MPCHECK - Check double LIBM functions
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


/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

#include "mpcheck.h"

void
test (const char *func,
      mp_exp_t e1, mp_exp_t e2,
      unsigned long N, unsigned long seed)
{
   unsigned long i, wrong, wrong_range, wrong_monoton, wrong_symm, tot;
   mpfr_t op1, op2, z, rmpfr, rmpfrref;
   fptype op1d, op2d, rlibmd, rmpfrd;
   fptype xmax, xmax_dir, tmax, tmax_dir, range_min, range_max;
   double u, umax, umax_dir, max_err_near, max_err_dir;
   int monoton, symm, num;
   gmp_randstate_t state;
   mp_rnd_t rnd;
   fptype (*testfun_libm) (fptype);
   int (*testfun_mpfr) (mpfr_ptr, mpfr_srcptr, mp_rnd_t);

   for (i = 0 ; i < FuncTabNum ; i++)
     if (strcmp (FuncTab[i].name, func) == 0)
       break;
   if (i == FuncTabNum)     
     {
      fprintf (stderr, "Unknown function: %s\n", func);
      return;
     }

   num          = FuncTab[i].NumArg;
   testfun_libm = (fptype (*)(fptype)) FuncTab[i].libm;
   testfun_mpfr = (int (*)(mpfr_ptr, mpfr_srcptr, mp_rnd_t)) FuncTab[i].mpfr;
   range_min    = FuncTab[i].min;
   range_max    = FuncTab[i].max;
   monoton      = FuncTab[i].monoton;
   symm         = FuncTab[i].symm;


   if (verbose >= 2)
     {
       if (num == 1)
	 printf ("Testing function %s for exponent %ld.\n", func, e1);
       else
	 printf ("Testing function %s for exponents %ld and %ld.\n", func, e1, e2);
     }
  
  mpfr_init2 (op1, FPPREC);
  mpfr_init2 (op2, FPPREC);
  mpfr_init2 (z, FPPREC);
  mpfr_init2 (rmpfr, FPPREC);
  mpfr_init2 (rmpfrref, 2*FPPREC);

  max_err_near = 0.0;
  max_err_dir = 0.0;

  for (rnd=0 ; rnd < MAX_RND ; rnd++)
    {
      LOG (3, printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd)));
      set_rnd_mode (rnd);
      
      /* reset the seed to test the same sequence of numbers with each
         rounding mode */
      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

      tot = 0;
      umax = 0.0;        /* (signed) maximal error in ulps */
      umax_dir = 0.0;    /* maximal error in ulps when wrong directed rounding */
      tmax = 0.0;
      tmax_dir = 0.0;
      wrong = 0;         /* number of wrong directed roundings */
      wrong_range = 0;   /* number of wrong results wrt range        */
      wrong_monoton = 0; /* number of wrong results wrt monotonicity */
      wrong_symm = 0;    /* number of wrong results wrt symmetry     */

      for (i = 0; i < N ; i++)
	{
	  /* Get a random number */
 	  mpfr_urandomb (op1, state);
	  MPFR_EXP (op1) = e1;
	  if (num == 2)
	    {
	      mpfr_urandomb (op2, state);
	      MPFR_EXP (op2) = e2;
	    }
	  op1d = (*get_fp) (op1, GMP_RNDN);
	  op2d = (*get_fp) (op2, GMP_RNDN);

	  /*mpfr_dump (op1);*/
	  /* Compute the MPFR / LIBM result */
	  if (num == 2)
	    {
	      (*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		 testfun_mpfr)) (rmpfr, op1, op2, rnd);
	      rlibmd = (*((fptype(*)(fptype, fptype))testfun_libm))
		(op1d, op2d);
	    }
	  else
	    {
	      (*testfun_mpfr) (rmpfr, op1, rnd);
	      rlibmd = (*testfun_libm) (op1d);
	    }
	  rmpfrd = (*get_fp) (rmpfr, rnd);
	  /* Check for correct rounding */
	  if (rmpfrd != rlibmd)
	    {
	      /*printf("recompute!\n");*/
	      /* Recompute MPFR result with more prec */
	      if (num == 2)
		(*((int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t))
		   testfun_mpfr)) (rmpfrref, op1, op2, rnd);
	      else
		(*testfun_mpfr) (rmpfrref, op1, rnd);

	      u = ulp (rlibmd, rmpfrref);
	      if (rnd != GMP_RNDN)
		{
		  mp_rnd_t rnd2;
		  
		  rnd2 = (rnd==GMP_RNDZ) ? (mpfr_sgn(rmpfr)>0 ? 
					    GMP_RNDD : GMP_RNDU): rnd;
		  
		  if ((rnd2 == GMP_RNDU && u < 0.0) ||
		      (rnd2 == GMP_RNDD && u > 0.0))
		    {
		      wrong++;
		      if (fabs(u) > fabs(umax_dir))
			{
			  umax_dir = u;
			  xmax_dir = op1d;
			  tmax_dir = op2d;
			}
		    }
		}
	      tot ++;
	      if (fabs(u) > fabs(umax))
		{
		  umax = u;
		  xmax = op1d;
		  tmax = op2d;
		}
	    }

	  /* Testing if the range is verified */
	  if (test_range && 
	      ((range_min/2.0 != range_min) || (range_max/2.0 != range_max) ))
	    {
	      /* one of the extremity is not infinite */
	      if ( (rlibmd < range_min) || (rlibmd > range_max))
		{
		  if (wrong_range == 0 && verbose >= 3)
		    {
		      printf ("      outside range for x=");
		      print_fp (op1d);
		      printf ("\n           f(x)=");
		      print_fp (rlibmd);
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
	  if (num == 1 && test_monotonicity && monoton != NO_MONOTON)
	    {
	      fptype xdplus, xdminus, rplus, rminus;

	      mpfr_set (z, op1, GMP_RNDN);
	      mpfr_nextbelow (z);
	      xdminus = (*get_fp) (z, GMP_RNDN);
	      rminus = (*testfun_libm) (xdminus);
	      
	      mpfr_set (z, op1, GMP_RNDN);
	      mpfr_nextabove (z);
	      xdplus = (*get_fp) (z, GMP_RNDN);
	      rplus = (*testfun_libm) (xdplus);
	      
	      if (monoton == INCREASING)
		{
		  if (rminus > rlibmd)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  printf ("      monotonicity not respected for x=");
			  print_fp (op1d);
			  printf ("\n            f(x-)=");
			  print_fp (rminus);
			  printf ("\n      not <= f(x)=");
			  print_fp (rlibmd);
			  printf (" \n");
			}
		      wrong_monoton++;
		    }
		  if (rplus < rlibmd)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  printf ("      monotonicity not respected for x=");
			  print_fp (op1d);
			  printf ("\n              f(x)=");
			  print_fp (rlibmd);
			  printf ("\n      not <= f(x+)=");
			  print_fp (rplus);
			  printf (" \n");
			}
		      wrong_monoton++;
		    }
		}
	      else 
		{/* monoton == DECREASING */
		  if (rminus < rlibmd)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  printf ("      monotonicity not respected for x=");
			  print_fp (op1d);
			  printf ("\n            f(x-)=");
			  print_fp (rminus);
			  printf ("\n      not >= f(x)=");
			  print_fp (rlibmd);
			  printf (" \n");
			}
		      wrong_monoton++;
		    }
		  if (rplus > rlibmd)
		    {
		      if (wrong_monoton == 0 && verbose >= 3)
			{
			  printf ("      monotonicity not respected for x=");
			  print_fp (op1d);
			  printf ("\n              f(x)=");
			  print_fp (rlibmd);
			  printf ("\n      not >= f(x+)=");
			  print_fp (rplus);
			  printf (" \n");
			}
		      wrong_monoton++;
		    }
		}
	    }
	  
	  /* Testing if the symmetry is verified */
	  if (num == 1 && test_symmetry && (rnd==GMP_RNDN) && (symm != NO_SYMM) )
	    {
	      fptype op1dopp, ropp;
	      op1dopp = -op1d;
	      ropp = (*testfun_libm) (op1dopp);
	      if (symm == ODD)
		{
		  if (ropp != -rlibmd)
		    {
		      if (wrong_symm == 0)
			{
			  printf ("      symmetry not respected for x=");
			  print_fp (op1d);
			  printf ("\n           f(x)= ");
			  print_fp (rlibmd);
			  printf ("\n      and f(-x)=");
			  print_fp (ropp);
			  printf (" \n");
			}
		      wrong_symm++;
		    }
		}
	      else
		{ /* symm == EVEN */
		  if (rlibmd != ropp) 
		    {
		      if (wrong_symm == 0) 
			{
			  printf ("      symmetry not respected for x=");
			  print_fp (op1d);
			  printf ("\n           f(x)=");
			  print_fp (rlibmd);
			  printf ("\n      and f(-x)=");
			  print_fp (ropp);
			  printf (" \n");
			}
		      wrong_symm++;
		    }
		}
	    }
	} /* for i */
      
      set_rnd_mode (GMP_RNDN);

      /* if any difference occured, prints the maximal ulp error */
      if (umax != 0.0 && verbose >= 3)
	{
	  printf ("      %f ulp(s) for x=", umax);
	  print_fp (xmax);
	  if (num == 2)
	    {      
	      printf (" t=");
	      print_fp (tmax);
	    }
	  printf ("\n");
	}
      
      if (wrong != 0 && verbose >= 3)
	{
	  printf ("      wrong directed rounding for x=");
	  print_fp (xmax_dir);
	  if (num == 2)
	    {
	      printf (" t=");
	      print_fp (tmax_dir);
	    }
	  printf (" [%f]\n", umax_dir);
	}
      
      umax = fabs(umax);
      
      if (verbose >= 3 && num == 1)
	{
	  if (rnd == GMP_RNDN)
	    printf ("   nb errors range/monotonicity/symmetry: %lu/%lu/%lu\n", 
		    wrong_range, wrong_monoton, wrong_symm);
	  else
	    printf ("   nb errors range/monotonicity: %lu/%lu\n", 
		    wrong_range, wrong_monoton);
	}
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
    } /* rnd */
  
  printf ("Max. errors for %s [exp. %ld]: %f (nearest), %f (directed)\n", 
	  func, e1, max_err_near, max_err_dir);

  if (verbose >= 3)
    printf ("\n");
  
  mpfr_clear (op1);
  mpfr_clear (op2);
  mpfr_clear (rmpfr);
  mpfr_clear (rmpfrref);
  mpfr_clear (z);

  if (max_err_near > MAX_ERR_NEAR)
    MAX_ERR_NEAR = max_err_near;
  
  if (max_err_dir > MAX_ERR_DIR)
    MAX_ERR_DIR = max_err_dir;
}

