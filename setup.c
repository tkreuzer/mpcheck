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

/* Check for Linux only? */
#if ((defined (__i386__) || defined (__i486__)))
# define _FPU_EXTENDED 0x0300
# define _FPU_DOUBLE   0x0200
# define _FPU_DEFAULT  0x137f
# define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
# define _fpu_ieee ((_FPU_DEFAULT & (~_FPU_EXTENDED)) | _FPU_DOUBLE)
# define set_double() __setfpucw(_fpu_ieee)
# define set_float()  __setfpucw(_FPU_DEFAULT & (~_FPU_EXTENDED))
#else
# define set_double()
# define set_float()
#endif

void
fprint_d (FILE *stream, fptype x)
{
  fprintf (stream, "%1.20e", x);
}

void
fprint_ld (FILE *stream, fptype x)
{
#if (FPPREC <= 64)
  fprintf (stream, "%1.24Le", (long double) x);
#else
  fprintf (stream, "%1.38Le", x);
#endif
}

int
mpfr_set_float (mpfr_ptr x, float y, mp_rnd_t rnd)
{
  return mpfr_set_d (x, (double) y, rnd);
}

float
mpfr_get_float (mpfr_srcptr x, mp_rnd_t rnd)
{
  return (float) mpfr_get_d (x, rnd);
}

void
setup ()
{
  fptype x, y;
  int j;

  /* sets the functions set_fp and get_fp */
#if (FPPREC == 24)
  set_fp = mpfr_set_float;
  get_fp = mpfr_get_float;
  fprint_fp = fprint_d;
#elif (FPPREC == 53)
  set_fp = mpfr_set_d;
  get_fp = mpfr_get_d;
  fprint_fp = fprint_d;
#else
  set_fp = mpfr_set_ld;
  get_fp = mpfr_get_ld;
  fprint_fp = fprint_ld;
#endif

  /* Special checking for double numbers */
#if (FPPREC == 53)
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
#elif (FPPREC == 24)
 set_float ();
#endif

  /* Checks that setting machine rounding mode works */
  set_rnd_mode (GMP_RNDD);
  x = 2.0; /* don't write x = sqrt (2.0) and y = sqrt (2.0) otherwise the
              compiler may optimize too much */
  x = sqrt (x);
  set_rnd_mode (GMP_RNDU);
  y = 2.0;
  y = sqrt (y);
  if (x == y)
    fprintf (stderr, "WARNING: setting rounding modes doesn't seem to affect libm.\n");

  /* Checks that the precision of fptype is FPPREC */
  set_rnd_mode (GMP_RNDN);
  x = 1.0;
  for (j=0; x + 1.0 != x; j++, x = 2.0 * x);
  if (j != FPPREC)
    {
      fprintf (stderr, "Precision of fptype is not %u but %u\n", 
	       FPPREC, j);
      exit (1);
    }

  /* Checks that the exponents EMIN and EMAX are correct */
  x = 0.5; /* EXP(x)=0 */
  for (j=0; x + x != x; j++, x = 2.0 * x);
  if (j-1 != EMAX)
    {
      fprintf (stderr, "Maximal exponent is %d instead of %d\n", j-1, EMAX);
      exit (1);
    }
  for (x=1.0, j=0; j<FPPREC; j++, x = x / 2.0);
  x = 0.5 + x;
  /* x = 0.5 + 2^(-FPPREC), EXP(x)=0 */
  for (j=0; ; j--)
    {
      y = x;
      x = x / 2.0;
      if (2.0 * x != y)
        break;
    }
  if (j != EMIN)
    {
      fprintf (stderr, "Minimal exponent is %d instead of %d\n", j, EMIN);
      exit (1);
    }

  /* sets the minimum and maximum exponents */
  mpfr_set_emax (EMAX);
  mpfr_set_emin (EMIN);

}
