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

/* computes the error in ulps between y [computed by libm]
   and z [computed by mpfr]. It assumes y <> z. */
double ulp (fptype y, mpfr_ptr z)
{
  mpfr_t yy;
  double u;
  
  mpfr_init2 (yy, 2*FPPREC);

  (*set_fp) (yy, y, GMP_RNDN); /* Exact */
  mpfr_sub (yy, yy, z, GMP_RNDN);
  MPFR_EXP (yy) += FPPREC - MPFR_EXP (z);
  u = mpfr_get_d (yy, GMP_RNDN);

  mpfr_clear (yy);
  return u;
}

