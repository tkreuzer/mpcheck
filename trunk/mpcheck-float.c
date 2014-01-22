/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2005 INRIA
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

#define fptype float 
#define NAME(x) x ## f

#undef HAVE_EXP
#define HAVE_EXP HAVE_EXPF
#undef HAVE_EXP2
#define HAVE_EXP2 HAVE_EXP2F
#undef HAVE_EXP10
#define HAVE_EXP10 HAVE_EXP10F
#undef HAVE_EXPM1
#define HAVE_EXPM1 HAVE_EXPM1F
#undef HAVE_SIN
#define HAVE_SIN HAVE_SINF
#undef HAVE_COS
#define HAVE_COS HAVE_COSF
#undef HAVE_TAN
#define HAVE_TAN HAVE_TANF
#undef HAVE_ASIN
#define HAVE_ASIN HAVE_ASINF
#undef HAVE_ACOS
#define HAVE_ACOS HAVE_ACOSF
#undef HAVE_ATAN
#define HAVE_ATAN HAVE_ATANF
#undef HAVE_ATAN2
#define HAVE_ATAN2 HAVE_ATAN2F
#undef HAVE_SINH
#define HAVE_SINH HAVE_SINHF
#undef HAVE_COSH
#define HAVE_COSH HAVE_COSHF
#undef HAVE_TANH
#define HAVE_TANH HAVE_TANHF
#undef HAVE_ASINH
#define HAVE_ASINH HAVE_ASINHF
#undef HAVE_ACOSH
#define HAVE_ACOSH HAVE_ACOSHF
#undef HAVE_ATANH
#define HAVE_ATANH HAVE_ATANHF
#undef HAVE_LOG
#define HAVE_LOG HAVE_LOGF
#undef HAVE_LOG2
#define HAVE_LOG2 HAVE_LOG2F
#undef HAVE_LOG10
#define HAVE_LOG10 HAVE_LOG10F
#undef HAVE_LOG1P
#define HAVE_LOG1P HAVE_LOG1PF
#undef HAVE_TGAMMA
#define HAVE_TGAMMA HAVE_TGAMMAF
#undef HAVE_SQRT
#define HAVE_SQRT HAVE_SQRTF
#undef HAVE_POW
#define HAVE_POW HAVE_POWF
#undef HAVE_CBRT
#define HAVE_CBRT HAVE_CBRTF
#undef HAVE_ERF
#define HAVE_ERF HAVE_ERFF
#undef HAVE_ERFC
#define HAVE_ERFC HAVE_ERFCF
#undef HAVE_J0
#define HAVE_J0 HAVE_J0F
#undef HAVE_J1
#define HAVE_J1 HAVE_J1F
#undef HAVE_JN
#define HAVE_JN HAVE_JNF
#undef HAVE_Y0
#define HAVE_Y0 HAVE_Y0F
#undef HAVE_Y1
#define HAVE_Y1 HAVE_Y1F
#undef HAVE_YN
#define HAVE_YN HAVE_YNF

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mpfr_set_d (dest, (double) (*(float*)fp), GMP_RNDN);
}
static void get_fp (void *fp, mpfr_srcptr src)
{
  /* There should be no rounding since src has the prec of a float */
  *(float*) fp = mpfr_get_d (src, GMP_RNDN);
}

#include "native.c"
