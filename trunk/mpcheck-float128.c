/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2012 INRIA
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

#if defined(HAVE_FLOAT128)

#define MPFR_WANT_FLOAT128 1
#include "quadmath.h"
#include "mpcheck.h"

#define fptype __float128
#define NAME(x) x ## q

#undef HAVE_EXP
#define HAVE_EXP HAVE_EXPL
#undef HAVE_EXP2
#undef HAVE_EXP10
#define HAVE_EXP10 HAVE_EXP10L
#undef HAVE_EXPM1
#define HAVE_EXPM1 HAVE_EXPM1L
#undef HAVE_SIN
#define HAVE_SIN HAVE_SINL
#undef HAVE_COS
#define HAVE_COS HAVE_COSL
#undef HAVE_TAN
#define HAVE_TAN HAVE_TANL
#undef HAVE_ASIN
#define HAVE_ASIN HAVE_ASINL
#undef HAVE_ACOS
#define HAVE_ACOS HAVE_ACOSL
#undef HAVE_ATAN
#define HAVE_ATAN HAVE_ATANL
#undef HAVE_ATAN2
#define HAVE_ATAN2 HAVE_ATAN2L
#undef HAVE_SINH
#define HAVE_SINH HAVE_SINHL
#undef HAVE_COSH
#define HAVE_COSH HAVE_COSHL
#undef HAVE_TANH
#define HAVE_TANH HAVE_TANHL
#undef HAVE_ASINH
#define HAVE_ASINH HAVE_ASINHL
#undef HAVE_ACOSH
#define HAVE_ACOSH HAVE_ACOSHL
#undef HAVE_ATANH
#define HAVE_ATANH HAVE_ATANHL
#undef HAVE_LOG
#define HAVE_LOG HAVE_LOGL
#undef HAVE_LOG2
#define HAVE_LOG2 HAVE_LOG2L
#undef HAVE_LOG10
#define HAVE_LOG10 HAVE_LOG10L
#undef HAVE_LOG1P
#define HAVE_LOG1P HAVE_LOG1PL
#undef HAVE_TGAMMA
#define HAVE_TGAMMA HAVE_TGAMMAL
#undef HAVE_SQRT
#define HAVE_SQRT HAVE_SQRTL
#undef HAVE_POW
#define HAVE_POW HAVE_POWL
#undef HAVE_CBRT
#define HAVE_CBRT HAVE_CBRTL
#undef HAVE_ERF
#define HAVE_ERF HAVE_ERFL
#undef HAVE_ERFC
#define HAVE_ERFC HAVE_ERFCL

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mpfr_set_float128 (dest, *(__float128 *)fp, GMP_RNDN);
}
static void get_fp (void *fp, mpfr_srcptr src)
{
  *(__float128 *) fp = mpfr_get_float128 (src, GMP_RNDN);
}

#include "native.c"

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "__float128 not available.\n");
  return 0;
}

#endif
