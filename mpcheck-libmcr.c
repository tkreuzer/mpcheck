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

#if defined(HAVE_LIBMCR)

#include "mpcheck.h"

#include "libmcr.h"

#define LIB_INIT() 
#define LIB_EXIT() 
#define fptype double

#define NAME(x) __libmcr_ ## x

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mpfr_set_d (dest, *(double*)fp, GMP_RNDN);
}
static void get_fp (void *fp, mpfr_srcptr src)
{
  *(double*) fp = mpfr_get_d (src, GMP_RNDN);
}

#undef HAVE_EXP
#define HAVE_EXP 1
#undef HAVE_EXP2
#undef HAVE_EXP10
#undef HAVE_EXPM1
#undef HAVE_SIN
#define HAVE_SIN 1
#undef HAVE_COS
#define HAVE_COS 1
#undef HAVE_TAN
#define HAVE_TAN 1
#undef HAVE_ASIN
#undef HAVE_ACOS
#undef HAVE_ATAN
#define HAVE_ATAN 1
#undef HAVE_SINH
#undef HAVE_COSH
#undef HAVE_TANH
#undef HAVE_ASINH
#undef HAVE_ACOSH
#undef HAVE_ATANH
#undef HAVE_LOG
#define HAVE_LOG 1
#undef HAVE_LOG2
#undef HAVE_LOG10
#undef HAVE_LOG1P
#undef HAVE_TGAMMA
#undef HAVE_SQRT
#undef HAVE_POW
#define HAVE_POW 1
#undef HAVE_CBRT
#undef HAVE_ERF
#undef HAVE_ERFC

#include "native.c"

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "LIBMCR not available.\n");
  return 0;
}

#endif
