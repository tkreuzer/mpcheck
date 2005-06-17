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

#if defined(HAVE_NTL)

#include "mpcheck.h"

#define NTL_STD_CXX
#include "NTL/RR.h"

using namespace std;
using namespace NTL;

#define MAX_PREC 8000000

static mp_prec_t prec = 0;
static unsigned long prec_ul;
static int count = 0;
static unsigned long stck;
static mpz_t z;
static ZZ zz;
static unsigned char buffer[MAX_PREC/8]; 
static long buffer_size;

extern "C" void *
new_fp (mp_prec_t p)
{
  RR *a;

  if (p > MAX_PREC)
    abort ();
  a = new RR;
  a->SetPrecision (p);
  return (void *) a;
}

extern "C" void 
del_fp (void *fp)
{
  RR *b = (RR *) fp;
  delete b;
}

extern "C" void 
set_fp (mpfr_ptr dest, const void *fp)
{
  RR *a = (RR*) fp;

  if (*a == 0)
    mpfr_set_ui (dest, 0, GMP_RNDN);
  else {
    long size = NumBytes (a->mantissa());
    BytesFromZZ (buffer, a->mantissa (), size);
    mpz_import (z, size, -1, 1, 0, 0, buffer);
    mpfr_set_z (dest, z, GMP_RNDN);
    mpfr_mul_2si (dest, dest, a->exponent(), GMP_RNDN);
    if (*a < 0)
      mpfr_neg (dest, dest, GMP_RNDN);
  }
}

extern "C" void 
get_fp (void *fp, mpfr_srcptr src)
{
  RR *a = (RR*) fp;
  if (mpfr_zero_p (src))
    *a = 0.0;
  else {
    long e;
    size_t size;

    e = mpfr_get_z_exp (z, src); /* src = z*2^e */
    mpz_export (buffer, &size, -1, 1, 0, 0, z);
    ZZFromBytes (zz, buffer, size);
    if (mpz_sgn (z) < 0)
      zz = -zz;
    *a = to_RR (zz);
    *a = *a * power2_RR (e);
  }
}

extern "C" int 
set_rnd_mode (mp_rnd_t rnd) {
  return rnd == GMP_RNDN;
}

/* Code to check */
extern "C" void 
my_add (void *dest, const void *src1, const void *src2) {
  *(RR*)dest = *(const RR*)src1 + *(const RR*)src2;
}

extern "C" void 
my_sub (void *dest, const void *src1, const void *src2) {
  *(RR*)dest = *(const RR*)src1 - *(const RR*)src2;
}

extern "C" void 
my_mul (void *dest, const void *src1, const void *src2) {
  *(RR*)dest = *(const RR*)src1 * *(const RR*)src2;
}

extern "C" void 
my_div (void *dest, const void *src1, const void *src2) {
  *(RR*)dest = *(const RR*)src1 / *(const RR*)src2;
}

static mpcheck_user_func_t tab[] = {
  {"add", my_add, 0, 0},
  {"add", my_add, LONG_MAX, LONG_MAX},
  {"sub", my_sub, LONG_MAX, LONG_MAX},
  {"sub", my_sub, 0, 0},
  {"mul", my_mul, 0, 0},
  {"mul", my_mul, LONG_MAX-1, LONG_MAX-1},
  {"div", my_div, 0, 0},
  {"div", my_div, LONG_MAX, LONG_MAX},
#if 0
  {"sqrt", my_sqrt, 0, 0},
  {"sqrt", my_sqrt, LONG_MAX, 0},
  {"sqrt", my_sqrt, LONG_MAX-2, 0},
  {"exp", my_exp, 0, 0},
  {"exp", my_exp, 9, 0}, /* FIXME: Improve by autodection of overflow? */
  {"log", my_log, 0, 0},
  {"log", my_log, LONG_MAX, 0},
  {"sin", my_sin, 0, 0},
  {"sin", my_sin, 10, 0},
  {"cos", my_cos, 0, 0},
  {"cos", my_cos, 10, 0},
  {"tan", my_tan, 0, 0},
  {"tan", my_tan, 10, 0},
  {"atan", my_atan, 0, 0},
  {"atan", my_atan, 53, 0},
  {"asin", my_asin, 0, 0},
  {"asin", my_asin, -10, 0},
  {"acos", my_acos, 0, 0},
  {"acos", my_acos, -10, 0},
  {"gamma", my_gamma, 0, 0},
  {"gamma", my_gamma, 10, 0},
  {"pow", my_pow, 0, 0},
  {"pow", my_pow, 5, 4},
  {"pow", my_pow, 16, 10},
  {"erfc", my_erfc, 0, 0},
  {"erfc", my_erfc, 2, 0},
#endif
  {NULL, NULL, 0, 0}
};

int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;

  mpz_init (z);

  /* Check if interface works */
  mpfr_inits2 (113, x, y, NULL);
  fp = new_fp (113);
  mpfr_random (x);
  mpfr_set (y, x, GMP_RNDN);
  get_fp (fp, x);
  set_fp (y, fp);
  if (mpfr_cmp (x, y) != 0) {
    printf ("set_fp(get_fp)) != Identity\n");
    exit (1);
  }
  del_fp (fp);

  /* Starting MPCHECK */
  mpcheck_init (argc, argv, 53, -1L<<20, 1L<<20,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		1, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);

  mpz_clear (z);
  return 0;
}

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "NTL not available.\n");
  return 0;
}

#endif
