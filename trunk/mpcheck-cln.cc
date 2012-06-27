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

#if defined(HAVE_CLN)

#include <cstring>
#include <string>
#include <iostream>
#include "mpcheck.h"

#include "cln/cln.h"

using namespace std;
using namespace cln;

#define MAX_SIZE 1000000

static cl_I integer;
static mpz_t z;
static char tampon[MAX_SIZE];

extern "C" void *
new_fp (mp_prec_t p)
{
  cl_F *a;
  /* note: float_format_t is an undocumented CLN feature that sets
     the bit precision (and not the digit precision as float_format) */
  cl_F b = cl_float (0.0, float_format_t (p));
  a = new cl_F (b);
  return (void *) a;
}

extern "C" void 
del_fp (void *fp)
{
  cl_F *a = (cl_F *) fp;
  delete a;
}

extern "C" void 
set_fp (mpfr_ptr dest, const void *fp)
{
  const cl_F *a = (const cl_F*) fp;

  if (*a == 0.0)
    mpfr_set_ui (dest, 0, GMP_RNDN);
  else {
    cl_idecoded_float idf = integer_decode_float (*a);
    /* idf.mantissa, idf.exponent and idf.sign */
    /* Undocumented function BUT this is what I need... */
    char *str = print_integer_to_string (16, idf.mantissa);
    //printf ("Str=%s\nPrec=%lu\n", str, mpfr_get_prec (dest));
    /* Warning: When converting to integer*2^N, the integer may be out of
       range! Fix it by temporary augment the exponent range. */
    mp_exp_t emax = mpfr_get_emax ();
    mpfr_set_emax (MPFR_EMAX_DEFAULT);
    int i = mpfr_strtofr (dest, str, NULL, 16, GMP_RNDN);
    /* if (i != 0) { // Sometimes we get more bit than needed...
      printf ("strtofr failed. Must be exact!\n");
      exit (1);
      } */
    mpfr_mul_2si (dest, dest, cl_I_to_long (idf.exponent), GMP_RNDN);
    if (idf.sign < 0)
      mpfr_neg (dest, dest, GMP_RNDN);
    mpfr_set_emax (emax);
    mpfr_check_range (dest, i, GMP_RNDN);
    free_hook (str);
  }
}

extern "C" void 
get_fp (void *fp, mpfr_srcptr src)
{
  cl_F *a = (cl_F*) fp;
  if (mpfr_zero_p (src))
    *a = 0.0;
  else {
    long e = mpfr_get_z_exp (z, src);
    if (MAX_SIZE < mpfr_get_prec (src)) {
      printf ("Abort: Too big prec\n");
      abort ();
    }
    mpz_get_str (tampon, 10, z);
    strcat (tampon, ".0"); /* Convertion must finish with .0 :( */
    *a = tampon;           /* should be exact */
    *a = scale_float (*a, e);
  }
}

extern "C" int 
set_rnd_mode (mp_rnd_t rnd) {
  return rnd == GMP_RNDN;
}

/* Code to check */
extern "C" void 
my_add (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = *(const cl_F*)src1 + *(const cl_F*)src2;
}

extern "C" void 
my_sub (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = *(const cl_F*)src1 - *(const cl_F*)src2;
}

extern "C" void 
my_mul (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = *(const cl_F*)src1 * *(const cl_F*)src2;
}

extern "C" void 
my_div (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = *(const cl_F*)src1 / *(const cl_F*)src2;
}

extern "C" void 
my_sqrt (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = sqrt (*(const cl_F*)src1);
}

extern "C" void 
my_exp (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = exp (*(const cl_F*)src1);
}

extern "C" void 
my_log (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = ln (*(const cl_F*)src1);
}

extern "C" void 
my_sin (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = sin (*(const cl_F*)src1);
}

extern "C" void 
my_cos (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = cos (*(const cl_F*)src1);
}

extern "C" void 
my_tan (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = tan (*(const cl_F*)src1);
}

extern "C" void 
my_asin (void *dest, const void *src1, const void *src2) {
  cl_N tmp = asin (*(const cl_F*)src1);
  *(cl_F*)dest = The(cl_F)(tmp); 
}

extern "C" void 
my_acos (void *dest, const void *src1, const void *src2) {
  cl_N tmp = acos (*(const cl_F*)src1);
  *(cl_F*)dest = The(cl_F)(tmp); 
}

extern "C" void 
my_atan (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = atan (*(const cl_F*)src1);
}


extern "C" void 
my_sinh (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = sinh (*(const cl_F*)src1);
}

extern "C" void 
my_cosh (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = cosh (*(const cl_F*)src1);
}

extern "C" void 
my_tanh (void *dest, const void *src1, const void *src2) {
  *(cl_F*)dest = tanh (*(const cl_F*)src1);
}

extern "C" void 
my_asinh (void *dest, const void *src1, const void *src2) {
  cl_N tmp = asinh (*(const cl_F*)src1);
  *(cl_F*)dest = The(cl_F)(tmp); 
}

extern "C" void 
my_acosh (void *dest, const void *src1, const void *src2) {
  cl_N tmp = acosh (*(const cl_F*)src1);
  *(cl_F*)dest = The(cl_F)(tmp); 
}

extern "C" void 
my_atanh (void *dest, const void *src1, const void *src2) {
  cl_N tmp = atanh (*(const cl_F*)src1);
  *(cl_F*)dest = The(cl_F)(tmp);
}

extern "C" void 
my_pow (void *dest, const void *src1, const void *src2) {
  cl_N tmp = expt (*(const cl_F*)src1, *(const cl_F*)src2);
  *(cl_F*)dest = The(cl_F)(tmp);
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
  {"atan", my_atan, 0, 0},
  {"atan", my_atan, 53, 0},
  {"asin", my_asin, 0, 0},
  {"asin", my_asin, -10, 0},
  //{"acos", my_acos, 0, 0}, /* Internal error: statement in file float/elem/cl_F_compare.cc, line 22 has been reached!! */
  //{"acos", my_acos, -10, 0},
  {"sinh", my_sinh, 0, 0},
  {"sinh", my_sinh, 9, 0}, /* TODO: Improve overflow detection */
  {"cosh", my_cosh, 0, 0},
  {"cosh", my_cosh, 9, 0}, /* TODO: Improve overflow detection */
  {"tanh", my_tanh, 0, 0},
  {"tanh", my_tanh, 4, 0}, 
  {"asinh", my_asinh, 0, 0},
  {"asinh", my_asinh, 9, 0}, /* TODO */
  //{"acosh", my_acosh, 1, 0}, * Internal error: statement in file float/elem/cl_F_compare.cc, line 22 has been reached!! */
  //{"acosh", my_acosh, 9, 0},
  {"atanh", my_atanh, 0, 0},
  {"atanh", my_atanh,-10, 0},
  {"pow", my_pow, 0, 0},
  {"pow", my_pow, 5, 4},
  {"pow", my_pow, 16, 10},

  {NULL, NULL, 0, 0}
};

int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;
  gmp_randstate_t rand;

  mpz_init (z);
  gmp_randinit_default(rand);

  /* Check if interface works */
  mpfr_inits2 (1130, x, y, NULL);
  fp = new_fp (1130);
  mpfr_urandom (x, rand, MPFR_RNDN);
  get_fp (fp, x);
  set_fp (y, fp);
  if (mpfr_cmp (x, y) != 0 || mpfr_nan_p (y)) {
    printf ("set_fp(get_fp)) != Identity\n");
    exit (1);
  }
  del_fp (fp);

  /* Starting MPCHECK */
  mpcheck_init (argc, argv, 53, -1L<<10, 1L<<10,
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
  fprintf (stderr, "CLN not available.\n");
  return 0;
}

#endif
