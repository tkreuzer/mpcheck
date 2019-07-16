/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2005-2019 INRIA
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

#if defined(HAVE_LIBBF)

#include "mpcheck.h"

#include "libbf.h"

static long bf_prec;
static bf_context_t context;
static bf_flags_t bf_rnd;
 
static void *my_bf_realloc(void *opaque, void *ptr, size_t size)
{
    return realloc(ptr, size);
}

static void *new_fp (mp_prec_t p)
{
  if (bf_prec == 0) {
    bf_prec = p;
  } else if (p != bf_prec) {
    fprintf (stderr, "ERROR: Requested prec: %lu. FPTYPE prec: %lu\n",
             p, bf_prec);
    abort ();
  }

  bf_t *r = malloc(sizeof(bf_t));
  if (!p) abort();

  bf_context_init(&context, my_bf_realloc, NULL);
  bf_init(&context, r);
  
  return (void *) r;
}

static void del_fp (void *fp)
{
  bf_delete((bf_t*)fp);
  free(fp);
}

static void set_fp (mpfr_ptr dest, const void *fp)
{
  const bf_t *r = fp;

  if (r->len == 0) {
    if (r->expn == BF_EXP_ZERO)
      mpfr_set_zero(dest, -r->sign);
    else if (r->expn == BF_EXP_INF)
      mpfr_set_inf(dest, -r->sign);
    else
      mpfr_set_nan(dest);
    return;
  }
  size_t s = r->len;
  size_t len = 1 + (dest->_mpfr_prec -1) / LIMB_BITS;
  len = len < s  ? len : s;
  dest->_mpfr_exp = r->expn;
  dest->_mpfr_sign = r->sign == 1 ? -1 : 1;
  memcpy(dest->_mpfr_d, r->tab, len * sizeof (limb_t) );
  /*
  mpfr_printf("RP=%Re and ", dest);
  char *buf;
  bf_ftoa(&buf, r, 10, 12, BF_FTOA_FORMAT_FREE);
  printf("BF=%s\n", buf);
  free(buf);
  */
}

static void get_fp (void *fp, mpfr_srcptr src)
{
  bf_t *r = fp;
  size_t len = 1 + (src->_mpfr_prec -1) / LIMB_BITS;
  if (mpfr_nan_p(src)) {
    bf_set_nan(r);
  } else if (mpfr_inf_p(src)) {
    bf_set_inf(r, MPFR_SIGN(src) < 0);
  } else if (mpfr_zero_p(src)) {
    bf_set_zero(r, MPFR_SIGN(src) < 0);
  } else {
    r->sign = MPFR_SIGN(src) < 0;
    r->expn = src->_mpfr_exp;
    bf_resize(r, len);
    memcpy(r->tab, src->_mpfr_d, len * sizeof (limb_t) );    
  }
  /*
  mpfr_printf("RP=%Re and ", src);
  char *buf;
  bf_ftoa(&buf, r, 10, 12, BF_FTOA_FORMAT_FREE);
  printf("BF=%s\n", buf);
  free(buf);
  */
}

static int set_rnd_mode (mpfr_rnd_t rnd) {
  switch (rnd) {
  case MPFR_RNDN: bf_rnd = BF_RNDN; break;
  case MPFR_RNDZ: bf_rnd = BF_RNDZ; break;
  case MPFR_RNDA: bf_rnd = BF_RNDU; break;
  case MPFR_RNDD: bf_rnd = BF_RNDD; break;
  case MPFR_RNDF: bf_rnd = BF_RNDF; break;
  default: return 0;
  }
  return 1;
}

/* Start Libbf Function to check */
void my_libbf_add (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_add ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_sub (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_sub ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_mul (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_mul ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_div (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_add ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_pow (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_pow ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}

void my_libbf_sqrt (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_sqrt ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_exp (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_exp ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_log (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_log ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_cos (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_cos ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_sin (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_sin ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_tan (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_tan ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_acos (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_acos ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_asin (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_asin ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_atan (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_atan ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}

static mpcheck_user_func_t tab[] = {
  {"add", my_libbf_add, 0, 0},
  {"add", my_libbf_add, LONG_MAX, LONG_MAX},
  {"add", my_libbf_add, LONG_MAX, 0},
  {"sub", my_libbf_sub, LONG_MAX, LONG_MAX},
  {"sub", my_libbf_sub, LONG_MAX, 0},
  {"sub", my_libbf_sub, 0, 0},
  {"mul", my_libbf_mul, 0, 0},
  {"mul", my_libbf_mul, LONG_MAX-1, LONG_MAX-1},
  {"div", my_libbf_div, 0, 0},
  {"div", my_libbf_div, LONG_MAX, LONG_MAX},
  {"sqrt", my_libbf_sqrt, 0, 0},
  {"sqrt", my_libbf_sqrt, LONG_MAX, 0},
  {"sqrt", my_libbf_sqrt, LONG_MAX-2, 0},
  {"exp", my_libbf_exp, 0, 0},
  {"exp", my_libbf_exp, 9, 0},
  {"log", my_libbf_log, 0, 0},
  {"log", my_libbf_log, LONG_MAX, 0},
  {"sin", my_libbf_sin, 0, 0},
  {"sin", my_libbf_sin, 10, 0},
  {"cos", my_libbf_cos, 0, 0},
  {"cos", my_libbf_cos, 10, 0},
  {"tan", my_libbf_tan, 0, 0},
  {"tan", my_libbf_tan, 10, 0},
  {"atan", my_libbf_atan, 0, 0},
  {"atan", my_libbf_atan, 53, 0},
  {"asin", my_libbf_asin, 0, 0},
  {"asin", my_libbf_asin, -10, 0},
  {"acos", my_libbf_acos, 0, 0},
  {"acos", my_libbf_acos, -10, 0},
  {"pow", my_libbf_pow, 0, 0},
  {"pow", my_libbf_pow, 5, 4},
  {"pow", my_libbf_pow, 16, 10},
  {NULL, NULL, 0, 0}
};


int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;
  gmp_randstate_t state;

  gmp_randinit_default (state);
  assert(GMP_NUMB_BITS == LIMB_BITS);

  /* Check if interface works */
  mpfr_inits2 (113, x, y, NULL);
  fp = new_fp (113);
  mpfr_urandomb (x, state);
  mpfr_set (y, x, GMP_RNDN);
  get_fp (fp, x);
  set_fp (y, fp);
  if (mpfr_cmp (x, y) != 0) {
    printf ("ERROR: set_fp(get_fp)) != Identity\n");
    exit (1);
  }
  mpfr_clears (x, y, NULL);
  del_fp (fp);
  bf_prec = 0;

  mpcheck_init (argc, argv, 64, -1L<<16, 1L<<16,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		ALL_RND, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
  gmp_randclear (state);
  return 0;
}


#else

#include <stdio.h>
int main () {
  fprintf (stderr, "LIBBF not available.\n");
  return 0;
}

#endif
