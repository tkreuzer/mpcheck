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
static int bf_status;

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
  if (s > len)
      s = len;
  dest->_mpfr_exp = r->expn;
  dest->_mpfr_sign = r->sign == 1 ? -1 : 1;
  /* BF floats can be shorter */
  memset(dest->_mpfr_d, 0, (len - s) * sizeof(limb_t));
  memcpy(dest->_mpfr_d + (len - s), r->tab, s * sizeof (limb_t));
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
  size_t len = 1 + (src->_mpfr_prec - 1) / LIMB_BITS;
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
    bf_flags_t v;
    
    switch (rnd) {
    case MPFR_RNDN: v = BF_RNDN; break;
    case MPFR_RNDZ: v = BF_RNDZ; break;
    case MPFR_RNDU: v = BF_RNDU; break;
    case MPFR_RNDD: v = BF_RNDD; break;
    case MPFR_RNDA: v = BF_RNDA; break;
    case MPFR_RNDF: v = BF_RNDF; break;
    default:
        return 0;
    }
    bf_rnd = (bf_rnd & ~BF_RND_MASK) | v;
  return 1;
}

static void bf_feclearexcept(void)
{
    bf_status = 0;
}

static int bf_fetestexcept(int flag)
{
    switch(flag) {
    case FE_INEXACT:
        return (bf_status & BF_ST_INEXACT) != 0;
    case FE_OVERFLOW:
        return (bf_status & BF_ST_OVERFLOW) != 0;
    case FE_UNDERFLOW:
        return (bf_status & BF_ST_UNDERFLOW) != 0;
    default:
        abort();
    }
}

/* Start Libbf Function to check */
void my_libbf_add (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_add ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_sub (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_sub ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_mul (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_mul ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_div (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_div ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}
void my_libbf_fmod (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_rem ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd, BF_RNDZ);
}
void my_libbf_remainder (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_rem ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd, BF_RNDN);
}
void my_libbf_pow (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_pow ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
}

void my_libbf_sqrt (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_sqrt ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_exp (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_exp ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_log (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_log ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_cos (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_cos ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_sin (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_sin ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_tan (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_tan ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_acos (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_acos ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_asin (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_asin ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_atan (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src2;
  (void) src3;
  bf_status |= bf_atan ((bf_t*) dest, (const bf_t*) src1, bf_prec, bf_rnd);
}
void my_libbf_atan2 (void *dest, const void *src1, const void *src2, const void *src3) {
  (void) src3;
  bf_status |= bf_atan2 ((bf_t*) dest, (const bf_t*) src1, (const bf_t*) src2, bf_prec, bf_rnd);
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
  {"atan2", my_libbf_atan2, 0, 0},
  {"atan2", my_libbf_atan2, 53, 0},
  {"asin", my_libbf_asin, 0, 0},
  {"asin", my_libbf_asin, -10, 0},
  {"acos", my_libbf_acos, 0, 0},
  {"acos", my_libbf_acos, -10, 0},
  {"pow", my_libbf_pow, 0, 0},
  {"pow", my_libbf_pow, 5, 4},
  {"pow", my_libbf_pow, 16, 10},
  {"fmod", my_libbf_fmod, 10, 0},
  {"remainder", my_libbf_remainder, 10, 0},
  {NULL, NULL, 0, 0}
};

int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;
  gmp_randstate_t state;
  int exp_bits, emin, emax, prec;
  int bf_errno_check = 0; /* libbf never sets errno */
  
  gmp_randinit_default (state);
  assert(GMP_NUMB_BITS == LIMB_BITS);

  bf_context_init(&context, my_bf_realloc, NULL);

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

  prec = 64;
  exp_bits = 16;
  assert(exp_bits >= BF_EXP_BITS_MIN && exp_bits <= BF_EXP_BITS_MAX);
  bf_rnd = bf_set_exp_bits(exp_bits) | BF_FLAG_SUBNORMAL;
  /* IEEE 754 convention for emin/emax */
  emin = -((slimb_t)1 << exp_bits) / 2 + 3;
  emax = ((slimb_t)1 << exp_bits) / 2;
  
  mpcheck_set_exception_functions(bf_feclearexcept, bf_fetestexcept,
                                  bf_errno_check);

  mpcheck_init (argc, argv, prec, emin, emax,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		ALL_RND, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
  gmp_randclear (state);
  bf_context_end (&context);
  bf_clear_cache (&context);
  return 0;
}


#else

#include <stdio.h>
int main () {
  fprintf (stderr, "LIBBF not available.\n");
  return 0;
}

#endif
