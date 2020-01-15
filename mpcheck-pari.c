/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2005, 2010 INRIA
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

#if defined(HAVE_PARI)

#include "mpcheck.h"

typedef unsigned long ulong;
#include "pari/pari.h"

static mp_prec_t prec = 0;
static unsigned long prec_ul;
static int count = 0;
static unsigned long stck;

static void *new_fp (mp_prec_t p)
{
  GEN a;
  char sprec[20];

  /* Setup internal prec */
  sprintf (sprec, "%lu", p);
  sd_realbitprecision (sprec, d_SILENT);

  if (prec == 0) {
    prec = p;
    prec_ul = (prec-1)/(8*sizeof (mp_limb_t))+1;
  } else if (p != prec) {
    fprintf (stderr, "ERROR: Requested prec: %lu. FPTYPE prec: %lu\n",
             p, prec);
    abort ();
  }

  /* Save Pari stack */
  if (count++ == 0)
    stck = avma;

  /* the 2nd argument (prec) of Pari mathematical functions is the number
     of words used to store the t_REAL type, i.e., 2 plus the number of
     words to store the significand. */
  a = gsqrt (stoi(3), nbits2prec (prec));

  return (void *) a;
}

static void del_fp (void *fp)
{
  if (--count == 0)
    avma = stck;
}

/* The philosophy of Pari is to hack ;)
   Do it. */
static void copy_reverse (mp_limb_t *dest, const mp_limb_t *src, mp_size_t n)
{
  mp_size_t i;
  for (i = 0; i < n ;i++)
    dest[i] = src[n-1-i];
}

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mp_prec_t p;
  long s = gsigne ((GEN) fp);

  if (s == 0)
    mpfr_set_ui (dest, 0, GMP_RNDN);
  else {
    mpfr_set_si (dest, s, GMP_RNDN);
    if (mpfr_set_exp (dest, gexpo ((GEN) fp)+1)) /* underflow/overflow */
      {
        if (gexpo ((GEN) fp)+1 < mpfr_get_emin ())
          mpfr_set_exp (dest, mpfr_get_emin ());
        else
          {
            mpfr_set_exp (dest, mpfr_get_emax ());
            mpfr_mul_2exp (dest, dest, 2, MPFR_RNDZ);
          }
      }
    else
      copy_reverse (dest->_mpfr_d, (mp_limb_t*) & ((long*)fp)[2], prec_ul);
    /* We may have copied more bit than necessary: reprec it */
    p = mpfr_get_prec (dest);
    dest->_mpfr_prec = prec_ul * (8*sizeof(mp_limb_t));
    mpfr_prec_round (dest, p, mpcheck_rnd_mode);
  }
  /* printf ("Set Fp="); output ((GEN) fp); */
}

static void get_fp (void *fp, mpfr_srcptr src)
{
  /*printf ("MPFR :");
    mpfr_out_str (stdout, 10, 0, src, GMP_RNDN); putchar ('\n'); */
  if (mpfr_nan_p (src) || mpfr_inf_p (src))
    {
      printf ("Error: Pari can't handle NAN or INF\n");
      exit (1);
    }
  if (!mpfr_zero_p (src)) {
    setsigne (fp, mpfr_sgn (src));
    setexpo (fp, mpfr_get_exp (src)-1);
    copy_reverse ((mp_limb_t*) & ((long*)fp)[2], src->_mpfr_d, prec_ul);
  } else {
    gsubz (fp, fp, fp); /* Set zero #top#! */
  }
  /*  printf ("Get Fp=");
      output ((GEN) fp); */
}

static int set_rnd_mode (mp_rnd_t rnd) {
  return rnd == GMP_RNDN; /* Pari has no Rounding Mode */
}

/* Start Pari Function to check */
void my_pari_add (void *dest, const void *src1, const void *src2) {
  gaddz ((GEN) src1, (GEN) src2, (GEN) dest);
}
void my_pari_sub (void *dest, const void *src1, const void *src2) {
  gsubz ((GEN) src1, (GEN) src2, (GEN) dest);
}
void my_pari_mul (void *dest, const void *src1, const void *src2) {
  gmulz ((GEN) src1, (GEN) src2, (GEN) dest);
}
void my_pari_div (void *dest, const void *src1, const void *src2) {
  gdivz ((GEN) src1, (GEN) src2, (GEN) dest);
}
void my_pari_sqrt (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gsqrt ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gsqrtz ((GEN) src1, (GEN) dest); */
}
void my_pari_exp (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gexp ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gexpz ((GEN) src1, (GEN) dest); */
}
void my_pari_log (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = glog ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* glogz ((GEN) src1, (GEN) dest); */
}

void my_pari_sin (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gsin ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gsinz ((GEN) src1, (GEN) dest); */
}
void my_pari_sinh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gsinh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_cos (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gcos ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gcosz ((GEN) src1, (GEN) dest); */
}
void my_pari_cot (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gcotan ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_coth (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gcotanh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_cosh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gch ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_tan (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gtan ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gtanz ((GEN) src1, (GEN) dest); */
}
void my_pari_tanh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gtanh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_asin (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gasin ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gasinz ((GEN) src1, (GEN) dest); */
}
void my_pari_asinh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gasinh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_acos (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gacos ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gacosz ((GEN) src1, (GEN) dest); */
}
void my_pari_acosh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gacosh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_atan (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gatan ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gatanz ((GEN) src1, (GEN) dest); */
}
void my_pari_atanh (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gatanh ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}

void my_pari_tgamma (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = ggamma ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* ggammaz ((GEN) src1, (GEN) dest); */
}
void my_pari_lngamma (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = glngamma ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
/* Doesn't work very well... Does a lot of assert failed */
void my_pari_erfc (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gerfc ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gerfcz ((GEN) src1, (GEN) dest); */
}
void my_pari_j0 (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN nu = gen_0;
  GEN tmp;
  tmp = jbessel (nu, (GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_j1 (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN nu = gen_1;
  GEN tmp;
  tmp = jbessel (nu, (GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_dilog (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = dilog ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_eint (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  /* Definition of eint in MPFR doesn't match eint1 in Pari:
     for x >= 0: mpfr_eint (x) = -Re(pari_eint1(-x)),
     for x < 0:  mpfr_eint (x) = pari_eint1(-x). */
  if (signe ((GEN) src1) >= 0)
    {
      tmp = gneg ((GEN) src1);
      tmp = eint1 (tmp, nbits2prec (prec));
      tmp = greal (tmp);
      tmp = gneg (tmp);
    }
  else
    {
      tmp = gneg ((GEN) src1);
      tmp = eint1 (tmp, nbits2prec (prec));
    }
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_expm1 (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gexpm1 ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_zeta (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp;
  tmp = gzeta ((GEN) src1, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_pow (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp = gpow ((GEN) src1, (GEN) src2, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
  /* gpowz ((GEN) src1, (GEN) src2, (GEN) dest); */
}
void my_pari_agm (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp = agm ((GEN) src1, (GEN) src2, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_gamma_inc (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp = incgam ((GEN) src1, (GEN) src2, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}
void my_pari_root (void *dest, const void *src1, const void *src2) {
  unsigned long st = avma;
  GEN tmp = gsqrtn ((GEN) src1, (GEN) src2, NULL, nbits2prec (prec));
  gaffect (tmp, (GEN) dest);
  avma = st;
}

static mpcheck_user_func_t tab[] = {
  {"add", my_pari_add, 0, 0},
  {"add", my_pari_add, LONG_MAX, LONG_MAX},
  {"add", my_pari_add, LONG_MAX, 0},
  {"sub", my_pari_sub, LONG_MAX, LONG_MAX},
  {"sub", my_pari_sub, LONG_MAX, 0},
  {"sub", my_pari_sub, 0, 0},
  {"mul", my_pari_mul, 0, 0},
  {"mul", my_pari_mul, LONG_MAX-1, LONG_MAX-1},
  {"div", my_pari_div, 0, 0},
  {"div", my_pari_div, LONG_MAX, LONG_MAX},
  {"sqrt", my_pari_sqrt, 0, 0},
  {"sqrt", my_pari_sqrt, LONG_MAX, 0},
  {"sqrt", my_pari_sqrt, LONG_MAX-2, 0},
  {"exp", my_pari_exp, 0, 0},
  {"exp", my_pari_exp, 9, 0},
  {"log", my_pari_log, 0, 0},
  {"log", my_pari_log, LONG_MAX, 0},
  {"sin", my_pari_sin, 0, 0},
  {"sin", my_pari_sin, 10, 0},
  {"cos", my_pari_cos, 0, 0},
  {"cos", my_pari_cos, 10, 0},
  {"tan", my_pari_tan, 0, 0},
  {"tan", my_pari_tan, 10, 0},
  {"atan", my_pari_atan, 0, 0},
  {"atan", my_pari_atan, 53, 0},
  {"asin", my_pari_asin, 0, 0},
  {"asin", my_pari_asin, -10, 0},
  {"acos", my_pari_acos, 0, 0},
  {"acos", my_pari_acos, -10, 0},
  {"tgamma", my_pari_tgamma, 0, 0},
  {"tgamma", my_pari_tgamma, 10, 0},
  {"pow", my_pari_pow, 0, 0},
  {"pow", my_pari_pow, 5, 4},
  {"pow", my_pari_pow, 16, 10},
  {"erfc", my_pari_erfc, 0, 0},
  {"erfc", my_pari_erfc, 2, 0},
  {"cosh", my_pari_cosh, 0, 0},
  {"sinh", my_pari_sinh, 0, 0},
  {"tanh", my_pari_tanh, 0, 0},
  {"acosh", my_pari_acosh, 1, 0},
  {"asinh", my_pari_asinh, 0, 0},
  {"atanh", my_pari_atanh, 0, 0},
  {"agm", my_pari_agm, 0, 0},
  {"j0", my_pari_j0, 0, 0},
  {"j1", my_pari_j1, 0, 0},
  {"cot", my_pari_cot, 0, 0},
  {"coth", my_pari_coth, 0, 0},
  {"dilog", my_pari_dilog, 0, 0},
  {"eint", my_pari_eint, 0, 0},
  {"expm1", my_pari_expm1, 0, 0},
  {"gamma_inc", my_pari_gamma_inc, 0, 0},
  {"lngamma", my_pari_lngamma, 0, 0},
  {"zeta", my_pari_zeta, 0, 0},
  {NULL, NULL, 0, 0}
};


int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;
  gmp_randstate_t state;

  /* Check if Pari may work */
  if (sizeof (long) != sizeof (unsigned long)
      || sizeof (void*) != sizeof (long)
      || sizeof (mp_limb_t) != sizeof (long)
      || sizeof (mp_exp_t) != sizeof (long))
    {
      printf ("Type Long is not the base of your system\n");
      exit (1);
    }

  pari_init (40000000, 10000);
  gmp_randinit_default (state);

  /* Check if interface works */
  mpfr_inits2 (113, x, y, NULL);
  fp = new_fp (113);
  mpfr_urandomb (x, state);
  mpfr_set (y, x, GMP_RNDN);
  get_fp (fp, x);
  set_fp (y, fp);
  if (mpfr_cmp (x, y) != 0) {
    printf ("set_fp(get_fp)) != Identity\n");
    exit (1);
  }
  /*printf ("X10="); mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
    printf ("\nFp ="); output ((GEN) fp);*/
  del_fp (fp);
  prec = 0;

#if 0
  /* we can test directly Pari functions like this */
  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str (x, "d.534cb60e42b88@1", 16, MPFR_RNDN);
  mpfr_out_str (stdout, 16, 0, x, MPFR_RNDN);
  printf ("\n");
  fp = new_fp (53);
  get_fp (fp, x);
  printf ("fp ="); output ((GEN) fp);
  GEN tmp = gerfc (fp, nbits2prec (prec));
  printf ("tmp ="); output ((GEN) tmp);
  del_fp (fp);
  prec = 0;
#endif

  mpcheck_init (argc, argv, 53, -1L<<16, 1L<<16,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		ALL_RND, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
  gmp_randclear (state);
  return 0;
}

int
is_accepted (char *name, int numarg, mpfr_t op1, mpfr_t op2)
{
  if (mpfr_nan_p (op1) || mpfr_inf_p (op1))
    return 0;

  if (numarg == 2 && (mpfr_nan_p (op2) || mpfr_inf_p (op2)))
    return 0;

  if (strcmp (name, "div") == 0 && mpfr_zero_p (op2))
    return 0;

  if (strcmp (name, "sqrt") == 0 && mpfr_cmp_ui (op1, 0) < 0)
    return 0;

  if (strcmp (name, "exp") == 0 && (mpfr_cmp_si (op1, 45426) > 0 ||
                                    mpfr_cmp_si (op1, -22749) < 0))
    return 0;

  if (strcmp (name, "expm1") == 0 && (mpfr_cmp_si (op1, 45426) > 0 ||
                                      mpfr_cmp_si (op1, -22749) < 0))
    return 0;

  if (strcmp (name, "log") == 0 && (mpfr_cmp_ui (op1, 0) <= 0 ||
                                    mpfr_get_exp (op1) <= -32820))
    return 0;

  if (strcmp (name, "sin") == 0 && mpfr_get_exp (op1) >= 32768)
    return 0;

  if (strcmp (name, "cos") == 0 && mpfr_get_exp (op1) >= 32768)
    return 0;

  if (strcmp (name, "tan") == 0 && mpfr_get_exp (op1) >= 32768)
    return 0;

  if (strcmp (name, "asin") == 0 && mpfr_get_exp (op1) > 0)
    return 0;

  if (strcmp (name, "acos") == 0 && mpfr_get_exp (op1) > 0)
    return 0;

  if (strcmp (name, "tgamma") == 0 && (mpfr_cmp_ui (op1, 0) <= 0 ||
                                      mpfr_get_exp (op1) >= 32768))
    return 0;
  
  if (strcmp (name, "pow") == 0 && (mpfr_zero_p (op1) && mpfr_zero_p (op2)))
    return 0;

  if (strcmp (name, "erfc") == 0 && mpfr_get_exp (op1) >= 32768)
    return 0;

  if (strcmp (name, "cosh") == 0 && (mpfr_cmp_si (op1, 45426) > 0 ||
                                     mpfr_cmp_si (op1, -22749) < 0))
    return 0;

  if (strcmp (name, "acosh") == 0 && mpfr_cmp_ui (op1, 1) < 0)
    return 0;

  if (strcmp (name, "atanh") == 0 && (mpfr_cmp_si (op1, -1) <= 0 ||
                                      mpfr_cmp_si (op1, 1) >= 0))
    return 0;

  if (strcmp (name, "j0") == 0 && mpfr_get_exp (op1) > 15)
    return 0;

  if (strcmp (name, "j1") == 0 && mpfr_get_exp (op1) > 15)
    return 0;

  if (strcmp (name, "cot") == 0 && (mpfr_cmp_ui (op1, 0) == 0 ||
                                    mpfr_get_exp (op1) > 67))
    return 0;

  if (strcmp (name, "coth") == 0 && mpfr_cmp_ui (op1, 0) == 0)
    return 0;

  if (strcmp (name, "eint") == 0 && (mpfr_get_exp (op1) > 17 ||
                                     mpfr_get_exp (op1) < -63))
    return 0;

  if (strcmp (name, "gamma_inc") == 0 && (mpfr_cmp_ui (op1, 0) == 0 ||
                                          mpfr_cmp_ui (op2, 0) == 0))
    return 0;

  if (strcmp (name, "lngamma") == 0 && mpfr_cmp_ui (op1, 0) <= 0)
    return 0;

  if (strcmp (name, "dilog") == 0 && mpfr_cmp_ui (op1, 1) > 0)
    return 0;

#if 0
  if (strcmp (name, "dilog") == 0)
    mpfr_printf ("op1=%Re\n", op1);
#endif

  return 1;
}

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "PARI not available.\n");
  return 0;
}

#endif
