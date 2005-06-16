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

#if defined(HAVE_PARI)

#include "mpcheck.h"

typedef unsigned long ulong;
#include "pari/pari.h"

static mp_prec_t prec = 0;
static int count = 0;
static unsigned long stck;

static void *new_fp (mp_prec_t p)
{
  GEN a;

  /* Setup internal prec */
  if (prec == 0)
    prec = p;
  else if (p != prec)
    {
      fprintf (stderr, "ERROR: Requested prec: %lu. FPTYPE prec: %lu\n",
	       p, prec);
      abort ();
    }

  /* Save Pari stack */
  if (count++ == 0)
    stck = avma;

  /* Since PARI can't handle prec bits exactly we request a few more */
  a = gsqrt (stoi(3), (p - 1)/BITS_IN_LONG + 1 + 2);

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

  if (s== 0)
    mpfr_set_ui (dest, 0, GMP_RNDN);
  else {
    mpfr_set_si (dest, s, GMP_RNDN);
    mpfr_set_exp (dest, gexpo ((GEN) fp)+1);
    copy_reverse (dest->_mpfr_d,
                  & ((long*)fp)[2],
                  (prec-1)/(8*sizeof (mp_limb_t))+1);
    /* We may have copied more bit than necessary: reprec it */
    p = mpfr_get_prec (dest);
    dest->_mpfr_prec = ((prec-1)/(8*sizeof (mp_limb_t))+1)
      * (8*sizeof(mp_limb_t));
    mpfr_prec_round (dest, p, mpcheck_rnd_mode);
  }
}

static void get_fp (void *fp, mpfr_srcptr src)
{
  if (mpfr_nan_p (src) || mpfr_inf_p (src))
    {
      printf ("Error: Pari can't handle NAN or INF\n");
      exit (1);
    }
  if (!mpfr_zero_p (src)) {
    setsigne (fp, mpfr_sgn (src));
    setexpo (fp, mpfr_get_exp (src)-1);
    copy_reverse (& ((long*)fp)[2], src->_mpfr_d,
                  (prec-1)/(8*sizeof (mp_limb_t))+1);
  } else {
    gsubz (fp, fp, fp);
  }
}

static int set_rnd_mode (mp_rnd_t rnd) {
  return 1; /* Pari has no Rounding Mode */
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

static mpcheck_user_func_t tab[] = {
  {"add", my_pari_add, 0, 0},
  {"add", my_pari_add, LONG_MAX, LONG_MAX},
  {"sub", my_pari_sub, LONG_MAX, LONG_MAX},
  {"sub", my_pari_sub, 0, 0},
  {"mul", my_pari_mul, 0, 0},
  {"mul", my_pari_mul, LONG_MAX-1, LONG_MAX-1},
  {"div", my_pari_div, 0, 0},
  {"div", my_pari_div, LONG_MAX, LONG_MAX},
  {NULL, NULL, 0, 0}
};


int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;

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
  /*printf ("X10="); mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
    printf ("\nFp ="); output ((GEN) fp);*/
  del_fp (fp);
  prec = 0;

  mpcheck_init (argc, argv, 53, -1021, 1021,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		1, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
  return 0;
}

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "PARI not available.\n");
  return 0;
}

#endif
