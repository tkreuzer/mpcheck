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

#if defined(HAVE_CRLIBM)

#include "mpcheck.h"

#include "crlibm.h"

#include "rnd_mode.c"

#define LIB_INIT() 
#define LIB_EXIT() 
#define fptype double

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mpfr_set_d (dest, *(double*)fp, GMP_RNDN);
}
static void get_fp (void *fp, mpfr_srcptr src)
{
  *(double*) fp = mpfr_get_d (src, GMP_RNDN);
}

static void my_exp (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode) 
    {
    case GMP_RNDN:
      *(double*) dest = exp_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = exp_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = exp_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = exp_rd (*(double*)a);
      break;
    default:
      abort ();
      break;
    }
  crlibm_exit (mode);
}

static void my_log (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = log_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = log_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = log_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = log_rd (*(double*)a);
      break;
    default:
      abort ();
      break;
    }
  crlibm_exit (mode);
}

static void my_sin (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = sin_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = sin_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = sin_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = sin_rd (*(double*)a);
      break;
    default:
      abort ();
      break;
    }
  crlibm_exit (mode);
}

static void my_cos (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = cos_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = cos_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = cos_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = cos_rd (*(double*)a);
      break;
    default:
      break;
    }
  crlibm_exit (mode);
}

static void my_tan (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = tan_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = tan_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = tan_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = tan_rd (*(double*)a);
      break;
    default:
      break;
    }
  crlibm_exit (mode);
}

static void my_atan (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = atan_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = atan_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = atan_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = atan_rd (*(double*)a);
      break;
    default:
      break;
    }
  crlibm_exit (mode);
}

static void my_sinh (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = sinh_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = sinh_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = sinh_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = sinh_rd (*(double*)a);
      break;
    default:
      break;
    }
  crlibm_exit (mode);
}

static void my_cosh (void *dest, const void *a, const void *b)
{
  int mode = crlibm_init ();
  switch (mpcheck_rnd_mode)
    {
    case GMP_RNDN:
      *(double*) dest = cosh_rn (*(double*)a);
      break;
    case GMP_RNDZ:
      *(double*) dest = cosh_rz (*(double*)a);
      break;
    case GMP_RNDU:
      *(double*) dest = cosh_ru (*(double*)a);
      break;
    case GMP_RNDD:
      *(double*) dest = cosh_rd (*(double*)a);
      break;
    default:
      break;
    }
  crlibm_exit (mode);
}

static mpcheck_user_func_t tab[] = {
  {"exp", my_exp, 0, 0},
  {"exp", my_exp, 9, 0}, 
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
  {"sinh", my_sinh, 0, 0},
  {"sinh", my_sinh, 9, 0}, 
  {"cosh", my_cosh, 0, 0},
  {"cosh", my_cosh, 9, 0}, 

  {NULL, NULL, 0, 0}
};

#include "native.c"

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "CRLIBM not available.\n");
  return 0;
}

#endif
