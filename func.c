/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2002, 2004, 2005, 2010, 2014 INRIA
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

int mpfr_Lgamma (mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  int s;
  return mpfr_lgamma (y, &s, x, r);
}

int mpfr_j17 (mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_jn (y, 17, x, r);
}

int mpfr_j42 (mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_jn (y, 42, x, r);
}

int mpfr_y17 (mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_yn (y, 17, x, r);
}

int mpfr_y42 (mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_yn (y, 42, x, r);
}

mpcheck_func_t  mpcheck_tab[] = {
  /* name,  MPFR name, number of arguments, left bound of output range,
     right bound of output range, monotonicity, symmetry, input range */
  {"add", mpfr_add, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POSNEG},
  {"sub", mpfr_sub, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POSNEG},
  {"mul", mpfr_mul, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POSNEG},
  {"div", mpfr_div, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POSNEG},

  {"sqrt", mpfr_sqrt, 1, -RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM, IN_POS},
  {"exp", mpfr_exp, 1, RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM, IN_POSNEG},
  {"log", mpfr_log, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM, IN_POSNEG},
  {"dilog", mpfr_li2, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM, IN_POSNEG},
  {"eint", mpfr_eint, 1, RANGE_ZERO, RANGE_INF, DECREASING, NO_SYMM, IN_POS},

  {"sin", mpfr_sin, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, ODD, IN_POSNEG},
  {"sincos1", mpfr_sin, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, ODD, IN_POSNEG},
  {"cos", mpfr_cos, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, EVEN, IN_POSNEG},
  {"sincos2", mpfr_cos, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, EVEN, IN_POSNEG},
  {"tan", mpfr_tan, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD, IN_POSNEG},
  {"asin", mpfr_asin, 1, -RANGE_PI2, RANGE_PI2, INCREASING, ODD, IN_POSNEG},
  {"acos", mpfr_acos, 1, RANGE_ZERO, RANGE_PI, DECREASING, NO_SYMM,IN_POSNEG},
  {"atan", mpfr_atan, 1, -RANGE_PI2, RANGE_PI2, INCREASING, ODD, IN_POSNEG},
  {"atan2", mpfr_atan2, 2, -RANGE_PI, RANGE_PI, NO_MONOTON, NO_SYMM,IN_POSNEG},

  {"sinh", mpfr_sinh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD, IN_POSNEG},
  {"cosh", mpfr_cosh, 1, RANGE_ONE, RANGE_INF, NO_MONOTON, EVEN, IN_POSNEG},
  {"tanh", mpfr_tanh, 1, -RANGE_ONE, RANGE_ONE, INCREASING, ODD, IN_POSNEG},
  {"asinh", mpfr_asinh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD, IN_POS},
  {"acosh", mpfr_acosh, 1, RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM, IN_POS},
  {"atanh", mpfr_atanh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD, IN_POS},
  {"cot", mpfr_cot, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD, IN_POSNEG},
  {"coth", mpfr_coth, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD, IN_POSNEG},
 
  {"exp2", mpfr_exp2, 1, RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM,IN_POSNEG},
  {"log2", mpfr_log2, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM,IN_POSNEG},
  {"expm1", mpfr_expm1, 1, -RANGE_ONE, RANGE_INF,INCREASING,NO_SYMM,IN_POSNEG},
  {"exp10", mpfr_exp10, 1, RANGE_ZERO, RANGE_INF,INCREASING,NO_SYMM,IN_POSNEG},
  {"log10", mpfr_log10, 1, -RANGE_INF, RANGE_INF,INCREASING,NO_SYMM,IN_POSNEG},
  {"log1p", mpfr_log1p, 1, -RANGE_INF, RANGE_INF,INCREASING,NO_SYMM,IN_POSNEG},

  {"gamma", mpfr_gamma, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM,IN_POS},
  {"gamma_inc", mpfr_gamma_inc, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"lgamma", mpfr_Lgamma, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM,IN_POSNEG},
  {"lngamma", mpfr_lngamma, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM,IN_POSNEG},
  {"cbrt", mpfr_cbrt, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD, IN_POSNEG},
  {"erf", mpfr_erf, 1, -RANGE_ONE, RANGE_ONE, INCREASING, ODD, IN_POSNEG},
  {"erfc", mpfr_erfc, 1, RANGE_ZERO, RANGE_TWO, DECREASING, NO_SYMM,IN_POSNEG},

  {"j0", mpfr_j0, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, EVEN, IN_POSNEG},
  {"j1", mpfr_j1, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD,  IN_POSNEG},
  {"j17", mpfr_j17, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD,  IN_POSNEG},
  {"j42", mpfr_j42, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, EVEN,  IN_POSNEG},
  {"y0", mpfr_y0, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"y1", mpfr_y1, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"y17", mpfr_y17, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"y42", mpfr_y42, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"zeta", mpfr_zeta, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POSNEG},

  {"hypot", mpfr_hypot, 2, -RANGE_INF, RANGE_INF,NO_MONOTON,NO_SYMM,IN_POSNEG},
  {"pow", mpfr_pow, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {"agm", mpfr_agm, 2, RANGE_ZERO, RANGE_INF, NO_MONOTON, NO_SYMM, IN_POS},
  {NULL, NULL, 0, 0, 0, 0, 0, 0}
};



