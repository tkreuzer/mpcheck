/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2002, 2004, 2005 INRIA
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

mpcheck_func_t  mpcheck_tab[] = {
  {"add", mpfr_add, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {"sub", mpfr_sub, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {"mul", mpfr_mul, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {"div", mpfr_div, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},

  {"sqrt", mpfr_sqrt, 1, -RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM},
  {"exp", mpfr_exp, 1, -RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM},
  {"log", mpfr_log, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM},

  {"sin", mpfr_sin, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, ODD},
  {"cos", mpfr_cos, 1, -RANGE_ONE, RANGE_ONE, NO_MONOTON, EVEN},
  {"tan", mpfr_tan, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, ODD},
  {"asin", mpfr_asin, 1, -RANGE_PI2, RANGE_PI2, INCREASING, ODD},
  {"acos", mpfr_acos, 1, -RANGE_ZERO, RANGE_PI2, DECREASING, NO_SYMM},
  {"atan", mpfr_atan, 1, -RANGE_PI2, RANGE_PI2, INCREASING, ODD},

  {"sinh", mpfr_sinh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD},
  {"cosh", mpfr_cosh, 1, RANGE_ONE, RANGE_INF, NO_MONOTON, EVEN},
  {"tanh", mpfr_tanh, 1, -RANGE_ONE, RANGE_ONE, INCREASING, ODD},
  {"asinh", mpfr_asinh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD},
  {"acosh", mpfr_acosh, 1, -RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM},
  {"atanh", mpfr_atanh, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD},
 
  {"exp2", mpfr_exp2, 1, RANGE_ZERO, RANGE_INF, INCREASING, NO_SYMM},
  {"log2", mpfr_log2, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM},
  {"expm1", mpfr_expm1, 1, -RANGE_ONE, RANGE_INF, INCREASING, NO_SYMM},
  {"log10", mpfr_log10, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM},
  {"log1p", mpfr_log1p, 1, -RANGE_INF, RANGE_INF, INCREASING, NO_SYMM},

  {"gamma", mpfr_gamma, 1, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {"cbrt", mpfr_cbrt, 1, -RANGE_INF, RANGE_INF, INCREASING, ODD},

  {"hypot", mpfr_hypot, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {"pow", mpfr_pow, 2, -RANGE_INF, RANGE_INF, NO_MONOTON, NO_SYMM},
  {NULL, NULL, 0, 0, 0, 0, 0}
};


