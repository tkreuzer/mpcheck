/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

#include "mpcheck.h"

fptype add (fptype x, fptype y)
{
  return x + y;
}

fptype sub (fptype x, fptype y)
{
  return x - y;
}

fptype mul (fptype x, fptype y)
{
  return x * y;
}

fptype my_div (fptype x, fptype y)
{
  return x / y;
}

#ifdef HAVE_TGAMMA
fptype my_gamma (fptype x)
{
  return libmname(tgamma)(x);
}
#endif

/* Tiny macro */
#define func(name) #name, (fptype (*)()) libmname(name), mpfr_ ## name

mpcheck_func_t  FuncTab[] = {
  {1, func(exp), 0.0, 1.0/0.0, INCREASING, NO_SYMM},
#ifdef HAVE_EXP2
  {1, func(exp2), 0.0, 1.0/0.0, INCREASING, NO_SYMM},
#endif
  {1, func(log), -1.0/0.0, 1.0/0.0, INCREASING, NO_SYMM},
#ifdef HAVE_LOG2
  {1, func(log2), -1.0/0.0, 1.0/0.0, INCREASING, NO_SYMM},
#endif
  {1, func(sin), -1.0, 1.0, NO_MONOTON, ODD},
  {1, func(cos), -1.0, 1.0, NO_MONOTON, EVEN},
  {1, func(tan), -1.0/0.0, 1.0/0.0, NO_MONOTON, ODD},
  {1, func(asin), -PiRoundedUp/2.0, PiRoundedUp/2.0, INCREASING, ODD},
  {1, func(acos), 0.0, PiRoundedUp, DECREASING, NO_SYMM},
  {1, func(atan), -PiRoundedUp/2.0, PiRoundedUp/2.0, INCREASING, ODD},
  {1, func(sqrt), 0.0, 1.0/0.0, INCREASING, NO_SYMM},
  {2, func(pow), -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
  {2, func(add), -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
  {2, func(sub), -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
  {2, func(mul), -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
  {2, "div",  (fptype (*)()) my_div, mpfr_div, -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
#ifndef MATHLIB
  {1, func(expm1), -1.0, 1.0/0.0, INCREASING, NO_SYMM},
  {1, func(log10), -1.0/0.0, 1.0/0.0, INCREASING, NO_SYMM},
  {1, func(log1p), -1.0/0.0, 1.0/0.0, INCREASING, NO_SYMM},
  {1, func(sinh), -1.0/0.0, 1.0/0.0, INCREASING, ODD},
  {1, func(cosh), 1.0, 1.0/0.0, NO_MONOTON, EVEN},
  {1, func(tanh), -1.0, 1.0, INCREASING, ODD},
  {1, func(asinh), -1.0/0.0, 1.0/0.0, INCREASING, ODD},
  {1, func(acosh), 0.0, 1.0/0.0, INCREASING, NO_SYMM},
  {1, func(atanh), -1.0/0.0, 1.0/0.0, INCREASING, ODD},
#ifdef HAVE_TGAMMA
  {1, "gamma", (fptype (*)()) my_gamma, mpfr_gamma, -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM},
#endif
  {1, func(cbrt), -1.0/0.0, 1.0/0.0, INCREASING, ODD},
  {2, func(hypot), -1.0/0.0, 1.0/0.0, NO_MONOTON, NO_SYMM}
#endif
};

int FuncTabNum = numberof (FuncTab);

