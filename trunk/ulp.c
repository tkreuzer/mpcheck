/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

#include "mpcheck.h"

/* computes the error in ulps between y [computed by libm]
   and z [computed by mpfr]. It assumes y <> z. */
double ulp (fptype y, mpfr_ptr z)
{
  mpfr_t yy;
  double u;
  
  mpfr_init2 (yy, 2*FPPREC);

  (*set_fp) (yy, y, GMP_RNDN); /* Exact */
  mpfr_sub (yy, yy, z, GMP_RNDN);
  MPFR_EXP (yy) += FPPREC - MPFR_EXP (z);
  u = mpfr_get_d (yy, GMP_RNDN);

  mpfr_clear (yy);
  return u;
}

