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

#ifndef __MPCHECK__
#define __MPCHECK__

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "gmp.h"
#include "mpfr.h"

#if defined (__cplusplus)
extern "C" {
#endif

  typedef enum {
    MPCHECK_TEST_RANGE=1, MPCHECK_TEST_MONOTON=2, MPCHECK_TEST_SYMM=4
  } mpcheck_test_e;
  
  typedef enum {
    RANGE_INF=1, RANGE_ZERO, RANGE_ONE, RANGE_PI2, RANGE_TWO, RANGE_PI
  } mpcheck_range_e;
  
  typedef struct {
    const char *const name;      /* Name of the function */
    int        (*mpfr)();        /* MPFR function */
    int         NumArg;          /* # of args */
    mpcheck_range_e min, max;    /* Range */
    int monoton;
    int  symm;
    int signed_input;
  } mpcheck_func_t;
  
  typedef struct {
    const char *const name;
    void (*func) (void*,const void*, const void*);
    mp_exp_t e1, e2;
  } mpcheck_user_func_t;

  extern mpcheck_func_t mpcheck_tab[];
  extern mp_rnd_t mpcheck_rnd_mode;

/* useful for testing monotonicity and symmetry */
#define NO_MONOTON 0
#define INCREASING 1
#define DECREASING -1

#define NO_SYMM 0
#define ODD 3
#define EVEN 4

#define IN_POS 0
#define IN_POSNEG 1

#define ALL_RND  15
#define ALL_TEST ((mpcheck_test_e) 7)

#define DEFAULT_N 10000 /* default number of tests per function */

  void mpcheck_init (int argc, const char *const argv[],
		     mp_prec_t prec2, mp_exp_t emin2, mp_exp_t emax2,
		     void *(*new2)(mp_prec_t), void (*del2)(void *),
		     void (*getfp2)(void *, mpfr_srcptr),
		     void (*setfp2)(mpfr_ptr, const void *),
		     int (*setrnd)(mp_rnd_t),
		     int rnd_mode2, mpcheck_test_e test2, unsigned long seed2,
		     unsigned long N2, int verbose2);
  void mpcheck_clear (FILE *out);
  
  void mpcheck (FILE *out, mp_exp_t e1, mp_exp_t e2,
		const char *const name, 
		void (*func) (void*,const void*, const void*));
  void mpcheck_check (FILE *out, mpcheck_user_func_t *tab);


#if defined (__cplusplus)
}
#endif

#endif 
