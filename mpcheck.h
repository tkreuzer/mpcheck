/*
 * MPCHECK - Check double LIBM functions
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

/* Define Precision to check */
#ifndef FPPREC
# define FPPREC 53
#endif

/* define MATHLIB to check library mathlib */
/* #define MATHLIB */

#if (FPPREC==113 && defined(__ia64))
#define __LONG_DOUBLE_MATH
#endif

/* Check MATHLIB */
#ifdef MATHLIB
# if (FPPREC != 53)
# error "MathLib only works for double precision"
# endif 
#endif 

/* Define libmname */
#ifdef MATHLIB
# define libmname(fct)  u ## fct
#elif (FPPREC==113 && defined(__ia64))
# define libmname(fct)  fct ## l
#else
# define libmname(fct)  fct
#endif

/* Define internal type to use */
#if (FPPREC == 24)
# define fptype float
# define EMIN -125
# define EMAX 128
# define PiRoundedUp 3.141592741
#elif (FPPREC == 53)
# define fptype double
# define EMIN -1021
# define EMAX 1024
# define PiRoundedUp 3.1415926535897935601
#elif (FPPREC == 64)
# if defined(__ia64)
#  define fptype extended
#  define _FPWIDETYPES /* for HP-UX */
# else
#  define fptype long double
# endif
# define EMIN -16381
# define EMAX 16384
# define PiRoundedUp 3.141592653589793238462643
#elif (FPPREC == 113)
# if defined(__ia64)
#  define _FPWIDETYPES /* for HP-UX */
#  define fptype long double
# else
#  define fptype quad
# endif
# define EMIN -16381
# define EMAX 16384
# define PiRoundedUp 3.141592653589793238462643383279502884197
#else
# error "FPPREC not yet implemented"
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "gmp.h"
#include "mpfr.h"

#ifdef MATHLIB
#include "MathLib.h"
#endif

#if defined(HAVE_TGAMMA) && !defined(tgamma)
fptype tgamma (fptype);
#endif
#if defined(HAVE_LOG2) && !defined(log2)
fptype log2 (fptype);
#endif
#if defined(HAVE_EXP2) && !defined(exp2)
fptype exp2 (fptype);
#endif

/* stolen from mpfr-impl.h (Ugly) */
#ifndef MPFR_EXP
# define MPFR_EXP(x) ((x)->_mpfr_exp)
#endif

/* MATHLIB only implements rounding to nearest */
#ifdef MATHLIB
# define MAX_RND 1
#else
# define MAX_RND 4
#endif

extern void test (const char *, mp_exp_t, mp_exp_t, unsigned long, unsigned long);
extern double ulp (fptype y, mpfr_ptr z);
extern void testall (unsigned long, unsigned long);
extern void setup (void);
extern void set_rnd_mode (mp_rnd_t r);

extern int (*set_fp) (mpfr_ptr, fptype, mp_rnd_t);
extern fptype (*get_fp) (mpfr_srcptr, mp_rnd_t);
extern void (*fprint_fp) (FILE *, fptype);
extern double MAX_ERR_NEAR;
extern double MAX_ERR_DIR;
extern int verbose ;
extern int test_monotonicity, test_range, test_symmetry ;

#define print_fp(x) (*fprint_fp)(stdout, x)

typedef struct {
  int         NumArg;      /* # of args */
  const char *name;        /* Name of the function */
  fptype     (*libm)();    /* Function to check */
  int        (*mpfr)();    /* Equivalent MPFR function */
  fptype      min, max;    /* Range */
  int         monoton;     /* Monotonicity property */
  int         symm;        /* symmetry property */
} mpcheck_func_t;

/* useful for testing monotonicity and symmetry */
#define NO_MONOTON 0
#define INCREASING 1
#define DECREASING -1

#define NO_SYMM 0
#define ODD 3
#define EVEN 4

#define numberof(x) (sizeof(x)/sizeof(x[0]))

extern mpcheck_func_t  FuncTab[];
extern int FuncTabNum;


/* Log level */
#define LOG(v, code) ( verbose >= (v) ? (void) (code) : (void) 0)

#endif 
