/*
 * MPCHECK - Set Machine Rounding Mode
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

/* Functions to set/get the machine rounding mode */
#if defined (USE_FENV_H)
# include <fenv.h>
#ifdef FE_TONEAREST
# define TONEAREST fesetround(FE_TONEAREST)
#endif
#ifdef FE_TOWARDZERO
# define TOZERO    fesetround(FE_TOWARDZERO)
#endif
#ifdef FE_UPWARD
# define TOINFP    fesetround(FE_UPWARD)
#endif
#ifdef FE_DOWNWARD
# define TOINFM    fesetround(FE_DOWNWARD)
#endif

#elif defined (USE_SYS_FPU_H)
# include <sys/fpu.h>
extern int swapRM();
# define TOZERO swapRM(ROUND_TO_ZERO)
# define TOINFP swapRM(ROUND_TO_PLUS_INFINITY)
# define TONEAREST swapRM(ROUND_TO_NEAREST)
# define TOINFM swapRM(ROUND_TO_MINUS_INFINITY)

#elif defined (USE_MATH_H)
# include <math.h>
# define TOZERO fpsetround(FP_RZ)
# define TOINFP fpsetround(FP_RP)
# define TONEAREST fpsetround(FP_RN)
# define TOINFM fpsetround(FP_RM)

#elif defined (USE_FLOATINGPOINT_H_1)
# include <floatingpoint.h>
# define TOZERO fpsetround(FP_RZ)
# define TOINFP fpsetround(FP_RP)
# define TONEAREST fpsetround(FP_RN)
# define TOINFM fpsetround(FP_RM)

#elif defined (USE_IEEEFP_H)
# include <ieeefp.h>
# define TOZERO fpsetround(FP_RZ)
# define TOINFP fpsetround(FP_RP)
# define TONEAREST fpsetround(FP_RN)
# define TOINFM fpsetround(FP_RM)

#elif defined (USE_FLOAT_H_1)
# include <float.h>
# define TOZERO write_rnd(FP_RND_RZ)
# define TOINFP write_rnd(FP_RND_RP)
# define TONEAREST write_rnd(FP_RND_RN)
# define TOINFM write_rnd(FP_RND_RM)

#elif defined (USE_FLOAT_H_2)
# include <float.h>
# define FP_RND_RZ       0
# define FP_RND_RN       1
# define FP_RND_RP       2
# define FP_RND_RM       3
# define TOZERO write_rnd(FP_RND_RZ)
# define TOINFP write_rnd(FP_RND_RP)
# define TONEAREST write_rnd(FP_RND_RN)
# define TOINFM write_rnd(FP_RND_RM)

#elif defined (USE_FLOAT_H_3)
# include <float.h>
# define FP_RND_RZ       0
# define FP_RND_RN       1
# define FP_RND_RP       2
# define FP_RND_RM       3
# define TOZERO fp_swap_rnd(FP_RND_RZ)
# define TOINFP fp_swap_rnd(FP_RND_RP)
# define TONEAREST fp_swap_rnd(FP_RND_RN)
# define TOINFM fp_swap_rnd(FP_RND_RM)

#elif defined (USE_FLOATINGPOINT_H_2)
# include <floatingpoint.h>
static char *out;
# define TOZERO ieee_flags("set","direction","tozero",&out)
# define TOINFP ieee_flags("set","direction","positive",&out)
# define TONEAREST ieee_flags("set","direction","nearest",&out)
# define TOINFM ieee_flags("set","direction","negative",&out)

#elif defined (USE_FPUCONTROL_H)
# include <fpu_control.h>
# define TOZERO _FPU_SETCW(_FPU_RC_ZERO)
# define TOINFP _FPU_SETCW(_FPU_RC_UP)
# define TOINFM _FPU_SETCW(_FPU_RC_DOWN)
# define TONEAREST _FPU_SETCW(_FPU_RC_NEAREST)

#elif defined (USE_GCC_I386ASM)
# define _FPU_EXTENDED 0x300
# define _FPU_DOUBLE   0x200
# define _FPU_DEFAULT  0x137f
# define _FPU_RC_NEAREST 0x0
# define _FPU_RC_DOWN    0x400
# define _FPU_RC_UP      0x800
# define _FPU_RC_ZERO    0xC00
# define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
# define _fpu_ieee ((_FPU_DEFAULT & (~_FPU_EXTENDED)) | _FPU_DOUBLE)
# define TOZERO __setfpucw(_fpu_ieee | _FPU_RC_ZERO)
# define TOINFP __setfpucw(_fpu_ieee | _FPU_RC_UP)
# define TOINFM __setfpucw(_fpu_ieee | _FPU_RC_DOWN)
# define TONEAREST __setfpucw(_fpu_ieee)

#else
# error "Can't set machine rounding modes."
#endif

/* sets the machine rounding mode to the value rnd_mode */
int 
set_rnd_mode(mp_rnd_t rnd_mode)
{
  switch (rnd_mode) {
  case GMP_RNDN: 
#ifdef TONEAREST
    TONEAREST;
    return 1;
#else
    return 0;
#endif
  case GMP_RNDZ: 
#ifdef TOZERO
    TOZERO;
    return 1;
#else
    return 0;
#endif
  case GMP_RNDU:
#ifdef TOINFP
    TOINFP;
    return 1;
#else
    return 0;
#endif
  case GMP_RNDD: 
#ifdef TOINFM
    TOINFM; 
    return 1;
#else
    return 0;
#endif
  default:
    fprintf(stderr, "invalid rounding mode\n"); 
    abort ();
  }
  return 1;
}

#ifdef WITHIN_CONFIGURE
int main ()
{
  set_rnd_mode (GMP_RNDD);
  set_rnd_mode (GMP_RNDU);
  set_rnd_mode (GMP_RNDZ);
  set_rnd_mode (GMP_RNDN);
  return 0;
}
#endif
