/*
 * MPCHECK - Check mathematical functions
 * Copyright (C) 2013 INRIA
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

#include "mpfr.h"

/*  Setup for double checks */
const mp_prec_t prec = 53;
const mp_exp_t  emin = -1021;
const mp_exp_t  emax = 1024;
typedef double fptype;
#define fpstr   "double"
#define FP_KERNEL_SETUP "#pragma OPENCL EXTENSION cl_khr_fp64: enable\n"

#include "mpcheck-opencl.c"
