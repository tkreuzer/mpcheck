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

#if defined(HAVE_OPEN_CL)

#include "mpcheck.h"

#include "CL/cl.h"

/*  Setup for float checks if double not already defined */
#ifndef fpstr
const mp_prec_t prec = 24;
const mp_exp_t  emin = -125;
const mp_exp_t  emax = 128;
typedef float fptype;
#define fpstr           "float"
#define FP_KERNEL_SETUP ""
#endif

cl_platform_id platform;
cl_device_id   device;
cl_context     context;
cl_program     program;
cl_command_queue queue;

/* Open CL embedded kernel source program (Initial kernel for all operators) */
/* Note: it seems that for Nvidia device, the more function are in the kernel,
   the slower the kernel is, even if the code is not used */
const char opencl_kernel_program[] = FP_KERNEL_SETUP
  "__kernel void opencl_add(__global " fpstr " *dest,                \n"
  "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
  "  dest[0] = a[0] + b[0];                                      \n"
  "}                                                             \n"
  "__kernel void opencl_sub(__global " fpstr " *dest,                \n"
  "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
  "  dest[0] = a[0] - b[0];                                      \n"
  "}                                                             \n"
  "__kernel void opencl_mul(__global " fpstr " *dest,                \n"
  "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
  "  dest[0] = a[0] * b[0];                                      \n"
  "}                                                             \n"
  "__kernel void opencl_div(__global " fpstr " *dest,                \n"
  "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
  "  dest[0] = a[0] / b[0];                                      \n"
  "}                                                             \n"
;

const char *opencl_kernel_program_ptr[1] = { opencl_kernel_program};

/* Run & Check an opencl function command */
#define RUN_CHECK_ERR(err, instruction)					\
  do {									\
    instruction;							\
    if ((err) != CL_SUCCESS) {						\
      fprintf (stderr, "ERROR: FILE = " __FILE__ " LINE = %d ERR=%d\n"	\
	       "Failure in calling : " #instruction "\n", __LINE__, (err)); \
      exit(1);								\
    }									\
  } while (0)

static void
init_opencl(void)
{
  cl_int err;
  cl_build_status build_status;
  char buffer[1024];

  /* Get platform */
  RUN_CHECK_ERR(err,
		err = clGetPlatformIDs(1, &platform, NULL));

  /* Get device */
  /* TODO: How to interface the way to select a specific device ? */
  RUN_CHECK_ERR(err,
		err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, NULL));
  RUN_CHECK_ERR(err,
		err = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof buffer -1,
				      buffer, NULL));
  buffer[sizeof buffer -1] = 0;
  printf ("[Selected Open CL device = %s]\n", buffer);				      

  /* Create context */
  RUN_CHECK_ERR(err,
		context = clCreateContext(NULL, 1, &device, NULL, NULL, &err));

  /* Build program */
  RUN_CHECK_ERR (err,
                 program = clCreateProgramWithSource(context, 1, 
						     opencl_kernel_program_ptr, NULL, &err));
  /* Do not RUN_CHECK_ERR: let the fallthrough handle the error */
  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

  /* Check success of build */
  do {
    RUN_CHECK_ERR (err,
		   err =  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS,
						sizeof(cl_build_status), &build_status,
						NULL));
  } while (build_status == CL_BUILD_IN_PROGRESS);
  if (build_status != CL_BUILD_SUCCESS)
    {
      /* Fail to build: display log info */
      size_t log_size;
      char *program_log;
      fprintf (stderr, "ERROR: Failed to build program (build_status=%d):\n", 
	       (int) build_status);
      /* Find size of log and print to std output */
      clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
			    0, NULL, &log_size);
      program_log = malloc(log_size + 1);
      if (program_log != NULL)
	{
	  program_log[log_size] = '\0';
	  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
				log_size + 1, program_log, NULL);
	  fprintf(stderr, "LOG[%d]=%s\n", (int) log_size, program_log);
	  free(program_log);
	}
      exit(1);
    }
  /* Create a command queue */
  RUN_CHECK_ERR (err,
		 queue = clCreateCommandQueue(context, device, 0, &err) );
}

static void
quit_opencl(void)
{
  clReleaseCommandQueue(queue);
  clReleaseProgram(program);
  clReleaseContext(context);
}

/* Execute the opencl function by giving the argument to the kernel function and calling it */
/* Note: This is far from being the more efficient way to do the test */
static void
run_opencl_function(const char *function_name,
		    void *dest, const void *a, const void *b)
{
  cl_int err;
  size_t global_size, local_size;
  cl_mem         input_buffer1, input_buffer2, output_buffer;
  cl_kernel kernel;

  /* Create data buffer */
  RUN_CHECK_ERR (err,
		 input_buffer1 = clCreateBuffer(context, CL_MEM_READ_ONLY |
						CL_MEM_COPY_HOST_PTR, sizeof(fptype), (void*) a, &err));
  RUN_CHECK_ERR (err, 
		 input_buffer2 = clCreateBuffer(context, CL_MEM_READ_ONLY |
						CL_MEM_COPY_HOST_PTR, sizeof(fptype), (void*) b, &err));
  RUN_CHECK_ERR (err,
		 output_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(fptype), dest, &err));

  /* Create a kernel */
  RUN_CHECK_ERR (err,
		 kernel = clCreateKernel(program, function_name, &err));

  /* Set up kernel arguments */
  RUN_CHECK_ERR (err,
		 err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buffer));
  RUN_CHECK_ERR (err,
		 err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_buffer1));
  RUN_CHECK_ERR (err,
		 err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &input_buffer2));

  /* Enqueue kernel */
  global_size = 1;
  local_size = 1;
  RUN_CHECK_ERR (err,
		 err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, 
					      &local_size, 0, NULL, NULL));

  /* Wait to read the output */
  RUN_CHECK_ERR (err,
		 err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0, 
					   sizeof(fptype), dest, 0, NULL, NULL));

  /* Deallocate ressources */
  clReleaseKernel(kernel);
  clReleaseMemObject(input_buffer1);
  clReleaseMemObject(input_buffer2);
  clReleaseMemObject(output_buffer);
}

/* Recreate a kernel for specific function call */
static void
build_and_run_opencl_function(const char *function_name, int arg_count,
			      void *dest, const void *a, const void *b)
{
  static char previous_function[50];
  static char opencl_name[5+sizeof previous_function];
  
  /* Cache mechanism */
  if (strcmp (function_name, previous_function) != 0) {
    char buffer[1024];
    cl_int err;
    cl_build_status build_status;
    const char *buffer_ptr[1] = { buffer };

    strncpy (previous_function, function_name, sizeof previous_function);

    if (arg_count == 1) {
      sprintf (buffer, FP_KERNEL_SETUP
	       "__kernel void opencl_%s(__global " fpstr " *dest,               \n"
	       "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
	       "  dest[0] = %s(a[0]);                                       \n"
	       "}                                                             \n",
	       function_name, function_name);
    } else {
      sprintf (buffer, FP_KERNEL_SETUP
	       "__kernel void opencl_%s(__global " fpstr " *dest,               \n"
	       "			  __global " fpstr " *a,__global " fpstr " *b) { \n"
	       "  dest[0] = %s(a[0], b[0]);                                 \n"
	       "}                                                           \n",
	       function_name, function_name);
    }

    clReleaseProgram(program);

    /* Build program */
    RUN_CHECK_ERR (err,
		   program = clCreateProgramWithSource(context, 1, 
						       buffer_ptr, NULL, &err));
    RUN_CHECK_ERR(err,
		  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL));

    /* Check success of build */
    do {
      RUN_CHECK_ERR (err,
		     err =  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS,
						  sizeof(cl_build_status), &build_status,
						  NULL));
    } while (build_status == CL_BUILD_IN_PROGRESS);
    if (build_status != CL_BUILD_SUCCESS)
      {
	size_t log_size;
	char *program_log;
	fprintf (stderr, "ERROR: Failed to build program (build_status=%d):\n", 
		 (int) build_status);
	/* Find size of log and print to std output */
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
			      0, NULL, &log_size);
	program_log = malloc(log_size + 1);
	if (program_log != NULL)
	  {
	    program_log[log_size] = '\0';
	    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
				  log_size + 1, program_log, NULL);
	    fprintf(stderr, "LOG[%d]=%s\n", (int) log_size, program_log);
	    free(program_log);
	    exit(1);
	  }
      }
    sprintf (opencl_name, "opencl_%s", function_name);
  } /* cached */

  run_opencl_function(opencl_name, dest, a, b);
}


/* Create a new FP number of prec p */
static void *new_fp (mp_prec_t p)
{
  void *fp;
  if (p != prec)
    {
      fprintf (stderr, "ERROR: Requested prec: %lu. FPTYPE prec: %lu\n",
	       p, prec);
      abort ();
    }
  fp = malloc (sizeof (fptype));
  if (fp == NULL)
    {
      fprintf (stderr, "ERROR: Can't allocate memory");
      abort ();
    }
  return fp;
}

static void del_fp (void *fp)
{
  free (fp);
}

static void set_fp (mpfr_ptr dest, const void *fp)
{
  mpfr_set_d (dest, (double) (*(fptype*)fp), GMP_RNDN);
}

static void get_fp (void *fp, mpfr_srcptr src)
{
  /* There should be no rounding since src has the prec of a float */
  *(fptype*) fp = mpfr_get_d (src, GMP_RNDN);
}

/* Even if the device supports other rounding mode,
   there is no syntax to select the rounding mode in the kernel with OpenCL. */
static int set_rnd_mode (mp_rnd_t rnd) {
  return rnd == GMP_RNDN;
}

/* Start Function to check */
void my_opencl_add (void *dest, const void *src1, const void *src2) {
  run_opencl_function("opencl_add", dest, src1, src2);
}
void my_opencl_sub (void *dest, const void *src1, const void *src2) {
  run_opencl_function("opencl_sub", dest, src1, src2);
}
void my_opencl_mul (void *dest, const void *src1, const void *src2) {
  run_opencl_function("opencl_mul", dest, src1, src2);
}
void my_opencl_div (void *dest, const void *src1, const void *src2) {
  run_opencl_function("opencl_div", dest, src1, src2);
}
void my_opencl_sqrt (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("sqrt", 1, dest, src1, src2);
}
void my_opencl_exp (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("exp", 1, dest, src1, src2);
}
void my_opencl_log (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("log", 1, dest, src1, src2);
}
void my_opencl_sin (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("sin", 1, dest, src1, src2);
}
void my_opencl_cos (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("cos", 1, dest, src1, src2);
}
void my_opencl_tan (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("tan", 1, dest, src1, src2);
}
void my_opencl_asin (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("asin", 1, dest, src1, src2);
}
void my_opencl_acos (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("acos", 1, dest, src1, src2);
}
void my_opencl_atan (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("atan", 1, dest, src1, src2);
}
void my_opencl_atan2 (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("atan2", 2, dest, src1, src2);
}
void my_opencl_sinh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("sinh", 1, dest, src1, src2);
}
void my_opencl_cosh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("cosh", 1, dest, src1, src2);
}
void my_opencl_tanh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("tanh", 1, dest, src1, src2);
}
void my_opencl_asinh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("asinh", 1, dest, src1, src2);
}
void my_opencl_acosh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("acosh", 1, dest, src1, src2);
}
void my_opencl_atanh (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("atanh", 1, dest, src1, src2);
}
void my_opencl_cbrt (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("cbrt", 1, dest, src1, src2);
}
void my_opencl_hypot (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("hypot", 2, dest, src1, src2);
}
void my_opencl_tgamma (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("tgamma", 1, dest, src1, src2);
}
void my_opencl_exp2 (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("exp2", 1, dest, src1, src2);
}
void my_opencl_expm1 (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("expm1", 1, dest, src1, src2);
}
void my_opencl_log10 (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("log10", 1, dest, src1, src2);
}
void my_opencl_log2 (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("log2", 1, dest, src1, src2);
}
void my_opencl_log1p (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("log1p", 1, dest, src1, src2);
}
void my_opencl_pow (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("pow", 2, dest, src1, src2);
}
void my_opencl_erf (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("erf", 1, dest, src1, src2);
}
void my_opencl_erfc (void *dest, const void *src1, const void *src2) {
  build_and_run_opencl_function("erfc", 1, dest, src1, src2);
}

static mpcheck_user_func_t tab[] = {
  {"add", my_opencl_add, 0, 0},
  {"add", my_opencl_add, LONG_MAX, LONG_MAX},
  {"sub", my_opencl_sub, 0, 0},
  {"sub", my_opencl_sub, LONG_MAX, LONG_MAX},
  {"mul", my_opencl_mul, 0, 0},
  {"mul", my_opencl_mul, LONG_MAX, LONG_MAX},
  {"div", my_opencl_div, 0, 0},
  {"div", my_opencl_div, LONG_MAX, LONG_MAX},
  {"sqrt", my_opencl_sqrt, 0, 0},
  {"sqrt", my_opencl_sqrt, LONG_MAX, LONG_MAX},
  {"sqrt", my_opencl_sqrt, LONG_MAX-2, 0},
  {"exp", my_opencl_exp, 0, 0},
  {"exp", my_opencl_exp, 9, 0}, 
  {"log", my_opencl_log, 0, 0},
  {"log", my_opencl_log, LONG_MAX, 0},
  {"sin", my_opencl_sin, 0, 0},
  {"sin", my_opencl_sin, 10, 0},
  {"cos", my_opencl_cos, 0, 0},
  {"cos", my_opencl_cos, 10, 0},
  {"tan", my_opencl_tan, 0, 0},
  {"tan", my_opencl_tan, 10, 0},
  {"atan", my_opencl_atan, 0, 0},
  {"atan", my_opencl_atan, 53, 0},
  {"atan2", my_opencl_atan2, 0, 0},
  {"atan2", my_opencl_atan2, 53, 0},
  {"asin", my_opencl_asin, 0, 0},
  {"asin", my_opencl_asin, -10, 0},
  {"acos", my_opencl_acos, 0, 0},
  {"acos", my_opencl_acos, -10, 0},
  {"sinh", my_opencl_sinh, 0, 0},
  {"sinh", my_opencl_sinh, 9, 0}, 
  {"cosh", my_opencl_cosh, 0, 0},
  {"cosh", my_opencl_cosh, 9, 0}, 
  {"tanh", my_opencl_tanh, 0, 0},
  {"tanh", my_opencl_tanh, 4, 0}, 
  {"asinh", my_opencl_asinh, 0, 0},
  {"asinh", my_opencl_asinh, 9, 0}, 
  {"acosh", my_opencl_acosh, 1, 0},
  {"acosh", my_opencl_acosh, 9, 0},
  {"atanh", my_opencl_atanh, 0, 0},
  {"atanh", my_opencl_atanh,-10, 0},
  {"cbrt", my_opencl_cbrt, 0, 0},
  {"cbrt", my_opencl_cbrt, LONG_MAX, 0},
  {"cbrt", my_opencl_cbrt, -1010, 0},
  {"hypot", my_opencl_hypot, 0, 0},
  {"hypot", my_opencl_hypot, LONG_MAX, LONG_MAX},
  {"hypot", my_opencl_hypot, -1010, -1010},
  {"tgamma", my_opencl_tgamma, 0, 0},
  {"tgamma", my_opencl_tgamma, 10, 0},
  {"exp2", my_opencl_exp2, 0, 0},
  {"exp2", my_opencl_exp2, 9, 0},
  {"log2", my_opencl_log2,  0, 0},
  {"log2", my_opencl_log2, LONG_MAX,0},
  {"expm1", my_opencl_expm1, 0, 0},
  {"expm1", my_opencl_expm1, -9, 0},
  {"log10", my_opencl_log10, 0, 0},
  {"log10", my_opencl_log10, LONG_MAX},
  {"log1p", my_opencl_log1p, 0, 0},
  {"log1p", my_opencl_log1p, LONG_MAX, 0},
  {"pow", my_opencl_pow, 0, 0},
  {"pow", my_opencl_pow, 5, 4},
  {"erf", my_opencl_erf, 0, 0},
  {"erf", my_opencl_erf, 9, 0},
  {"erfc", my_opencl_erfc, 0, 0},
  {"erfc", my_opencl_erfc, 2, 0},
  {NULL, NULL, 0, 0}
};


int main (int argc, const char *argv[])
{
  mpfr_t x,y;
  void *fp;
  gmp_randstate_t state;

  /* Init Open CL */
  init_opencl();

  /* Check if interface works */
  gmp_randinit_default (state);
  mpfr_inits2 (prec, x, y, NULL);
  fp = new_fp (prec);
  mpfr_urandomb (x, state);
  mpfr_set (y, x, GMP_RNDN);
  get_fp (fp, x);
  set_fp (y, fp);
  if (mpfr_cmp (x, y) != 0) {
    fprintf (stderr, "ERROR: set_fp(get_fp)) != Identity\n");
    exit (1);
  }

  /* Start test */
  mpcheck_init (argc, argv, prec, emin, emax,
		new_fp, del_fp, get_fp, set_fp, set_rnd_mode,
		ALL_RND, ALL_TEST, 0, 10000, 2);
  mpcheck_check (stdout, tab);
  mpcheck_clear (stdout);
  gmp_randclear (state);

  quit_opencl();
  return 0;
}

#else

#include <stdio.h>
int main () {
  fprintf (stderr, "OPENCL not available.\n");
  return 0;
}

#endif
