/**************************************************************/
/*                                                            */
/*                                                            */
/*     Testing the accuracy of elementary functions           */
/*                                                            */
/*             Projects Arenaire and Spaces                   */
/*                                                            */
/**************************************************************/

#include "mpcheck.h"

int    (*set_fp) (mpfr_ptr, fptype, mp_rnd_t) = NULL;
fptype (*get_fp) (mpfr_srcptr, mp_rnd_t) = NULL;
void   (*fprint_fp) (FILE *, fptype) = NULL;

double MAX_ERR_NEAR = 0.0;
double MAX_ERR_DIR = 0.0;

/* Global options */
int verbose = 3;
int test_monotonicity = 1;
int test_range = 1;
int test_symmetry = 1;

void
testall (unsigned long N, unsigned long seed)
{
  test ("exp",    0, 0, N, seed);
  test ("exp",    9, 0, N, seed);
  test ("exp2",   0, 0, N, seed);
  test ("exp2",   9, 0, N, seed);
  test ("log",    0, 0, N, seed);
  test ("log", EMAX, 0, N, seed);
  test ("log2",   0, 0, N, seed);
  test ("log2", EMAX,0,  N, seed);
  test ("sin",    0, 0, N, seed);
  test ("sin",   10, 0, N, seed); /* mpfr-2.0.1 is too slow for 1024 */
  test ("cos",    0, 0, N, seed);
  test ("cos",   10, 0, N, seed);
  test ("tan",    0, 0, N, seed);
  test ("tan",   10, 0, N, seed);
  test ("asin",   0, 0, N, seed);
  test ("asin", -10, 0, N, seed); /* mpfr-2.0.1 is too slow for -1021 */
  test ("acos",   0, 0, N, seed);
  test ("acos", -10, 0, N, seed); /* mpfr-2.0.1 is too slow for -1021 */
  test ("atan",   0, 0, N, seed);
  test ("atan",  53, 0, N, seed); /* why 53 ??? mpfr too slow ? */
  test ("sqrt",  0, 0, N, seed);
  test ("sqrt",  EMAX, 0, N, seed);
  test ("sqrt",  EMIN, 0, N, seed);
  test ("pow", 0, 0, N, seed);
#if (FPPREC <= 53)
  test ("pow", 8, 7, N, seed);
#else
  test ("pow", 16, 10, N, seed);
#endif
  test ("add", 0, 0, N, seed);
  test ("add", EMAX-1, EMAX-1, N, seed);
  test ("sub", EMAX, EMAX, N, seed);
  test ("sub", 0, 0, N, seed);
  test ("mul", 0, 0, N, seed);
  test ("mul", EMAX/2, EMAX/2, N, seed);
  test ("div", 0, 0, N, seed);
  test ("div", EMAX, EMAX, N, seed);

#ifndef MATHLIB /* mathlib does not implement those functions */
  test ("expm1",  0, 0, N, seed);
  test ("expm1", -9, 0, N, seed);
  test ("log10",  0, 0, N, seed);
  test ("log10", EMAX, 0, N, seed);
  test ("log1p",  0, 0, N, seed);
  test ("log1p", EMAX, 0, N, seed);
  test ("sinh",   0, 0, N, seed);
  test ("sinh",   9, 0, N, seed);
  test ("cosh",   0, 0, N, seed);
  test ("cosh",   9, 0, N, seed);
  test ("tanh",   0, 0, N, seed);
  test ("tanh",   4, 0, N, seed);
  test ("asinh",  0, 0, N, seed);
  test ("asinh", EMAX, 0, N, seed);
  test ("acosh",  1, 0, N, seed);
  test ("acosh", EMAX,0, N, seed);
  test ("atanh",  0, 0, N, seed);
  test ("atanh", -10, 0, N, seed);
  test ("gamma",  0, 0, N, seed);
#if (FPPREC <= 53)
  test ("gamma",  7, 0, N, seed);
#else
  test ("gamma", 10, 0, N, seed);
#endif
  test ("cbrt",  0, 0, N, seed);
  test ("cbrt",  EMAX, 0, N, seed);
  test ("cbrt",  EMIN, 0, N, seed);
  test ("hypot", 0, 0, N, seed);
  test ("hypot", EMAX-1, EMAX-1, N, seed);
  test ("hypot", EMIN, EMIN, N, seed);
#endif

  printf ("Maximal errors for all functions: %f (nearest), %f (directed)\n",
          MAX_ERR_NEAR, MAX_ERR_DIR);
}

void 
usage ()
{
  fprintf (stderr, "Usage: mpcheck [options]\n");
  fprintf (stderr, "Usage: mpcheck [options] <function> <exponent>\n");
  fprintf (stderr, "Usage: mpcheck [options] <function> <exp1> <exp2>\n");
  fprintf (stderr, "where options are:\n");
  fprintf (stderr, "-seed s: set random seed to s [default 1]\n");
  fprintf (stderr, "-verb k: set verbose level to k [default 3]\n");
  fprintf (stderr, "-range : do not test output range\n");
  fprintf (stderr, "-mono  : do not test monotonicity\n");
  fprintf (stderr, "-num  n: N\n");
  exit (1);
}


int
main (int argc, char *argv[])
{
  mp_exp_t exponent = 1, exp2 = 1;
  unsigned long N = 1000, seed = 1;
  int i, expc = 0;
  const char *func = NULL;

  setup ();

  for (i = 1 ; i < argc ; i++)
    {
      if (argv[i][0] == '-')
	{
	  if (strcmp (argv[1], "-mono") == 0)
	    test_monotonicity = 0;
	  else if (strcmp (argv[i], "-rang") == 0)
	    test_range = 0;
	  else if (strcmp (argv[i], "-symm") == 0)
	    test_symmetry = 0;
	  else if (strcmp (argv[i], "-seed") == 0)
	    {
	      seed = atoi (argv[i+1]);
	      i++;
	    }
	  else if (strcmp (argv[i], "-verb") == 0)
	    {
	      verbose = atoi (argv[i+1]);
	      i++;
	    }
	  else if (strcmp (argv[i], "-num") == 0)
	    {
	      N = atoi (argv[i+1]);
	      i++;
	    }
	  else
	    usage ();

	}
      else if (isdigit(argv[i][0]))
	{
	  if (expc == 0)
	    exponent = atol (argv[i]);
	  else if (expc == 1)
	    exp2 = atol (argv[i]);
	  else 
	    usage ();
	  expc ++;
	}
      else
	{
	  if (func != NULL)
	    usage ();
	  func = argv[i];
	}
    }

  printf ("*******************************************************************\n");
  printf ("*                                                                 *\n");
  printf ("* MpCheck version 1.1 (c) INRIA 2002, 2004 (Arenaire & Spaces)    *\n");
  printf ("*                                                                 *\n");
  printf ("*******************************************************************\n");

#ifdef MATHLIB
  printf ("Testing MathLib ");
#else
  printf ("Testing libm ");
#endif

  printf ("[precision=%u, seed=%lu]\n", FPPREC, seed);

  if (func == NULL)
    testall (N, seed);
  else
    test (func, exponent, exp2, N, seed);

  return 0;
}
