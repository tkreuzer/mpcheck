bin_PROGRAMS = mpcheck-double mpcheck-float mpcheck-longdouble mpcheck-crlibm mpcheck-mathlib mpcheck-libmcr mpcheck-pari mpcheck-ntl mpcheck-cln mpcheck-float128 mpcheck-opencl mpcheck-opencl-double mpcheck-libbf

mpcheck_double_SOURCES = check.c func.c mpcheck.h mpcheck-double.c
mpcheck_float_SOURCES = check.c func.c mpcheck.h mpcheck-float.c
mpcheck_longdouble_SOURCES = check.c func.c mpcheck.h mpcheck-longdouble.c
mpcheck_crlibm_SOURCES = check.c func.c mpcheck.h mpcheck-crlibm.c
mpcheck_mathlib_SOURCES = check.c func.c mpcheck.h mpcheck-mathlib.c
mpcheck_libmcr_SOURCES = check.c func.c mpcheck.h mpcheck-libmcr.c
mpcheck_pari_SOURCES = check.c func.c mpcheck.h mpcheck-pari.c
mpcheck_libbf_SOURCES = check.c func.c mpcheck.h mpcheck-libbf.c
mpcheck_ntl_SOURCES = check.c func.c mpcheck.h mpcheck-ntl.cc
mpcheck_cln_SOURCES = check.c func.c mpcheck.h mpcheck-cln.cc
mpcheck_float128_SOURCES = check.c func.c mpcheck.h mpcheck-float128.c
mpcheck_opencl_SOURCES = check.c func.c mpcheck.h mpcheck-opencl.c
mpcheck_opencl_double_SOURCES = check.c func.c mpcheck.h mpcheck-opencl-double.c

EXTRA_DIST = VERSION rnd_mode.c native.c ulps.h

check: mpcheck-double mpcheck-float mpcheck-longdouble
	./mpcheck-double
	./mpcheck-float
	./mpcheck-longdouble

