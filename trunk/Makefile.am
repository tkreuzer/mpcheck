bin_PROGRAMS = mpcheck-double mpcheck-float mpcheck-longdouble mpcheck-crlibm

mpcheck_double_SOURCES = check.c func.c mpcheck.h mpcheck-double.c
mpcheck_float_SOURCES = check.c func.c mpcheck.h mpcheck-float.c
mpcheck_longdouble_SOURCES = check.c func.c mpcheck.h mpcheck-longdouble.c
mpcheck_crlibm_SOURCES = check.c func.c mpcheck.h mpcheck-crlibm.c

check: mpcheck-double mpcheck-float mpcheck-longdouble
	./mpcheck-double
	./mpcheck-float
	./mpcheck-longdouble

