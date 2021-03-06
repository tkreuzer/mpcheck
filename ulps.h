/* this file was produced using the parse.sage script from the
   sysdeps/x86_64/fpu/libm-test-ulps from glibc revision a8f0fc4,
   and edited to match worst cases found by mpcheck:
   * for lines without "was", it is the original glibc bound
   * for lines with "was", the original glibc bound is indicated */
typedef struct {
  char fn[16];
  char type[16];
  char rnd[16];
  long err;
} entry_t;
entry_t max_ulps[] = {
{"yn", "ldouble", "RNDZ", 5},
{"cabs", "double", "RNDU", 1},
{"tgamma", "float128", "RNDZ", 6}, /* was 5 */
{"carg", "double", "RNDD", 1},
{"y1", "double", "RNDN", 3},
{"expm1", "ldouble", "RNDD", 4},
{"j1", "float128", "RNDN", 4},
{"carg", "float", "RNDD", 2},
{"cabs", "ldouble", "RNDZ", 1},
{"tgamma", "float", "RNDZ", 5},
{"acos", "ldouble", "RNDU", 3}, /* was 2 */
{"y1", "ldouble", "RNDU", 7},
{"sin", "ldouble", "RNDZ", 3}, /* was 2 */
{"lgamma", "float128", "RNDD", 8},
{"sinh", "ldouble", "RNDD", 5},
{"log", "float128", "RNDU", 2}, /* was 1 */
{"expm1", "float", "RNDZ", 2},
{"sinh", "float", "RNDU", 3},
{"jn", "double", "RNDZ", 5},
{"acos", "double", "RNDU", 1},
{"tanh", "float", "RNDU", 3},
{"lgamma", "ldouble", "RNDN", 4},
{"cosh", "ldouble", "RNDD", 3},
{"carg", "ldouble", "RNDD", 1},
{"atan2", "float128", "RNDU", 3}, /* was 2 */
{"j1", "double", "RNDD", 3},
{"carg", "float128", "RNDN", 2},
{"tgamma", "ldouble", "RNDU", 6}, /* was 5 */
{"hypot", "ldouble", "RNDU", 1},
{"hypot", "double", "RNDZ", 1},
{"tan", "ldouble", "RNDD", 3},
{"sincos", "float128", "RNDZ", 2},
{"atan2", "float", "RNDU", 3}, /* was 2 */
{"log1p", "float128", "RNDU", 2},
{"asinh", "float", "RNDD", 3},
{"log", "ldouble", "RNDU", 1},
{"erf", "float", "RNDZ", 1},
{"y0", "double", "RNDU", 3},
{"pow", "float", "RNDD", 1},
{"log1p", "double", "RNDN", 1},
{"gamma", "double", "RNDN", 4},
{"atanh", "double", "RNDD", 3},
{"exp2", "double", "RNDZ", 1},
{"exp10", "ldouble", "RNDZ", 2},
{"cabs", "double", "RNDD", 1},
{"asin", "float", "RNDD", 1},
{"y0", "float", "RNDU", 5},
{"yn", "double", "RNDD", 3},
{"log1p", "float", "RNDZ", 2},
{"exp2", "float128", "RNDZ", 1},
{"y0", "float128", "RNDD", 4},
{"cbrt", "ldouble", "RNDD", 2}, /* was 1 */
{"cosh", "float", "RNDN", 1},
{"j0", "double", "RNDN", 2},
{"exp2", "float", "RNDZ", 1},
{"atanh", "float", "RNDU", 3},
{"atanh", "float128", "RNDD", 4},
{"asinh", "float128", "RNDD", 4},
{"log", "float128", "RNDD", 2}, /* was 1 */
{"log10", "float128", "RNDZ", 1},
{"sincos", "double", "RNDZ", 1},
{"exp", "double", "RNDZ", 1},
{"j1", "float", "RNDN", 2},
{"erfc", "float", "RNDZ", 5}, /* was 4 */
{"atan", "double", "RNDZ", 1},
{"tanh", "float128", "RNDN", 2},
{"cosh", "ldouble", "RNDU", 3},
{"atan", "float128", "RNDZ", 2}, /* was 1 */
{"carg", "ldouble", "RNDU", 1},
{"log10", "float128", "RNDD", 2}, /* was 1 */
{"exp10", "float", "RNDZ", 1},
{"asinh", "ldouble", "RNDZ", 5}, /* was 4 */
{"jn", "float", "RNDN", 4},
{"gamma", "double", "RNDD", 5},
{"expm1", "double", "RNDN", 1},
{"log10", "double", "RNDD", 3}, /* was 2 */
{"sinh", "float", "RNDN", 2},
{"erf", "ldouble", "RNDZ", 1},
{"acos", "double", "RNDZ", 1},
{"atan2", "float", "RNDD", 3}, /* was 2 */
{"exp", "ldouble", "RNDD", 2}, /* was 1 */
{"asin", "float128", "RNDD", 2},
{"erfc", "double", "RNDD", 6}, /* was 5 */
{"y1", "float", "RNDU", 2},
{"log1p", "ldouble", "RNDZ", 4},
{"exp10", "float128", "RNDN", 2},
{"j1", "float", "RNDU", 5},
{"cabs", "float128", "RNDD", 1},
{"log10", "ldouble", "RNDZ", 2},
{"gamma", "float", "RNDZ", 4},
{"sin", "float", "RNDZ", 1},
{"log2", "float", "RNDN", 1},
{"cbrt", "double", "RNDU", 5},
{"y0", "float", "RNDD", 4},
{"log2", "ldouble", "RNDN", 1},
{"jn", "float128", "RNDN", 7},
{"acosh", "double", "RNDU", 4}, /* was 2 */
{"tgamma", "float128", "RNDU", 6}, /* was 4 */
{"cbrt", "float128", "RNDZ", 1},
{"acosh", "float128", "RNDN", 4}, /* was 2 */
{"asinh", "float128", "RNDU", 4},
{"tan", "float128", "RNDZ", 2}, /* was 1 */
{"y1", "ldouble", "RNDZ", 5},
{"expm1", "float", "RNDD", 1},
{"acosh", "ldouble", "RNDZ", 5}, /* was 4 */
{"lgamma", "double", "RNDN", 4},
{"carg", "ldouble", "RNDZ", 1},
{"exp", "float128", "RNDN", 1},
{"exp", "float128", "RNDZ", 1}, /* was unspecified */
{"exp", "float128", "RNDU", 0}, /* was unspecified */
{"exp", "float128", "RNDD", 1}, /* was unspecified */
{"log", "ldouble", "RNDZ", 3}, /* was 2 */
{"tanh", "ldouble", "RNDN", 3},
{"sin", "float128", "RNDU", 3},
{"tan", "float", "RNDN", 1},
{"acos", "ldouble", "RNDD", 3}, /* was 2 */
{"exp10", "ldouble", "RNDN", 1},
{"sinh", "ldouble", "RNDN", 2},
{"expm1", "float128", "RNDD", 2},
{"tgamma", "ldouble", "RNDZ", 6}, /* was 5 */
{"acosh", "float128", "RNDU", 4}, /* was 2 */
{"y1", "float", "RNDZ", 2},
{"exp10", "double", "RNDN", 2},
{"sincos", "float128", "RNDU", 3},
{"exp", "float", "RNDZ", 1},
{"tanh", "float", "RNDD", 3},
{"j0", "ldouble", "RNDN", 2},
{"tgamma", "double", "RNDN", 6}, /* was 5 */
{"atan", "float", "RNDZ", 1},
{"cos", "float128", "RNDN", 1},
{"acosh", "double", "RNDD", 4}, /* was 2 */
{"exp2", "double", "RNDU", 1},
{"asinh", "double", "RNDN", 2}, /* was 1 */
{"tgamma", "float128", "RNDD", 6}, /* was 5 */
{"atan2", "ldouble", "RNDN", 1},
{"expm1", "ldouble", "RNDN", 2},
{"exp2", "float128", "RNDU", 2},
{"yn", "float128", "RNDD", 5},
{"asinh", "float", "RNDU", 3},
{"exp10", "float", "RNDU", 1},
{"exp10", "ldouble", "RNDD", 2},
{"asin", "float", "RNDZ", 1},
{"sincos", "double", "RNDU", 1},
{"log2", "ldouble", "RNDU", 1},
{"sinh", "double", "RNDD", 4}, /* was 3 */
{"log2", "double", "RNDD", 3},
{"tan", "float128", "RNDU", 3}, /* was 1 */
{"hypot", "ldouble", "RNDN", 1},
{"erf", "float128", "RNDN", 1},
{"tanh", "float128", "RNDU", 3},
{"sincos", "float", "RNDU", 1},
{"atanh", "float128", "RNDZ", 3}, /* was 2 */
{"cbrt", "ldouble", "RNDU", 1},
{"erfc", "ldouble", "RNDN", 3},
{"cos", "float", "RNDZ", 1},
{"atanh", "float", "RNDD", 3},
{"asinh", "ldouble", "RNDU", 5},
{"expm1", "float128", "RNDU", 3},
{"acos", "float", "RNDZ", 1},
{"yn", "double", "RNDZ", 3},
{"erf", "ldouble", "RNDU", 2}, /* was 1 */
{"asin", "ldouble", "RNDN", 1},
{"tanh", "ldouble", "RNDD", 4},
{"log10", "float", "RNDN", 2},
{"pow", "float", "RNDZ", 1},
{"y1", "ldouble", "RNDD", 7},
{"cbrt", "float", "RNDU", 1},
{"sincos", "ldouble", "RNDN", 1},
{"lgamma", "float128", "RNDN", 5},
{"cabs", "double", "RNDN", 1},
{"sin", "double", "RNDD", 1},
{"log10", "double", "RNDU", 3}, /* was 2 */
{"exp2", "float128", "RNDD", 1},
{"atan", "ldouble", "RNDZ", 1},
{"exp10", "double", "RNDU", 4}, /* was 2 */
{"cbrt", "float128", "RNDU", 1},
{"carg", "float", "RNDZ", 2},
{"erf", "float", "RNDD", 2}, /* was 1 */
{"cosh", "ldouble", "RNDN", 2},
{"asinh", "float128", "RNDZ", 3}, /* was 2 */
{"log", "float128", "RNDN", 1},
{"erfc", "ldouble", "RNDD", 5}, /* was 4 */
{"log2", "float128", "RNDZ", 3}, /* was 1 */
{"sincos", "double", "RNDD", 1},
{"j1", "double", "RNDN", 1},
{"log2", "double", "RNDU", 3},
{"j1", "float", "RNDD", 3},
{"jn", "ldouble", "RNDZ", 5},
{"erfc", "float128", "RNDZ", 4},
{"cbrt", "double", "RNDD", 4},
{"sin", "double", "RNDN", 1},
{"cbrt", "float", "RNDN", 1},
{"tanh", "float128", "RNDD", 4},
{"sincos", "float", "RNDD", 1},
{"tan", "ldouble", "RNDN", 2},
{"asinh", "float", "RNDN", 2}, /* was 1 */
{"expm1", "double", "RNDD", 1},
{"acos", "float128", "RNDN", 1},
{"sinh", "double", "RNDU", 3},
{"exp10", "float128", "RNDU", 3},
{"atanh", "double", "RNDZ", 3}, /* was 2 */
{"expm1", "float", "RNDU", 1},
{"pow", "ldouble", "RNDN", 1},
{"erf", "ldouble", "RNDD", 2}, /* was 1 */
{"j1", "ldouble", "RNDU", 3},
{"exp", "ldouble", "RNDN", 1},
{"tan", "double", "RNDD", 1},
{"jn", "ldouble", "RNDD", 4},
{"erfc", "double", "RNDN", 3},
{"y0", "double", "RNDN", 2},
{"cos", "double", "RNDN", 1},
{"j0", "ldouble", "RNDU", 4},
{"cabs", "float128", "RNDN", 1},
{"j1", "ldouble", "RNDD", 4},
{"yn", "ldouble", "RNDD", 5},
{"tgamma", "float", "RNDD", 5},
{"jn", "float128", "RNDD", 8},
{"pow", "double", "RNDN", 1},
{"y1", "float128", "RNDN", 2},
{"cbrt", "float128", "RNDD", 1},
{"exp", "float", "RNDD", 1},
{"acosh", "float128", "RNDD", 5}, /* was 3 */
{"sinh", "float128", "RNDN", 2},
{"erf", "float", "RNDU", 2}, /* was 1 */
{"jn", "float128", "RNDU", 7},
{"sincos", "float128", "RNDD", 3},
{"y1", "float128", "RNDU", 5},
{"erfc", "ldouble", "RNDZ", 4},
{"exp2", "ldouble", "RNDN", 1},
{"atan2", "float128", "RNDZ", 3},
{"log10", "double", "RNDN", 2},
{"lgamma", "float", "RNDN", 4},
{"exp2", "double", "RNDD", 1},
{"atan2", "float", "RNDZ", 3}, /* was 2 */
{"cos", "ldouble", "RNDZ", 3}, /* was 2 */
{"asin", "float128", "RNDN", 1},
{"atanh", "ldouble", "RNDZ", 5}, /* was 4 */
{"cos", "float", "RNDD", 1},
{"cosh", "double", "RNDN", 1},
{"cbrt", "ldouble", "RNDZ", 1},
{"log10", "float128", "RNDU", 2}, /* was 1 */
{"gamma", "double", "RNDU", 5},
{"exp10", "float", "RNDD", 1},
{"cosh", "float128", "RNDN", 1},
{"expm1", "float128", "RNDZ", 4},
{"log2", "float", "RNDD", 3},
{"log2", "ldouble", "RNDD", 1},
{"atan", "float128", "RNDN", 1},
{"atan", "ldouble", "RNDU", 1},
{"j0", "ldouble", "RNDD", 4},
{"cos", "float128", "RNDD", 3},
{"expm1", "float", "RNDN", 1},
{"erf", "double", "RNDN", 1},
{"sin", "float128", "RNDZ", 3}, /* was 2 */
{"yn", "ldouble", "RNDN", 4},
{"asinh", "ldouble", "RNDD", 5},
{"asinh", "double", "RNDD", 3},
{"tgamma", "float128", "RNDN", 6}, /* was 4 */
{"yn", "float", "RNDN", 3},
{"yn", "float128", "RNDU", 5},
{"j0", "float", "RNDN", 2},
{"y0", "float128", "RNDZ", 3},
{"yn", "ldouble", "RNDU", 4},
{"erfc", "ldouble", "RNDU", 6}, /* was 5 */
{"acos", "float", "RNDD", 1},
{"acos", "ldouble", "RNDN", 2}, /* was 1 */
{"atan", "float", "RNDU", 2},
{"log2", "double", "RNDN", 2},
{"hypot", "double", "RNDN", 1},
{"gamma", "float", "RNDD", 4},
{"y0", "ldouble", "RNDN", 1},
{"y0", "float", "RNDZ", 3},
{"erf", "float128", "RNDD", 2},
{"pow", "double", "RNDZ", 1},
{"exp10", "double", "RNDD", 4}, /* was 3 */
{"tanh", "float", "RNDN", 2},
{"gamma", "double", "RNDZ", 5},
{"carg", "float128", "RNDZ", 3},
{"atan", "float", "RNDD", 2},
{"log2", "float", "RNDU", 3},
{"log1p", "float", "RNDN", 1},
{"acosh", "float", "RNDZ", 3}, /* was 2 */
{"atan", "ldouble", "RNDD", 1},
{"acosh", "ldouble", "RNDD", 5}, /* was 4 */
{"atan2", "ldouble", "RNDD", 1},
{"sinh", "double", "RNDZ", 4}, /* was 2 */
{"log10", "float", "RNDD", 3},
{"hypot", "float128", "RNDN", 1},
{"sinh", "ldouble", "RNDZ", 4},
{"j0", "float128", "RNDZ", 2},
{"cos", "float128", "RNDU", 3}, /* was 2 */
{"log2", "float128", "RNDD", 3},
{"sincos", "ldouble", "RNDD", 3},
{"j1", "float", "RNDZ", 2},
{"erfc", "float", "RNDN", 3}, /* was 2 */
{"cosh", "float128", "RNDU", 3},
{"expm1", "double", "RNDU", 1},
{"yn", "float", "RNDZ", 3},
{"acosh", "double", "RNDZ", 4}, /* was 2 */
{"exp2", "float128", "RNDN", 1},
{"log", "float", "RNDN", 1},
{"gamma", "ldouble", "RNDZ", 7},
{"tan", "double", "RNDU", 1},
{"y1", "double", "RNDZ", 3},
{"j0", "double", "RNDZ", 3},
{"cabs", "ldouble", "RNDN", 1},
{"log1p", "float", "RNDU", 2},
{"tanh", "double", "RNDN", 2},
{"atan2", "double", "RNDZ", 1},
{"tgamma", "float", "RNDN", 5},
{"sincos", "double", "RNDN", 1},
{"j1", "ldouble", "RNDN", 1},
{"atanh", "float", "RNDN", 2},
{"gamma", "float", "RNDU", 5},
{"cbrt", "double", "RNDN", 3},
{"cbrt", "float", "RNDD", 1},
{"jn", "double", "RNDN", 4},
{"hypot", "float128", "RNDU", 1},
{"log1p", "float128", "RNDD", 3},
{"exp", "float", "RNDU", 1},
{"cos", "ldouble", "RNDU", 3}, /* was 2 */
{"asin", "ldouble", "RNDD", 2},
{"j1", "float128", "RNDZ", 4},
{"carg", "float128", "RNDU", 2},
{"jn", "float", "RNDZ", 5},
{"atan2", "float128", "RNDD", 3}, /* was 2 */
{"tan", "float", "RNDZ", 2}, /* was 1 */
{"sin", "ldouble", "RNDN", 1},
{"pow", "ldouble", "RNDD", 4},
{"acosh", "float", "RNDU", 3}, /* was 2 */
{"erf", "ldouble", "RNDN", 1},
{"acosh", "ldouble", "RNDU", 5}, /* was 3 */
{"exp10", "ldouble", "RNDU", 2},
{"acosh", "float128", "RNDZ", 4}, /* was 2 */
{"tanh", "double", "RNDU", 3},
{"cos", "double", "RNDD", 1},
{"sincos", "ldouble", "RNDU", 3},
{"erf", "float", "RNDN", 1},
{"cos", "float", "RNDU", 1},
{"log2", "float", "RNDZ", 2},
{"asin", "double", "RNDZ", 1},
{"tgamma", "float", "RNDU", 5},
{"yn", "float", "RNDU", 5},
{"log2", "ldouble", "RNDZ", 1},
{"carg", "float", "RNDU", 2},
{"lgamma", "double", "RNDZ", 6}, /* was 5 */
{"gamma", "ldouble", "RNDU", 6},
{"cbrt", "float128", "RNDN", 1},
{"lgamma", "float", "RNDU", 5},
{"atanh", "ldouble", "RNDD", 5},
{"sinh", "float128", "RNDD", 3},
{"sincos", "float128", "RNDN", 1},
{"exp2", "ldouble", "RNDD", 1},
{"exp2", "float", "RNDN", 1},
{"j1", "ldouble", "RNDZ", 4},
{"cosh", "double", "RNDU", 2},
{"acos", "float128", "RNDD", 2}, /* was 1 */
{"tgamma", "ldouble", "RNDN", 6}, /* was 5 */
{"lgamma", "ldouble", "RNDU", 6},
{"hypot", "double", "RNDU", 1},
{"exp2", "double", "RNDN", 1},
{"log", "float", "RNDU", 2},
{"asinh", "double", "RNDU", 3},
{"exp10", "double", "RNDZ", 4}, /* was 3 */
{"erfc", "float128", "RNDD", 5},
{"tan", "double", "RNDZ", 1},
{"asin", "ldouble", "RNDU", 1},
{"tanh", "float128", "RNDZ", 3},
{"j1", "float128", "RNDU", 3},
{"log", "ldouble", "RNDN", 1},
{"tgamma", "double", "RNDZ", 6}, /* was 5 */
{"sinh", "ldouble", "RNDU", 5},
{"atan2", "double", "RNDU", 1},
{"cosh", "float", "RNDZ", 2}, /* was 1 */
{"pow", "float128", "RNDZ", 2},
{"acosh", "float", "RNDD", 3}, /* was 2 */
{"acos", "float", "RNDU", 1},
{"expm1", "ldouble", "RNDZ", 4},
{"sinh", "float", "RNDZ", 3}, /* was 2 */
{"y1", "float", "RNDN", 2},
{"j0", "float128", "RNDD", 4},
{"erf", "float128", "RNDU", 2},
{"sincos", "float", "RNDZ", 1},
{"exp2", "float", "RNDU", 1},
{"log1p", "ldouble", "RNDN", 2},
{"cos", "float128", "RNDZ", 3}, /* was 1 */
{"exp10", "float128", "RNDZ", 3},
{"erf", "double", "RNDD", 2}, /* was 1 */
{"asinh", "ldouble", "RNDN", 3},
{"expm1", "double", "RNDZ", 1},
{"gamma", "ldouble", "RNDD", 7},
{"lgamma", "float", "RNDD", 4},
{"sinh", "float128", "RNDU", 4},
{"atan2", "ldouble", "RNDU", 1},
{"j0", "float", "RNDD", 4},
{"log10", "float", "RNDU", 2},
{"cosh", "double", "RNDD", 2},
{"asin", "double", "RNDU", 1},
{"acos", "float128", "RNDU", 1},
{"acos", "float", "RNDN", 1},
{"log2", "float128", "RNDU", 2}, /* was 1 */
{"tan", "float128", "RNDN", 1},
{"pow", "float128", "RNDN", 2},
{"hypot", "double", "RNDD", 1},
{"sin", "float128", "RNDD", 3},
{"cosh", "float128", "RNDD", 2},
{"gamma", "float", "RNDN", 4},
{"lgamma", "double", "RNDU", 5},
{"y0", "ldouble", "RNDD", 5},
{"erfc", "float128", "RNDU", 5},
{"hypot", "ldouble", "RNDZ", 1},
{"jn", "float128", "RNDZ", 8},
{"j1", "float128", "RNDD", 4},
{"log1p", "ldouble", "RNDU", 3},
{"tanh", "ldouble", "RNDZ", 4}, /* was 3 */
{"y1", "ldouble", "RNDN", 2},
{"lgamma", "ldouble", "RNDZ", 7},
{"log1p", "float", "RNDD", 2},
{"lgamma", "ldouble", "RNDD", 7},
{"atan", "ldouble", "RNDN", 1},
{"acosh", "ldouble", "RNDN", 3}, /* was 2 */
{"asin", "ldouble", "RNDZ", 1},
{"tanh", "float", "RNDZ", 3}, /* was 2 */
{"hypot", "float128", "RNDD", 1},
{"cosh", "ldouble", "RNDZ", 2},
{"y0", "ldouble", "RNDU", 3},
{"cos", "ldouble", "RNDD", 3},
{"carg", "ldouble", "RNDN", 1},
{"log2", "float128", "RNDN", 3}, /* was 2 */
{"sincos", "ldouble", "RNDZ", 2},
{"erf", "double", "RNDU", 2}, /* was 1 */
{"carg", "float128", "RNDD", 2},
{"sin", "float", "RNDU", 2}, /* was 1 */
{"pow", "ldouble", "RNDU", 4},
{"log10", "ldouble", "RNDD", 2},
{"tan", "ldouble", "RNDZ", 3},
{"y1", "double", "RNDD", 3},
{"j0", "float", "RNDU", 2},
{"carg", "float", "RNDN", 1},
{"cabs", "ldouble", "RNDD", 1},
{"log1p", "double", "RNDD", 2},
{"tanh", "double", "RNDD", 3},
{"asin", "double", "RNDD", 1},
{"cos", "double", "RNDU", 1},
{"atanh", "double", "RNDN", 2},
{"cosh", "float", "RNDU", 2},
{"j1", "double", "RNDZ", 3},
{"exp", "double", "RNDD", 1},
{"j0", "ldouble", "RNDZ", 5},
{"jn", "ldouble", "RNDN", 4},
{"atan", "float", "RNDN", 1},
{"log1p", "double", "RNDU", 2},
{"asin", "float", "RNDN", 1},
{"log", "double", "RNDU", 1},
{"cbrt", "float", "RNDZ", 1},
{"jn", "double", "RNDD", 5},
{"asinh", "double", "RNDZ", 2},
{"atan", "float128", "RNDD", 2},
{"atan2", "ldouble", "RNDZ", 1},
{"atanh", "ldouble", "RNDU", 5},
{"j0", "double", "RNDD", 2},
{"tgamma", "double", "RNDU", 6}, /* was 5 */
{"cabs", "float128", "RNDZ", 1},
{"atanh", "float128", "RNDN", 3},
{"acos", "float128", "RNDZ", 2}, /* was 1 */
{"pow", "float128", "RNDU", 2},
{"exp2", "ldouble", "RNDU", 1},
{"erfc", "float", "RNDD", 6},
{"y1", "float128", "RNDZ", 2},
{"log10", "float128", "RNDN", 2}, /* was 1 */
{"exp", "ldouble", "RNDZ", 2},
{"expm1", "ldouble", "RNDU", 4},
{"sinh", "float128", "RNDZ", 3},
{"log", "float", "RNDD", 2},
{"erf", "float128", "RNDZ", 1},
{"pow", "double", "RNDD", 1},
{"y0", "double", "RNDD", 3},
{"pow", "float", "RNDN", 1},
{"jn", "float", "RNDD", 5},
{"sin", "float", "RNDD", 2}, /* was 1 */
{"atan2", "double", "RNDD", 1},
{"sinh", "float", "RNDD", 3},
{"lgamma", "float128", "RNDZ", 6}, /* was 5 */
{"cabs", "double", "RNDZ", 1},
{"acos", "double", "RNDD", 1},
{"atan2", "float", "RNDN", 1},
{"yn", "double", "RNDN", 3},
{"asin", "float128", "RNDZ", 1},
{"atanh", "ldouble", "RNDN", 3},
{"yn", "float128", "RNDZ", 5},
{"cabs", "ldouble", "RNDU", 1},
{"log1p", "float128", "RNDN", 2},
{"j0", "float128", "RNDU", 5},
{"cbrt", "ldouble", "RNDN", 1},
{"carg", "double", "RNDZ", 1},
{"cosh", "float", "RNDD", 2}, /* was 1 */
{"j1", "double", "RNDU", 3},
{"log10", "float", "RNDZ", 2},
{"exp2", "float", "RNDD", 1},
{"log1p", "double", "RNDZ", 2},
{"cosh", "float128", "RNDZ", 2},
{"asinh", "float128", "RNDN", 3},
{"atan2", "float128", "RNDN", 2}, /* was 1 */
{"sin", "double", "RNDU", 1},
{"yn", "float128", "RNDN", 5},
{"tgamma", "ldouble", "RNDD", 6}, /* was 5 */
{"sin", "ldouble", "RNDD", 3},
{"jn", "double", "RNDU", 5},
{"tan", "ldouble", "RNDU", 3}, /* was 2 */
{"erfc", "float128", "RNDN", 3}, /* was 2 */
{"j0", "double", "RNDU", 3},
{"log", "ldouble", "RNDD", 2},
{"tanh", "ldouble", "RNDU", 4},
{"pow", "float", "RNDU", 1},
{"y0", "float128", "RNDN", 3},
{"log1p", "float128", "RNDZ", 3},
{"pow", "float128", "RNDD", 2},
{"log", "float128", "RNDZ", 2},
{"tan", "float", "RNDU", 3}, /* was 1 */
{"asin", "float", "RNDU", 1},
{"yn", "double", "RNDU", 4},
{"log10", "ldouble", "RNDN", 1},
{"lgamma", "double", "RNDD", 5},
{"y0", "float128", "RNDU", 3},
{"sin", "double", "RNDZ", 1},
{"atan", "double", "RNDD", 1},
{"pow", "double", "RNDU", 1},
{"log1p", "ldouble", "RNDD", 4},
{"atanh", "float128", "RNDU", 4},
{"jn", "float", "RNDU", 5},
{"acos", "ldouble", "RNDZ", 3}, /* was 2 */
{"atan", "float128", "RNDU", 2},
{"asinh", "float", "RNDZ", 2},
{"erf", "double", "RNDZ", 1},
{"y0", "float", "RNDN", 1},
{"jn", "ldouble", "RNDU", 5},
{"pow", "ldouble", "RNDZ", 4},
{"erfc", "double", "RNDU", 7}, /* was 5 */
{"carg", "double", "RNDU", 1},
{"j0", "float", "RNDZ", 2},
{"erfc", "double", "RNDZ", 5}, /* was 3 */
{"yn", "float", "RNDD", 4},
{"cos", "double", "RNDZ", 1},
{"sin", "ldouble", "RNDU", 3},
{"log2", "double", "RNDZ", 2},
{"lgamma", "float128", "RNDU", 8},
{"acosh", "float", "RNDN", 2},
{"log10", "ldouble", "RNDU", 2}, /* was 1 */
{"y0", "ldouble", "RNDZ", 5},
{"exp", "ldouble", "RNDU", 2}, /* was 1 */
{"y0", "double", "RNDZ", 3},
{"asin", "float128", "RNDU", 2},
{"y1", "float", "RNDD", 2},
{"y1", "float128", "RNDD", 4},
{"j0", "float128", "RNDN", 2},
{"exp10", "float128", "RNDD", 3},
{"y1", "double", "RNDU", 7},
{"cabs", "float128", "RNDU", 1},
{"cosh", "double", "RNDZ", 2},
{"tan", "float", "RNDD", 3}, /* was 2 */
{"exp", "double", "RNDU", 1},
{"exp2", "ldouble", "RNDZ", 1},
{"log", "float", "RNDZ", 2},
{"expm1", "float128", "RNDN", 1},
{"gamma", "ldouble", "RNDN", 4},
{"tan", "float128", "RNDD", 3}, /* was 1 */
{"hypot", "ldouble", "RNDD", 1},
{"atan", "double", "RNDU", 1},
{"log10", "double", "RNDZ", 3}, /* was 2 */
{"lgamma", "float", "RNDZ", 4},
{"sinh", "double", "RNDN", 2},
{"hypot", "float128", "RNDZ", 1},
{"tanh", "double", "RNDZ", 3}, /* was 2 */
{"cos", "ldouble", "RNDN", 1},
{"tgamma", "double", "RNDD", 6}, /* was 5 */
{"atanh", "float", "RNDZ", 2},
{"cbrt", "double", "RNDZ", 4}, /* was 3 */
{"atanh", "double", "RNDU", 3},
{"sin", "float128", "RNDN", 1},
{"acosh", "double", "RNDN", 2},
{"erfc", "float", "RNDU", 6},
{"","","",0},
};
