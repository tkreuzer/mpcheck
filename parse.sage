# parse file libm-test-ulps

# Example:
# T = parse("libm-test-ulps")
def parse(file):
   f = open(file,"r")
   current_function = None
   rnd = None
   T = dict()
   while true:
      s = f.readline()
      if s == '':
         break
      if s[0] == '#' or s == '\n':
         continue
      s = s.split(':')
      if s[0] == 'Function':
         s = s[1][1:].strip('"')
         if s[:4] == 'Real' or s[:9] == 'Imaginary' or 'vlen' in s:
            current_function = None # skip
         elif s[-8:] == 'downward':
            rnd = 'RNDD'
            current_function = s[:-9]
         elif s[-10:] == 'towardzero':
            rnd = 'RNDZ'
            current_function = s[:-11]
         elif s[-6:] == 'upward':
            rnd = 'RNDU'
            current_function = s[:-7]
         else:
            rnd = 'RNDN'
            current_function = s
      elif s[0] in ["float","double","ldouble","float128"]:
         if current_function is not None:
            assert rnd is not None
            type = s[0]
            maxerr = ZZ(s[1])
            key = (current_function,type,rnd)
            assert not T.has_key(key) 
            T[key] = maxerr
   f.close()
   return T

# export T in .h file
# Example:
# T = parse("libm-test-ulps")
# export(T,"ulps.h")
def export(T,file):
   fp = open(file,"w")
   fp.write("typedef struct {\n")
   fp.write("  char fn[16];\n")
   fp.write("  char type[16];\n")
   fp.write("  char rnd[16];\n")
   fp.write("  int err;\n")
   fp.write("} entry_t;\n")
   fp.write("entry_t T[] = {\n")
   for f,t,r in T.keys():
      e = T[(f,t,r)]
      fp.write("{\"" + f + "\", \"" + t + "\", \""+r+"\", " + str(e) + "},\n")
   fp.write("{\"\",\"\",\"\",0},\n")
   fp.write("};\n")
   fp.close()

