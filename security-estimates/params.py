from MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import *

d = 4096                  #choose 4096 for 64-bit modulus
#q = 1073479681            #32-bit modulus
q = 4611686018326724609  #64-bit modulus
l = 1
k = 1

LC = 100
n = LC+8

beta = 2**32 - d -1
beta = 8*d*beta
print(log2(beta), log2(q), beta < q)
msis = MSISParameterSet(d, l+n, k, beta,q, norm="linf")
print(MSIS_summarize_attacks(msis))

#mlwe = MLWEParameterSet(d,l,n+l+k,2,q, distr="binomial")
#print(MLWE_summarize_attacks(mlwe))
