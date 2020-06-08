import numpy as np
#import scipy.stats.stats as stats

#assuming all variables scalar or numpy array as required

def inverse_proportion_normalize( f, l ):          # a(j) = f(j)/l(j)/(sum_k=1^m-1(f(k)/l(k))), j=1,2,...,m-1, Eqn. (1)
  return (f[:-1]/l)/((f[:-1]/l).sum())          

def bootstrap_standard_error( f, l, bts ):         # se(a(j)) = stderr( bts_a[:,j] ), Eqn. (11)
  bts_a = np.array( [ inverse_proportion_normalize( bt, l ) for bt in bts ] )
  return np.std( bts_a, axis=0 )

def weighted_average( a, l ):                      # Eqn. (14)
  return np.dot(a,l)

"""
def asymptotic_standard_error( f, l, Prg ):
  Io = np.array( [[0.]*len(l)]*len(l) )
  for x in range(0, len(l)):
    for y in range(0, len(l)):
      Io[x,y] = 0
      for i in range

"""


