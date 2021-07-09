from sympy import *
import numpy as np

x,y = symbols('x y', real=True)

def L (k, x):
  if k == 0:
    return 1
  elif k == 1:
    return x
  else:
    return (2*k-1)/k * x * L(k-1,x) - (k-1)/k * L(k-2,x)


def l (k, x):
  if k == 0:
    return (1-x)/2
  elif k == 1:
    return (1+x)/2
  else:
    return integrate(L(k-1,y),(y,-1,x))

def phi (k, x):
  return l(k+2,x) / l(0,x) / l(1,x)


# transform coordinates in [0,1] to [-1,1]
def x_ (x):
  return 2*x-1

def l_ (k, x):
  return l(k,x_(x))

def phi_ (k, x):
  return l_(k+2,x).as_poly(x) / l_(0,x).as_poly(x) / l_(1,x).as_poly(x)


# Compute and print the lobatto kernel functions without the sqrt prefactor
for k in range(0,4):
  K = sympify(k)
  #f = 2/sqrt(2/(2*(K+2)-1))

  l_k = simplify(phi_(K,x)/2).as_poly(x)
  print('l(',k,', x): ',l_k)
  # print('l(',k,', x): ',l_k.all_coeffs())


# test the lobatto polynomials by evaluation
for k in range(2,4):
  K = sympify(k)
  f = 2/sqrt(2/(2*(K+2)-1))
  l_k = simplify(f*l_(K,x)/2).as_poly(x)

  l_k_0 = N(l_k.subs(x,0.0))
  l_k_1 = N(l_k.subs(x,0.25))
  l_k_2 = N(l_k.subs(x,0.5))
  l_k_3 = N(l_k.subs(x,0.75))
  l_k_4 = N(l_k.subs(x,1.0))

  print('l(',k,', x): ',[l_k_0,l_k_1,l_k_2,l_k_3,l_k_4])