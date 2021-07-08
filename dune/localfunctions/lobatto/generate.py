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
    return 1/sqrt(2/(2*k-1)) * integrate(L(k-1,y),(y,-1,x))

def phi (k, x):
  return l(k+2,x) / l(0,x) / l(1,x)


# transform coordinates in [0,1] to [-1,1]
def x_ (x):
  return 2*x-1

def l_ (k, x):
  return l(k,x_(x))

def phi_ (k, x):
  return l_(k+2,x).as_poly(x) / l_(0,x).as_poly(x) / l_(1,x).as_poly(x)

factors = [sympify('sqrt(6)'),sympify('sqrt(10)'),sympify('sqrt(14)'),
           sympify('sqrt(2)'),sympify('sqrt(22)'),sympify('sqrt(26)'),
           sympify('sqrt(30)'),sympify('sqrt(34)'),sympify('sqrt(38)'),
           sympify('sqrt(42)'),sympify('sqrt(46)'),sympify('sqrt(2)'),
           sympify('sqrt(6)'),sympify('sqrt(58)'),sympify('sqrt(62)'),
           sympify('sqrt(66)'),sympify('sqrt(70)'),sympify('sqrt(74)'),
           sympify('sqrt(78)'),sympify('sqrt(82)'),sympify('sqrt(86)')]

for k in range(21):
  f = factors[k]
  l_k = simplify(phi_(sympify(k),x)/f).as_poly(x)
  print('l(',k,', x): ',l_k.all_coeffs())
