from sympy import *
import numpy as np

x,y = symbols('x y', real=True)

def L (k, x):
  if k == 0:
    return 1
  elif k == 1:
    return x
  else:
    _k = sympify(k)
    return (2*_k-1)/_k * x * L(k-1,x) - (_k-1)/_k * L(k-2,x)


def l (k, x):
  if k == 0:
    return (1-x)/2
  elif k == 1:
    return (1+x)/2
  else:
    _k = sympify(k)
    return 1/sqrt(2/(2*_k-1)) * integrate(L(k-1,y),(y,-1,x))

def phi (k, x):
  return l(k+2,x) / ((1-x)/2 * (1+x)/2)


def l_ (k, x):
  return l(k,2*x-1)

def phi_ (k, x):
  return phi(k, 2*x-1)

def evaluate (expr):
  return [N(expr.subs(x,0.0)), N(expr.subs(x,0.25)), N(expr.subs(x,0.5)),
          N(expr.subs(x,0.75)), N(expr.subs(x,1.0))]


# Compute and print the Lobatto kernel functions without the sqrt prefactor
for k in range(0,18):
  _k = sympify(k)
  f = 2/sqrt(2/(2*(_k+2)-1))
  phi_k_x = simplify(phi(k,x)/f)
  print('phi(',k,', x): ',phi_k_x.as_poly(x).all_coeffs(),', f=',f)

print("--------------------")

for k in range(0,18):
  _k = sympify(k)
  f = 2/sqrt(2/(2*(_k+2)-1))
  phi_k_x = simplify(phi_(k,x)/f)
  print('phi_(',k,', x): ',phi_k_x.as_poly(x).all_coeffs())

print("--------------------")


# test the Lobatto polynomials by evaluation
for k in range(2,7):
  phi_k_x = simplify(phi(k,x))
  dphi_k_x = simplify(phi(k,x).diff(x))
  d2phi_k_x = simplify(phi(k,x).diff(x,x))
  print('phi(',k,', x): ',evaluate(phi_k_x))
  print('dphi(',k,', x): ',evaluate(dphi_k_x))
  print('d2phi(',k,', x): ',evaluate(d2phi_k_x))

for k in range(2,7):
  phi_k_x = simplify(phi_(k,x))
  dphi_k_x = simplify(phi_(k,x).diff(x))
  d2phi_k_x = simplify(phi_(k,x).diff(x,x))
  print('phi_(',k,', x): ',evaluate(phi_k_x))
  print('dphi_(',k,', x): ',evaluate(dphi_k_x))
  print('d2phi_(',k,', x): ',evaluate(d2phi_k_x))

  l_k_x = simplify(l_(k,x))
  dl_k_x = simplify(l_(k,x).diff(x))
  d2l_k_x = simplify(l_(k,x).diff(x,x))
  print('l(',k,', x):   ',evaluate(l_k_x))
  print('dl(',k,', x):   ',evaluate(dl_k_x))
  print('d2l(',k,', x):   ',evaluate(d2l_k_x))
