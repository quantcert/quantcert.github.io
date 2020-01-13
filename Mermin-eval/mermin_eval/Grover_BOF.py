#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module aimed to reproduce the work done in the article from Batle et al.
"Global versus local quantum correlations in the Grover search algorithm" where
they could not find proof of non locality during Grover's algorithm. On our end,
we did find proof of non locality, so it would require more work to understand 
why their result is different to ours.

Instead of simulation, their work used the explicit expression of the states, 
hence the presence of functions ``a`` and ``b``. (See article for more details)
"""
from sage.all import *
from mermin_eval import *

n = 6
k = 1
nu = arcsin(sqrt(k/2**n))

def a(_j):
  return sin((2*_j+1)/nu)/sqrt(k)

def b(_j):
  return cos((2*_j+1)/nu)/sqrt(2**n-k)

# this considers that the searched state is |0>, and thus that k=1, so it is not
# very interesting to keep k before, but oh well ...
def phi(_a, _b, _n):
  return vector([_a]+[_b]*(2**_n-1))

def rho(_a, _b, _n):
  return matrix(
    [[_a**2]+[_a*_b]*(2**_n-1)] + \
    [[_a*_b]+[_b**2]*(2**_n-1)]*(2**_n-1))

# results have been computed for iterations up to 2
for j in range(3,7):
  def M_func(_a_a_prime):
    _a_a_prime_packed = [[_a_a_prime[3*_i], _a_a_prime[3*_i+1], _a_a_prime[3*_i+2]] 
                for _i in range(len(_a_a_prime)/3)]
    _a = _a_a_prime_packed[:len(_a_a_prime_packed)/2]
    _a_prime = _a_a_prime_packed[len(_a_a_prime_packed)/2:]
    return M_eval_all(n, _a, _a_prime, rho(a(j), b(j), n))

  r = optimize(M_func, [1]*3*n*2, 5, 10**(-1), 10**2)
  print(CC(r[1]).real())