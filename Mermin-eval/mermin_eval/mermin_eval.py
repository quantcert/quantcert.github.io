#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is a SageMath module aimed at computing Mermin operators optimized 
to detect a specific quantum state.
"""
from sage.all import *

import csv

from opti import *
from run_circuit import *
import grover

def mermin_operator_opti(target_state_vector, precomputed_filename=None, verbose=False):
  r"""
  Computes the pseudo optimal operator used to perform the Mermin evaluation 
  during Grover's algorithm.
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm (only 
    single item searches are supported for now).
  :param str precomputed_filename: File where the precomputed coefficients for 
    optimal Mermin operator are stored, if left empty,  precomputation will
    not be used. If precomputation is used and the searched state is not in 
    this database, once the optimization done, the result will be added to 
    the file.
  :param bool verbose: If *verbose* then extra run information will be displayed in 
    terminal.
  :returns: matrix, real -- The Mermin operator satisfying the required conditions
    and the value reached.
  """
  N = len(target_state_vector)
  n = log(N)/log(2)
  if precomputed_filename:
    target_state_string = grover.target_state_ket_vector_to_string(target_state_vector)
    if not os.path.exists(os.path.dirname(precomputed_filename)) and os.path.dirname(precomputed_filename) != '':
      os.makedirs(os.path.dirname(precomputed_filename))
    try:
      precomputed_file = open(precomputed_filename, "a+b")
    except IOError:
      if verbose:
        print "Could not open file, falling back to default method"
      (a,b,c,m,p,q), _ = mermin_coef_opti(target_state_vector, verbose)
    with precomputed_file:
      reader = csv.DictReader(precomputed_file, delimiter=',', quotechar="\'")
      found = False
      for row in reader:
        if row['ket'] == target_state_string:
          found = True
          a,b,c,m,p,q = RR(row['a']), RR(row['b']), RR(row['c']), \
                        RR(row['m']), RR(row['p']), RR(row['q']), 
      if not found:
        (a,b,c,m,p,q), _ = mermin_coef_opti(target_state_vector, verbose)
        writer = csv.DictWriter(precomputed_file, fieldnames=reader.fieldnames)
        writer.writerow({'ket':target_state_string, 'a':a, 'b':b, 'c':c, 'm':m, 'p':p, 'q':q})
  else:
    (a,b,c,m,p,q), _ = mermin_coef_opti(target_state_vector, verbose)

  return M_from_coef(n,a,b,c,m,p,q)


def mermin_coef_opti(target_state, verbose=False):
  r"""
  Returns the Mermin operator maximizing the measure for a given input.

  :param vector[int] target_state_vector: State searched by Grover's algorithm (only 
    single item searches are supported for now).
  :param bool verbose: If *verbose* then extra run information will be displayed in 
    terminal.
  :returns: matrix, real -- Coefficients of the optimal Mermin operator for 
    *target_state* in th Grover algorithm and the value reached.
  """
  n = log(len(target_state))/log(2)
  plus = vector([1,1])/sqrt(2)
  plus_n = kronecker_power(plus, n)
  phi = (target_state + plus_n).normalized()

  def M_phi(_abcmpq):
    a,b,c,m,p,q = _abcmpq
    return M_eval(a,b,c,m,p,q, phi)

  abcmpq,v = optimize(M_phi, (1,1,1,1,1,1), 5, 10**(-2), 10**2)
  if verbose:
    print("Optimization completed result: optimum reached at (a,b,c,m," + \
      "p,q) = " + str(abcmpq) + " with a value of v = " + \
      str(CC(v).real()))
  return abcmpq,v


def M_eval(a, b, c, m, p, q, phi):
  r"""
  This function evaluates <*phi*|M_n|*phi*> with *(a,b,c,m,p,q)* describing 
  M_n, the Mermin operator.

    M_n traditionally uses two families of operators, a_n and a'_n, in our 
    case, a_n = a*X+b*Y+c*Z and a'_n = m*X+p*Y+q*Z.

  :param real a,b,c,m,p,q: Coefficients for the Mermin operator, used as 
    described above.
  :param vector[complex] phi: Vector to be evaluated with M.
  :returns: complex -- <*phi*|M_n|*phi*>
  """
  phi_matrix = matrix(phi).transpose()
  n = log(len(phi))/log(2)
  return abs((phi_matrix.transpose().conjugate() * M_from_coef(n,a,b,c,m,p,q) * 
    phi_matrix)[0][0])


def M_from_coef(n,a,b,c,m,p,q):
  r"""
  Returns the Mermin operator for a given size *n* and coefficients *a* 
  through *q*.

    M traditionally uses two families of operators, a_n and a'_n, in our 
    case, a_n = a*X+b*Y+c*Z and a'_n = m*X+p*Y+q*Z.
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param real a,b,c,m,p,q: Coefficients for the Mermin operator, used as described 
    above.
  :returns: matrix -- The Mermin operator M_n.
  """
  X = matrix([[0, 1],
        [1, 0]])
  Y = matrix([[0, -i],
        [i, 0]])
  Z = matrix([[1, 0],
        [0, -1]])
  a1 = vector([a,b,c]).normalized()
  a2 = vector([m,p,q]).normalized()
  A1 = a1[0]*X+a1[1]*Y+a1[2]*Z
  A2 = a2[0]*X+a2[1]*Y+a2[2]*Z
  return M(n, A1, A2)


def M(n, a, a_prime):
  r"""
  M_n is defined as such:
    M_n = (1/2)*(M_(n-1).tensor(a + a') + M'_(n-1).tensor(a - a'))
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param matrix a,a_prime: Size 2 hermitian operators, defining M as
    given above.
  :returns: matrix -- A size 2^n operator, following the definition given 
    above.
  """
  if n == 1:
    return a
  return (kronecker(M(n-1, a, a_prime), a + a_prime) + 
    kronecker(M_prime(n-1, a, a_prime), a - a_prime))/2


def M_prime(n, a, a_prime):
  r"""
  M'_n is defined as such:
    M'_n = (1/2)*(M'_(n-1).tensor(a + a') + M_(n-1).tensor(a' - a))
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param matrix a,a_prime: Size 2 hermitian operators, defining M as
    given above.
  :returns: matrix -- A size 2^n operator, following the definition given 
    above.
  """
  if n == 1:
    return a_prime
  return (kronecker(M(n-1, a, a_prime), a_prime - a) + 
    kronecker(M_prime(n-1, a, a_prime), a + a_prime))/2


def mermin_coef_opti_all(phi, verbose=False):
  r"""
  Returns the Mermin operator' coefficients maximizing tr(M_n * rho) for a 
  given input phi (where `rho = |phi><phi|` is the density matrix corresponding
  to the state phi).

  :param vector[complex] phi: State used for the optimization of tr(M_n * rho)
    (only single item searches are supported for now).
  :param bool verbose: If *verbose* then extra run information will be displayed in 
    terminal.
  :returns: list[real], real -- Coefficients of the optimal Mermin operator for 
    *target_state* in th Grover algorithm and the value reached.
  """
  n = log(len(phi))/log(2)
  rho = (matrix(phi).transpose())*(matrix(phi).conjugate())

  def M_func(_a_a_prime_coefs):
    _a_coefs, _a_prime_coefs = coefficients_packing(_a_a_prime_coefs)
    return M_eval_all(n, _a_coefs, _a_prime_coefs, rho)

  def normalizer_func(_a_a_prime_coefs):
    _a_coefs, _a_prime_coefs = coefficients_packing(_a_a_prime_coefs)
    _a_coefs_normalized = []
    for _a_coef in _a_coefs:
      _a_coefs_normalized.append(list(vector(_a_coef).normalized()))
    _a_prime_coefs_normalized = []
    for _a_prime_coef in _a_prime_coefs:
      _a_prime_coefs_normalized.append(list(vector(_a_prime_coef).normalized()))
    return coefficients_unpacking(_a_coefs_normalized,_a_prime_coefs_normalized)

  _a_a_prime_coefs_opti,v = optimize(M_func, [1]*3*n*2, 5, 10**(-2), 10**3)
  # _a_a_prime_coefs_opti,v = optimize_normalized(M_func, normalizer_func, [1]*3*n*2, 5, 10**(-2), 10**4)
  if verbose:
    print("Optimization completed result: optimum reached at (a, a') = " + \
      str(_a_a_prime_coefs_opti) + " (packed coefficients) with a value of " + \
      "v = " + str(v))
  return _a_a_prime_coefs_opti,v


def M_eval_all(_n, _a_coefs, _a_prime_coefs, _rho):
  r"""
  This function evaluates tr(M_n * rho) with *a* and *a'* describing Mn, the 
  Mermin operator.

  The coefficients must be given in the following shape:
  >>> [[a,b,c], [d,e,f], ...]
  and will result in the following family of observables:
  >>> a[0] = a*X + b*Y + c*Z
  >>> a[1] = d*X + e*Y + f*Z
  >>> ...

  :param int _n: Size of the system.
  :param list[list[real]] _a_coefs, _a_prime_coefs: Coefficients for the Mermin 
    operator, used as described above.
  :param matrix[complex] rho: Density matrix of the state to be evaluated with 
    M_n
  :returns: complex
  """
  return abs((M_from_coef_all(_n, _a_coefs, _a_prime_coefs)*_rho).trace())


def M_from_coef_all(_n, _a_coefs, _a_prime_coefs):
  r"""
  Returns the Mermin operator for a given size *n* and coefficients of *a* 
    and *a'*

  The coefficients must be given in the following shape:
  >>> [[a,b,c], [d,e,f], ...]
  and will result in the following family of observables:
  >>> a[0] = a*X + b*Y + c*Z
  >>> a[1] = d*X + e*Y + f*Z
  >>> ...
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param list[list[real]] _a_coefs, _a_prime_coefs: Coefficients for the Mermin 
    operator, used as described above.
  :returns: matrix -- The Mermin operator M_n.
  """
  X = matrix([[0, 1],
        [1, 0]])
  Y = matrix([[0, -i],
        [i, 0]])
  Z = matrix([[1, 0],
        [0, -1]])

  _a = []
  for _a_coef in _a_coefs:
    _a_n = vector(_a_coef).normalized() 
    _a.append(_a_n[0]*X + _a_n[1]*Y + _a_n[2]*Z)
  _a_prime = []
  for _a_prime_coef in _a_prime_coefs:
    _a_prime_n = vector(_a_prime_coef).normalized() 
    _a_prime.append(_a_prime_n[0]*X + _a_prime_n[1]*Y + _a_prime_n[2]*Z)

  return M_all(_n, _a, _a_prime)


def M_all(_n, _a, _a_prime):
  r"""
  M_n is defined as such:
    M_n = (1/2)*(M_(n-1).tensor(a_n + a_n') + M'_(n-1).tensor(a-n - a-n'))
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param list[matrix] a,a_prime: List of size 2 hermitian operators, defining M
    as given above.
  :returns: matrix -- A size 2^n operator, following the definition given 
    above.
  """
  if _n == 1:
    return _a[0]
  return (kronecker(M_all(_n-1, _a, _a_prime), _a[_n-1] + _a_prime[_n-1]) + 
    kronecker(M_prime_all(_n-1, _a, _a_prime), _a[_n-1] - _a_prime[_n-1]))/2


def M_prime_all(_n, _a, _a_prime):
  r"""
  M'_n is defined as such:
    M'_n = (1/2)*(M'_(n-1).tensor(a_n + a_n') + M_(n-1).tensor(a_n' - a_n))
  
  :param int n: Iteration for the Mermin operator (determines its size).
  :param list[matrix] a,a_prime: List of size 2 hermitian operators, defining M' 
    as given above.
  :returns: matrix -- A size 2^n operator, following the definition given 
    above.
  """
  if _n == 1:
    return _a_prime[0]
  return (kronecker(M_all(_n-1, _a, _a_prime), _a_prime[_n-1] - _a[_n-1]) + 
    kronecker(M_prime_all(_n-1, _a, _a_prime), _a[_n-1] + _a_prime[_n-1]))/2


def coefficients_packing(_a_a_prime_coefs):
  r"""
  Packs a list of elements in two lists of lists of three elements

  Example:
  >>>  coefficients_packing([1,2,3,4,5,6,7,8,9,10,11,12])
  ([[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]])

  This function is used to interface above *..._all* functions and the 
    *optimize* function.

  :param list[any] _a_a_prime_coefs: List of elements.
  :returns: tuple[list[list[any]]] -- Lists of lists of elements as described 
    above. 
  """
  _a_a_prime_coefs_packed = [[_a_a_prime_coefs[3*_i], _a_a_prime_coefs[3*_i+1], 
    _a_a_prime_coefs[3*_i+2]] for _i in range(len(_a_a_prime_coefs)/3)]
  _a_coefs = _a_a_prime_coefs_packed[:len(_a_a_prime_coefs_packed)/2]
  _a_prime_coefs = _a_a_prime_coefs_packed[len(_a_a_prime_coefs_packed)/2:]
  return _a_coefs, _a_prime_coefs


def coefficients_unpacking(_a_coefs, _a_prime_coefs):
  r"""
  Unpacks two lists of lists of three elements to one list of elements

  Example:
  >>>  coefficients_unpacking([[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]])
  [1,2,3,4,5,6,7,8,9,10,11,12]

  This function is used to interface above *..._all* functions and the 
    *optimize* function.

  :param tuple[list[list[any]]] _a_coefs, _a_prime_coefs: Lists of lists of 
    elements as described above. 
  :returns: list[any] -- List of elements.
  """
  _a_a_prime_coefs = [k for element in _a_coefs+_a_prime_coefs for k in element]
  return _a_a_prime_coefs
  