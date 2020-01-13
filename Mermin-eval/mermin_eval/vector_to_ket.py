#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is an ease of life adding to SageMath: it complements the ``latex``
method to quantum ket notations. 
"""
from sage.all import *

def vector_to_ket(v):
  r"""Transform a sage vector to a Latex ket notation
  
  Example:
    >>> v = vector(SR,[1,0,0,1])
    >>> vector_to_ket(v)
    1 \ket{00} + 1 \ket{11}

  :param vector v: the vector corresponding to a pure state (size must be a power of 
    two) 
  :returns: str -- a string corresponding to the Latex code for the ket 
    notation of the input vector
  """
  nb_qubits = log(len(v))/log(2)
  if not nb_qubits.is_integer():
    raise ValueError("v should be a vector with a dimension equal to a " + 
      "power of two but the given vector has a size of " + str(len(v)))
  latex_str = ""
  first = True
  index = 0;
  for i in range(len(v)):
    if v[i] != 0:
      if first:
        latex_str += latex(v[i]) + int_index_to_ket(i, nb_qubits)
        first = False
      else:
        latex_str += "+" + latex(v[i]) + int_index_to_ket(i, nb_qubits)
  return latex_str

def int_index_to_ket(index, register_size):
  r"""Makes a int from a register into a binary ket notation.
  IMPORTANT: this notation requires the *physics* package in Latex

  Example:
    >>> int_index_to_ket(2,2)
    '\\ket{10}'

  :param int index: the index to be transformed 
  :param int register_size: the number of qubits in the system 
  :returns: string -- a string in latex format
  """
  ket = ""
  for k in range(register_size):
    ket = str((index & ( 1 << k )) >> k) + ket
  return  "\\ket{" + ket + "}"