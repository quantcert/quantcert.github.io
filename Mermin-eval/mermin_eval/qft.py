#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is aimed to provide the user means to study the QFT and especially 
how entanglement changes during its execution using the Mermin evaluation.
"""
from sage.all import *

import csv
import time

from mermin_eval import *
from run_circuit import *


def qft_main(state, verbose=False, file_name=None):
  r""" 
  Prints in terminal or in file ``file_name`` the Mermin evaluation of each 
  step of the QFT. 

  Example:
    >>> grover(periodic_state(0,11,4))0,1.99823485241887
    1.99933058720672
    1.99761683801972
    1.84600827724524
    1.66149864543535
    1.66010335826574
    1.66127978203391
    2.17163453049907
    2.17187381801221
    1.54230165695219
    1.54129902140376
    1.54169705834015

  :param vector[int] state: State on which the operator wants the QFT performed 
    of, this is usually a periodic state.
  :param bool verbose: If *verbose* then extra run information will be displayed 
    in terminal.
  :param str file_name: File name for the registration of the Mermin evaluation 
    for each step of the algorithm, in csv format.
  :returns: any -- Result of this function depends on file_name. If a file name
    is given `qft_main` returns None, otherwise, it returns an array with the 
    evaluation values .
  """
  states = qft_run(state, verbose)

  return qft_evaluate(states, verbose, file_name)


def qft_run(state, verbose=False):
  r""" 
  Runs a simulation of the QFT.

  Example:
    >>> qft_run(periodic_state(2,1,3))
    [(0, 0, 1/6*sqrt(6), 1/6*sqrt(6), 1/6*sqrt(6), 1/6*sqrt(6), 1/6*sqrt(6), 1/6*sqrt(6)),
     (1/12*sqrt(6)*sqrt(2), 1/12*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), -1/12*sqrt(6)*sqrt(2), -1/12*sqrt(6)*sqrt(2), 0, 0),
     (1/12*sqrt(6)*sqrt(2), 1/12*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), -1/12*sqrt(6)*sqrt(2), -1/12*sqrt(6)*sqrt(2), 0, 0),
     (1/12*sqrt(6)*sqrt(2), 1/12*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), 1/6*sqrt(6)*sqrt(2), -1/12*sqrt(6)*sqrt(2), -(1/12*I + 1/12)*sqrt(6), 0, 0),
     (1/4*sqrt(6), 1/4*sqrt(6), -1/12*sqrt(6), -1/12*sqrt(6), -1/12*sqrt(6), -(1/24*I + 1/24)*sqrt(6)*sqrt(2), -1/12*sqrt(6), -(1/24*I + 1/24)*sqrt(6)*sqrt(2)),
     (1/4*sqrt(6), 1/4*sqrt(6), -1/12*sqrt(6), -1/12*I*sqrt(6), -1/12*sqrt(6), -(1/24*I + 1/24)*sqrt(6)*sqrt(2), -1/12*sqrt(6), -(1/24*I - 1/24)*sqrt(6)*sqrt(2)),
     (1/4*sqrt(6)*sqrt(2), 0, -(1/24*I + 1/24)*sqrt(6)*sqrt(2), (1/24*I - 1/24)*sqrt(6)*sqrt(2), -1/24*sqrt(6)*sqrt(2) - (1/24*I + 1/24)*sqrt(6), -1/24*sqrt(6)*sqrt(2) + (1/24*I + 1/24)*sqrt(6), -1/24*sqrt(6)*sqrt(2) - (1/24*I - 1/24)*sqrt(6), -1/24*sqrt(6)*sqrt(2) + (1/24*I - 1/24)*sqrt(6)),
     (1/4*sqrt(6)*sqrt(2), -1/24*sqrt(6)*sqrt(2) - (1/24*I + 1/24)*sqrt(6), -(1/24*I + 1/24)*sqrt(6)*sqrt(2), -1/24*sqrt(6)*sqrt(2) - (1/24*I - 1/24)*sqrt(6), 0, -1/24*sqrt(6)*sqrt(2) + (1/24*I + 1/24)*sqrt(6), (1/24*I - 1/24)*sqrt(6)*sqrt(2), -1/24*sqrt(6)*sqrt(2) + (1/24*I - 1/24)*sqrt(6))]
  
  :param vector[int] state: State on which the operator wants the QFT performed 
    of, this is usually a periodic state.
  :param bool verbose: If ``verbose`` then extra run information will be 
    displayed in terminal.
  :returns: list[vectors] -- the list of states after each step.
  """
  layers = qft_layers(state)
  if verbose:
    print "Layers completed"

  states, _ = run(layers, state)
  if verbose:
    print "Run completed"

  return states


def qft_layers(state):
  r"""
  Computes the circuit for the QFT using the circuit format  used for the 
  ``run`` function from ``run_circuit.sage``.

  Example:
    >>> qft_layers(periodic_state(2,1,3))
    [
      [
        [ 1/2*sqrt(2)  1/2*sqrt(2)]  [1 0]  [1 0]
        [ 1/2*sqrt(2) -1/2*sqrt(2)], [0 1], [0 1]
      ],[
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 I 0]
        [0 0 0 0 0 0 0 I]
      ],[
        [1 0 0 0 0                     0 0                     0]
        [0 1 0 0 0                     0 0                     0]
        [0 0 1 0 0                     0 0                     0]
        [0 0 0 1 0                     0 0                     0]
        [0 0 0 0 1                     0 0                     0]
        [0 0 0 0 0 (1/2*I + 1/2)*sqrt(2) 0                     0]
        [0 0 0 0 0                     0 1                     0]
        [0 0 0 0 0                     0 0 (1/2*I + 1/2)*sqrt(2)]
      ],[
        [1 0]  [ 1/2*sqrt(2)  1/2*sqrt(2)]  [1 0]
        [0 1], [ 1/2*sqrt(2) -1/2*sqrt(2)], [0 1]
      ],[
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 I 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 I]
      ],[
        [1 0]  [1 0]  [ 1/2*sqrt(2)  1/2*sqrt(2)]
        [0 1], [0 1], [ 1/2*sqrt(2) -1/2*sqrt(2)]
      ],[
        [1 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 1 0]
        [0 1 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 1]
      ]
    ]

  
  :param vector[int] state: State on which the operator wants the QFT performed 
    of, this is usually a periodic state.
  :returns: (list[list[matrix]],int) -- Circuit for the QFT.
  """
  H = matrix(field, [[1,  1], 
                     [1, -1]])/sqrt(2)

  I2 = matrix.identity(field, 2)

  def swap(wire1=0,wire2=1,size=2):
    result = matrix(field, 2**size)
    # matrix uses little endian and bit values modification uses big endian so 
    # we need to convert between the two
    w1_big = size-1-wire1
    w2_big = size-1-wire2
    for i in range(2**size):
      # we exchange the bit corresponding to each wire
      result[
        set_bit_value(
          set_bit_value(i,w1_big,bit_value(i,w2_big)),
          w2_big,bit_value(i,w1_big)
        ),i] = 1
    return result

  def R(k,target,control,size):
    cR = matrix.identity(field, 4)
    cR[3,3] = exp(2*i*pi/(2**k))
    result = kronecker(cR, matrix.identity(field, 2**(size-2)))
    result = \
      swap(target, 0, size) * swap(control, 1, size) * \
      result *\
      swap(1, control, size) * swap(0, target, size)
    return result

  nWires = log(len(state))/log(2)

  layers = []

  for wire in range(nWires):
    layers.append([I2]*wire + [H] + [I2]*(nWires-wire-1))
    for k in range(2, nWires-(wire-1)):
      layers.append([R(k, wire, k+(wire-1), nWires)])

  global_swap = matrix.identity(field, 2**nWires)
  for wire in range(nWires/2):
    global_swap *= swap(wire, nWires-1-wire, nWires)
  layers.append([global_swap])

  return layers


def qft_evaluate(states, verbose=False, file_name=None):
  r""" 
  Computes the Mermin Evaluation for each state in ``states``.

  Example:
    >>> qft_evaluate([periodic_state(1,5,4)])
    [1.33201398251237]
  
  :param list[vector] states: The states after each step of the QFT.
  :param bool verbose: If `verbose` then extra run information will be displayed 
    in terminal.
  :param str file_name: File name for the registration of the Mermin evaluation 
    for each step of the algorithm, in csv format.
  :returns: any -- Result of this function depends on file_name. If a file name
    is given `qft_main` returns None, otherwise, it returns an array with the 
    evaluation values.
  """
  result = []
  for state in states:
    # we are setting the calculations in floating numbers (instead on symbolic
    # to speed them up)
    result.append(mermin_coef_opti_all(vector(CC,state), verbose)[1])
  if file_name == None:
    return result
  else :
    lines, i = [["iteration", "intricationValue"]], 0
    for value in result:
      lines.append([i, value])
      i += 1
    if not os.path.exists(os.path.dirname(file_name)) and os.path.dirname(file_name) != '':
      os.makedirs(os.path.dirname(file_name))
    with open(file_name, 'w') as writeFile:
      writer = csv.writer(writeFile)
      writer.writerows(lines)
    writeFile.close()
    return None


def qft_matrix(size):
  r""" 
  This function should compute a matrix equivalent to the whole QFT operation. 
  It has not been tester much though and should not be relied on for now.

  Example:
    >>> qft_matrix(3)
    [   1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)    1/4*sqrt(2)]
    [   1/4*sqrt(2)    1/4*I + 1/4  1/4*I*sqrt(2)    1/4*I - 1/4   -1/4*sqrt(2)   -1/4*I - 1/4 -1/4*I*sqrt(2)   -1/4*I + 1/4]
    [   1/4*sqrt(2)  1/4*I*sqrt(2)   -1/4*sqrt(2) -1/4*I*sqrt(2)    1/4*sqrt(2)  1/4*I*sqrt(2)   -1/4*sqrt(2) -1/4*I*sqrt(2)]
    [   1/4*sqrt(2)    1/4*I - 1/4 -1/4*I*sqrt(2)    1/4*I + 1/4   -1/4*sqrt(2)   -1/4*I + 1/4  1/4*I*sqrt(2)   -1/4*I - 1/4]
    [   1/4*sqrt(2)   -1/4*sqrt(2)    1/4*sqrt(2)   -1/4*sqrt(2)    1/4*sqrt(2)   -1/4*sqrt(2)    1/4*sqrt(2)   -1/4*sqrt(2)]
    [   1/4*sqrt(2)   -1/4*I - 1/4  1/4*I*sqrt(2)   -1/4*I + 1/4   -1/4*sqrt(2)    1/4*I + 1/4 -1/4*I*sqrt(2)    1/4*I - 1/4]
    [   1/4*sqrt(2) -1/4*I*sqrt(2)   -1/4*sqrt(2)  1/4*I*sqrt(2)    1/4*sqrt(2) -1/4*I*sqrt(2)   -1/4*sqrt(2)  1/4*I*sqrt(2)]
    [   1/4*sqrt(2)   -1/4*I + 1/4 -1/4*I*sqrt(2)   -1/4*I - 1/4   -1/4*sqrt(2)    1/4*I - 1/4  1/4*I*sqrt(2)    1/4*I + 1/4]

  :param int size: The size of the QFT (number of qubits).
  :returns: matrix -- The matrix corresponding to the QFT.
  """
  N = 2**size
  w = exp(2*i*pi/N)
  result = matrix(SR, N)
  for k in range(N):
    for l in range(N):
      result[k,l] = w**(k*l)
  return result/sqrt(N)


def periodic_state(l,r,nWires):
  r""" 
  Returns the periodic state `|\varphi^{l,r}>` of size `2^{nWires}`. We have:
    
    `|\varphi^{l,r}> = \sum_{i=0}^{A-1}|l+ir>/sqrt(A)` with
    `A = floor((2^{nWires}-l)/r)+1`

    In this definition, ``l`` is the shift of the state, and ``r`` is the period 
    of the state.

  Example:
    Since
    `|\varphi^{1,5}> = (|1>+|6>+|11>)/sqrt(3)=(|0001>+|0110>+|1011>)/sqrt(3)`,

    >>> periodic_state(1,5,4)
    (0, 1/3*sqrt(3), 0, 0, 0, 0, 1/3*sqrt(3), 0, 0, 0, 0, 1/3*sqrt(3), 0, 0, 0, 0)

  :param int l: The shift of the state.
  :param int r: The period of the state.
  :param int nWires: The size of the system (number of qubits).
  :returns: vector -- The state defined by ``l``, ``r`` and ``nWires`` according
    to the definition given above.
  """
  N = 2**nWires
  result = vector(SR, N)
  for k in range(ceil((N-l)/r)):
    result[l+k*r] = 1
  return result.normalized()


def bit_value(n, k, base=2):
  """ 
  Returns the value of the ``k`` `^{th}` digit of the integer ``n`` given in 
  base ``base``.

  Example:
    >>> bit_value(5, 0)
    1
    >>> bit_value(5, 1)
    0
    >>> bit_value(15, 0, base=10)
    5

  :param int n: Integer studied.
  :param int k: Digit desired.
  :param int base: Base in which ``n`` is studied.
  :returns: int -- Value of the ``k`` `^{th}` digit of the integer ``n`` given 
    in base ``base``.
  """
  return floor(abs(n)/base**(k)) - floor(abs(n)/base**(k+1))*base


def set_bit_value(n, k, value, base=2):
  """ 
  Returns ``n`` with its ``k`` `^{th}` bit set to ``value`` in base ``base``.

  Example:
    >>> set_bit_value(5, 0, 0)
    4
    >>> set_bit_value(5, 1, 0)
    5
    >>> set_bit_value(15, 0, 9, base=10)
    19

  :param int n: Integer modified.
  :param int k: Digit to change.
  :param int value: Value wanted for the ``k`` `^{th}` digit of ``n`` (must be 
    between 0 and ``base-1``).
  :param int base: Base in which ``n`` is modified.
  :returns: int -- ``n`` with its ``k`` `^{th}` bit set to ``value`` in base 
    ``base``.
  """
  return n - bit_value(n, k, base)*(2**k) + value*(2**k)
