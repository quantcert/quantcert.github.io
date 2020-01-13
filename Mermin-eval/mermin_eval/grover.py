#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is aimed to provide the user means to study the Grover algorithm and
especially how the entanglement changes during its execution using the Mermin 
evaluation.
"""
from sage.all import *

import csv
import time
from copy import copy
import os

from mermin_eval import *
from run_circuit import *


def target_state_ket_list_to_vector(target_state_ket):
  r""" 
  Converts the target state from a digit list to a vector.

  Example:
    `|5>_4 = |0101> = (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)` so :
    
    >>> target_state_ket_list_to_vector([0,1,0,1])
    (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

  :param list[int] target_state_ket: List of digit of the boolean expression of 
    the sate looked for.
  :returns: vector -- Vector corresponding to the target state.
  """
  n = len(target_state_ket)
  t = 0
  target_state_ket_reversed = copy(target_state_ket)
  target_state_ket_reversed.reverse()
  for i in range(n):
    t += target_state_ket_reversed[i]*2**i
  result = zero_vector(field, 2**n)
  result[t] = 1
  return result


def target_state_ket_vector_to_string(target_state_ket):
  r""" 
  Converts the target state from a vector to a string containing a digit list.

  Example:
    `|5>_4 = |0101> = (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)` so :
    
    >>> v = vector((0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    >>> target_state_ket_vector_to_string(v)
    '0101'

  :param vector target_state_ket: Vector corresponding to the target state.
  :returns: string -- List of digit of the boolean expression of the searched 
    sate.
  """
  N = len(target_state_ket)
  n = log(N)/log(2)
  i = 0
  while i < N and target_state_ket[i] != 1:
    i += 1
  if target_state_ket[i] != 1:
    raise ValueError("argument target_state_ket is not in the appropriate" +
      " shape: it should be a vector of size 2^n with only zeros except" +
      " for one coefficient being 1, given : " + str(target_state_ket))
  result = "{0:b}".format(i)
  while len(result) < n:
    result = "0" + result
  return result


def target_state_ket_string_to_vector(target_state_ket):
  r""" 
  Converts the target state from a string containing a digit list to a vector.

  Example:
    ``|5>_4 = |0101> = (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)`` so :
    
    >>> target_state_ket_string_to_vector("0101")
    (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

  
  :param str target_state_ket: List of digit of the boolean expression of the 
    searched sate.
  :returns: vector -- Vector corresponding to the target state.
  """
  n = len(target_state_ket)
  t = int(target_state_ket, 2)
  result = zero_vector(field, 2**n)
  result[t] = 1
  return result


def grover(target_state_vector, verbose=False, file_name=None, use_precomputed=False, artificial=False):
  r""" 
  Prints in terminal or in file ``file_name`` the Mermin evaluation of each step 
  of the Grover algorithm. 

  Example:
    >>> grover(target_state_ket_string_to_vector("0000"))
    [0.168628057515893, 1.32148634347176, 1.01189404012546, 0.469906068135870]

  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :param bool verbose: If ``verbose`` then extra run information will be 
    displayed in terminal.
  :param str file_name: File name for the registration of the Mermin evaluation 
    for each step of the algorithm, in csv format.
  :param bool use_precomputed: for some states, the optimal Mermin polynomial 
    has been precomputed, use this option to speed up the overall computation. 
  :param bool artificial: Due to technological limits, it is not always possible 
    to compute the states of Grover's algorithm in a realistic simulator. This 
    parameters allows the user to use a non realistic simulator which is way 
    quicker.
  :returns: any -- Result of this function depends on file_name. If a file name
    is given ``qft_main`` returns None, otherwise, it returns an array with the 
    evaluation values.
  """
  end_loop_states = time_step("Run", grover_run, 
    (target_state_vector, artificial, verbose), verbose)

  M_opt = time_step("Optimization", grover_optimize, 
    (target_state_vector, use_precomputed, verbose), verbose)

  return time_step("Evaluation", grover_evaluate, 
    (end_loop_states, M_opt, file_name), verbose)


def time_step(step_name, step_function, step_args, verbose=False):
  r""" 
  Times a function, with time information displayed if ``verbose``.
  
  Example:
    >>> result = time_step("Add_10", lambda a : a + 10, [5], True)
    Add_10 starts now
    Add_10 complete, took 1.81198120117e-05s
    >>> result
    15

  :param str setp_name: Name of the step, used for display.
  :param function step_function: Function to be timed.
  :param tuple step_args: Arguments of the function.
  :param bool verbose: If `verbose` then extra run information will be displayed 
    in terminal.
  :returns: any -- The output of ``step_function(*step_args)```.
  """
  if verbose:
    print(step_name + " starts now")
  start_time = time.time()

  result = step_function(*step_args)

  end_time = time.time()
  interval = end_time - start_time
  if verbose:
    print(step_name + " complete, took " + str(interval) + "s")

  return result


def grover_run(target_state_vector, artificial=False, verbose=False):
  r""" 
  Runs a simulation of the grover algorithm.
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :param bool artificial: due to technological limits, it is not always possible 
    to compute the states of Grover's algorithm in a realistic simulator. This 
    parameters allows the user to use a non realistic simulator which is way 
    quicker.
  :param bool verbose: If `verbose` then extra run information will be displayed 
    in terminal.
  :returns: list[vectors] -- the states at the end of each loop.
  """
  if artificial:
    end_loop_states = grover_artifical(target_state_vector)
  else:
    end_loop_states = grover_vanilla(target_state_vector, verbose)

  return end_loop_states


def grover_optimize(target_state_vector, use_precomputed=False, verbose=False):
  r""" 
  Computes the optimal Mermin operator to evaluate the entanglement in the
  Grover algorithm
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :param bool use_precomputed: For some states, the optimal Mermin polynomial 
    has been precomputed, use this option to speed up the overall computation.
  :param bool verbose: If `verbose` then extra run information will be displayed 
    in terminal.
  :returns: matrix -- The optimal Mermin operator.
  """
  if use_precomputed:
    precomputed_filename = "precomputed_opti.csv"
  else: 
    precomputed_filename = None
  M_opt = mermin_operator_opti(target_state_vector, 
    precomputed_filename=precomputed_filename, verbose=verbose)

  return M_opt


def grover_evaluate(end_loop_states, M_opt, file_name):
  r""" 
  Uses the previously found optimal mermin operator to evaluate the entanglement 
  for each state in ``end_loop_states``.
  
  :param list[matrix] end_loop_states: The states at the end of each loop, these
    are vectors transformed to column matrix.
  :param matrix M_opt: The optimal Mermin operator.
  :param str file_name: File name for the registration of the Mermin evaluation 
    for each step of the algorithm, in csv format.
  :returns: any -- Result of this function depends on file_name. If a file name
    is given ``qft_main`` returns None, otherwise, it returns an array with the 
    evaluation values.
  """
  if file_name == None:
    result = []
    for state in end_loop_states:
      result.append(abs((state.transpose().conjugate()*M_opt*state)[0][0]))
    return result
  else :
    lines, i = [["iteration", "intricationValue"]], 0
    for state in end_loop_states:
      lines.append([i, (state.transpose().conjugate()*M_opt*state)[0][0].real()])
      i += 1
    if not os.path.exists(os.path.dirname(file_name)) and os.path.dirname(file_name) != '':
      os.makedirs(os.path.dirname(file_name))
    with open(file_name, 'wb+') as writeFile:
      writer = csv.writer(writeFile)
      writer.writerows(lines)
    writeFile.close()
    return None


def grover_vanilla(target_state_vector, verbose=False):
  r"""
  Computes the successive states at the end of the loop using a realistic
  simulation of the execution of Grover's algorithm (using the simulator
  from ``run_circuit.sage``).

  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :param bool verbose: If `verbose` then extra run information will be displayed 
    in terminal.
  :returns: list[vector] -- List of states after each application of the 
    diffusion operator (as well as the initial state). 
  """
  N = len(target_state_vector)
  if verbose:
    print("Layers preparation")
  layers, k_opt = grover_layers_kopt(target_state_vector)

  V0 = vector([0, 1] + [0]*(2*N-2))
  if verbose:
    print("Layers preparation done")

  states = run(layers, V0)[0]
  end_loop_states = [matrix(states[1]).transpose()]
  for i in range(k_opt):
    end_loop_states.append(matrix(states[4*(i)+5]).transpose())
  for i in range(len(end_loop_states)):
    halved = vector([end_loop_states[i][2*k,0] for k in 
      range(end_loop_states[i].nrows()/2)])
    end_loop_states[i] = matrix(halved.normalized()).simplify().transpose()
  return end_loop_states


def grover_layers_kopt(target_state_vector):
  r"""
  Computes the circuit for Grover's algorithm using the circuit format used for 
  the ``run`` function from ``run_circuit.sage``.
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :returns: (list[list[matrix]],int) -- Circuit for Grover's algorithm and 
    optimal value found for ``k_opt``.
  """
  X = matrix([[0, 1],
        [1, 0]])
  H = matrix(field, [[1, 1], 
             [1, -1]])/sqrt(2)
  I2 = matrix.identity(2)

  N = len(target_state_vector)
  n = log(N)/log(2)

  ket_zero = matrix([1]+[0]*(N-1))
  symmetry = 2*ket_zero.transpose()*ket_zero-matrix.identity(N)

  loop = [[oracle(target_state_vector)], [H]*n+[I2],[symmetry,I2],[H]*n+[I2]]

  k_opt = round((pi/4)*sqrt(N))

  layers = [[H]*(n+1)]
  for _ in range(k_opt):
    layers += loop

  return layers, k_opt


def oracle(target_state_vector):
  r""" 
  The oracle O satisfies `O|x,y> = |x,f(x)+y>` where `f(x)=1` if `x` is the 
  element looked for and `f(x)=0` otherwise.

  Example:
    The oracle flips the qubit where the target is, so :

    >>> oracle(target_state_ket_string_to_vector("101"))
    [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
    
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :returns: matrix -- The matrix corresponding to the oracle.
  """
  N = len(target_state_vector)
  state_matrix = matrix(target_state_vector)
  I2n = matrix.identity(N*2)
  I2 = matrix.identity(2)
  X = matrix([[0, 1],
        [1, 0]])
  return I2n + kronecker(state_matrix.transpose()*state_matrix, X - I2)


def grover_artifical(target_state_vector):
  r"""
  To accelerate the computation of the states, the alternative to 
  ``grover_vanilla`` is to work directly on the vectors. INdeed, we know what 
  the effects of each steps are, so we don't need to apply the gates we can 
  simply change the vectors accordingly. For big dimensions, this is a way more 
  efficient method.
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :returns: list[vector] -- List of states after each application of the 
    diffusion operator (as well as the initial state). 
  """
  N = len(target_state_vector)
  n = log(N)/log(2)
  k_opt = round((pi/4)*sqrt(N))

  V0 = vector([1]+[0]*(N-1))

  H = matrix(field, [[1, 1], 
             [1, -1]])/sqrt(2)
  hadamard_layer = kronecker_power(H, n)
  V = hadamard_layer * V0
  end_loop_states = [matrix(V).transpose()]

  for k in range(k_opt):
    V = oracle_artificial(target_state_vector, V)
    V = diffusion_artificial(V)
    end_loop_states.append(matrix(V).transpose())

  return end_loop_states


def oracle_artificial(target_state_vector, V):
  r"""
  Flips the coefficient of ``V`` corresponding to the state 
  ``target_state_vector``.
  
  :param vector[int] target_state_vector: State searched by Grover's algorithm 
    (only single item searches are supported for now).
  :param vector V: The running state in Grover's algorithm.
  :returns: vector -- ``V`` with the correct coefficient flipped.
  """
  result = copy(V)
  for i in range(len(target_state_vector)):
    if target_state_vector[i] == 1:
      result[i] = -V[i]
  return result


def diffusion_artificial(V):
  r"""
  Performs the inversion about the mean for ``V``.
  
  :param vector V: The running state in Grover's algorithm .
  :returns: vector -- ``V`` Inverted about the mean.
  """
  result = vector(field, [0]*len(V))
  m = mean(V)
  for i in range(len(V)):
    result[i] = 2 * m - V[i]
  return result
