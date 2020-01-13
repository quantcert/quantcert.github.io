#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module provides a simple quantum circuit simulator.
"""
from sage.all import *

import numpy as np
import copy as cp

from vector_to_ket import *

field = SR
# field = UniversalCyclotomicField()
# e8 = field.gen(8)
# sqrt2 = e8 + conjugate(e8) # fix for UniversalCyclotomicField


def kronecker(a, b):
  r"""
  Computes the Kronecker product of ``a`` and ``b``.
  
  Example:
    >>> a = matrix([[1,0],[0,1]])
    >>> b = matrix([[1,2],[3,4]])
    >>> kronecker(a,b)
    [1 2 0 0]
    [3 4 0 0]
    [0 0 1 2]
    [0 0 3 4]

  :param matrix,vector a,b: Operands for the kronecker operator. 
  :returns: matrix *or* vector --  The kronecker product of ``a`` and ``b`` 
    (return type is the same a type of ``a``).
  """
  a_matrix = matrix(a)
  b_matrix = matrix(b)
  ma, na = a_matrix.nrows(), a_matrix.ncols()
  mb, nb = b_matrix.nrows(), b_matrix.ncols()
  mr, nr = ma * mb, na * nb
  result = Matrix(field, mr, nr)
  for i in range(0, mr):
    for j in range(0, nr):
      result[i, j] = a_matrix[i // mb, j // nb] * b_matrix[i % mb, j % nb]
  if "Vector" in str(parent(a)):
    return vector(result)
  elif "Matrix" in str(parent(a)):
    return matrix(result)
  else:
    raise ValueError("Input type of first parameter must be Vector or " + 
      "Matrix but " + str(parent(a)) + " was given")


def kronecker_power(a, n):
  r"""
  Computes the ``n`` `^{th}` Kronecker power of ``a``.

  Example:
    >>> a = matrix([[1,0],[0,2]])
    >>> kronecker_power(a,2)
    [1 0 0 0]
    [0 2 0 0]
    [0 0 2 0]
    [0 0 0 4]
  
  :param matrix,vector a: The matrix (or vector) to be elevated to the 
    ``n`` `^{th}` power.
  :param int n: The power ``a`` has to be elevated to. 
  :returns: matrix *or* vector -- The ``n`` `^{th}` kronecker power of ``a`` 
    (return type is the same a type of ``a``).
  """
  result = a
  for _ in range(n-1):
    result = kronecker(result, a)
  if "Vector" in str(parent(a)):
    return vector(result)
  elif "Matrix" in str(parent(a)):
    return matrix(result)
  else:
    raise ValueError("Input type must be Vector or Matrix but " + 
      str(parent(a)) + " was given")


def output_commant(command_name, command, output=False, output_to_file=False, w_file=None):
  r""" 
  This function is used to print useful informations in various ways, see 
  arguments details for more information.

  Example:
    >>> output_commant("test",vector(SR,[1,0,0,1]), output=True)
    \newcommand{\test}{1 \ket{00} + 1 \ket{11}}

  :param str command_name: The command name. 
  :param any command: The command content (if type is ``str``, will be printed 
    as such; if ``vector``, ``vector_to_ket`` will be called on it and if 
    ``matrix``, ``latex`` method from SageMath will be called on it).
  :param bool output: disables or enables the output.
  :param bool output_to_file: Whether the output should be in the standard 
    output or in a file.
  :param file w_file: The opened file to write the output to.
  :returns: None
  """
  if output:
    if "Vector" in str(parent(command)):
      latex_command = vector_to_ket(command)
    elif "str" in str(parent(a)):
      latex_command = command
    else:
      latex_command = latex(command)
    newcommand = "\\newcommand{\\" + command_name + "}{" + str(latex_command) + "}\n"
    if output_to_file:
      w_file.write(newcommand)
    else:
      print(newcommand)


def run(matrix_layers, V_init, output=False, output_file=False, file=None, vector_name="V", matrix_name="M"):
  r""" 
  Runs the algorithm specified by the matrix_layers.

  Example:  
    >>> I2 = matrix.identity(2)
    >>> I4 = matrix.identity(4)
    >>> H = matrix([[1,1],[1,-1]])/sqrt(2)
    >>> swap = matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    >>> v = vector([1, 0, 0, 0, 0, 0, 0, 0])
    >>> layers = [[I4,H],[H,H,I2],[I2,swap]]
    >>> run(layers,v)
    ([(1, 0, 0, 0, 0, 0, 0, 0),
      (1/2*sqrt(2), 1/2*sqrt(2), 0, 0, 0, 0, 0, 0),
      (1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2)),
      (1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2), 1/4*sqrt(2))],
     [
    [ 1/2*sqrt(2)  1/2*sqrt(2)            0            0            0            0            0            0]
    [ 1/2*sqrt(2) -1/2*sqrt(2)            0            0            0            0            0            0]
    [           0            0  1/2*sqrt(2)  1/2*sqrt(2)            0            0            0            0]
    [           0            0  1/2*sqrt(2) -1/2*sqrt(2)            0            0            0            0]
    [           0            0            0            0  1/2*sqrt(2)  1/2*sqrt(2)            0            0]
    [           0            0            0            0  1/2*sqrt(2) -1/2*sqrt(2)            0            0]
    [           0            0            0            0            0            0  1/2*sqrt(2)  1/2*sqrt(2)]
    [           0            0            0            0            0            0  1/2*sqrt(2) -1/2*sqrt(2)],
    [ 1/2    0  1/2    0  1/2    0  1/2    0]  [1 0 0 0 0 0 0 0]
    [   0  1/2    0  1/2    0  1/2    0  1/2]  [0 0 1 0 0 0 0 0]
    [ 1/2    0 -1/2    0  1/2    0 -1/2    0]  [0 1 0 0 0 0 0 0]
    [   0  1/2    0 -1/2    0  1/2    0 -1/2]  [0 0 0 1 0 0 0 0]
    [ 1/2    0  1/2    0 -1/2    0 -1/2    0]  [0 0 0 0 1 0 0 0]
    [   0  1/2    0  1/2    0 -1/2    0 -1/2]  [0 0 0 0 0 0 1 0]
    [ 1/2    0 -1/2    0 -1/2    0  1/2    0]  [0 0 0 0 0 1 0 0]
    [   0  1/2    0 -1/2    0 -1/2    0  1/2], [0 0 0 0 0 0 0 1]
    ])

  
  :param array[array[matrix]] matrix_layers: Algorithm described as it would 
    be in a circuit.
  :param vector V_init: Initial input of the algorithm.
  :param bool output: If true, outputs are enabled.
  :param bool output_to_file: If true, trace is returned in file, else it is 
    printed.
  :param file file: File to output the commands to.
  :param str vector_name: Base name used for the vectors commands. 
  :param str matrix_name: Base name used for the matrices commands. 
  :returns: vector -- The states along the execution of the algorithm as well as
    the matrix corresponding to each layer.
  """
  V_running = cp.deepcopy(V_init)
  vectors_list = [V_running]
  output_commant(vector_name + int_name(0), V_running, output, output_file, file)
  matrices_list = []
  for i in range(len(matrix_layers)):
    layerMatrix = Matrix(field, [[1]])
    for running_matrix in matrix_layers[i]:
      layerMatrix = kronecker(layerMatrix, running_matrix)
    V_running = layerMatrix*V_running
    matrices_list.append(layerMatrix)
    output_commant(matrix_name + int_name(i+1), layerMatrix, output, output_file, file)
    vectors_list.append(V_running)
    output_commant(vector_name + int_name(i+1), V_running, output, output_file, file)
  return vectors_list, matrices_list


def int_name(num):
  r""" 
  Converts a number to a string composed of the list of its digits in English.
  
  Example:
    >>> int_name(152)
    'onefivetwo'
    
  :param int num: Number to be converted to a list of digits.
  :returns: str -- List of digits of ``num`` in base 10 concatenated.
  """
  digit_str = ["zero", "one", "two", "three", "four", 
        "five", "six", "seven", "eight", "nine"]

  if num == 0:
    return digit_str[0]
  k = 0
  result = ""
  while num+1 > 10**k:
    k += 1
  k -= 1
  while k >= 0:
    result += digit_str[digit(num, k)]
    k -= 1
  return result


def digit(n, k, base = 10):
  r"""
  Computes the digit ``k`` for the integer ``n`` in its representation in base
  ``base``.

  Example:
    >>> digit(152, 1)
    5
    >>> digit(152, 0)
    2

  :param int n: The integer for which the digit is needed.
  :param int k: The number of the digit.
  :param int base: The base used for the representation of ``n``. 
  :returns: int -- The ``k`` `^{th}` digit of ``n`` in base ``base``.
  """
  return floor(abs(n)/base**(k)) - floor(abs(n)/base**(k+1))*base


def print_circuit(name_layers, to_latex=False):
  r"""
  Each name is a tuple with the name of the gate and its dimension 
  
  Example:
    >>> circuit = [[('I',2),('H',1)],[('H',1),('H',1),('I',1)],[('I',1),('S',2)]]
    >>> print_circuit(circuit, to_latex=False)
    ---H---
    ---H-S-
    -H---|-
    >>> print_circuit(circuit, to_latex=True)
    \begin{align*}
     \Qcircuit @C=1em @R=.7em {
      & \qw       & \gate{H}  & \qw               & \qw\\
      & \qw       & \gate{H}  & \multigate{1}{S}  & \qw\\
      & \gate{H}  & \qw       & \ghost{S}         & \qw
     }
    \end{align*}
  
  ``@multi-gost_``, ``@multi-source_`` and ``@multi-size_`` are reserved names, 
  they should not be used as a gate name.

  :param list[List[(str,int)]] name_layers: Circuit name formalism in the shape 
    of the example given above. 
    First layer should not be empty.
    sum([item[1] for item in layer]) should be constant for layer in layers
  :param bool latex: If true, a string is returned containing the circuit in
    the format given by the LaTeX package ``qcircuit``. Otherwise, the string
    returned is under a custom format meant to be easily readable in the 
    terminal.
  :returns: str -- The circuit in the format described above.
  """
  def gate_name_to_string(name, _to_latex=False):
    r""" 
    Converts the gate name to the printed result, the reserved names are used 
    as follow :

    For a multi gate, the first wire will have as gate name 
    ``@multi-source_*name*@multi-size_*n*`` 
    where *name* is the actual name of the gate and *n* is the number of other 
    wires affected by this gate
    the *n* other wires will have as gate name
    ``@multi-ghost_*name*``

    Some examples :

    >>> gate_name_to_string('I')
    '--'
    >>> gate_name_to_string('I', True)
    '& \\qw '
    >>> gate_name_to_string('H', False)
    'H-'
    >>> gate_name_to_string('H', True)
    '& \\gate{H}'
    >>> gate_name_to_string('@multi-source_S@multi-size_1')
    'S-'
    >>> gate_name_to_string('@multi-source_S@multi-size_1', True)
    '& \\multigate{1}{S} '
    >>> gate_name_to_string('@multi-ghost_S')
    '|-'
    >>> gate_name_to_string('@multi-ghost_S', True)
    '& \\ghost{S} '
    """
    if _to_latex:
      if name == "I":
        return "& \\qw "
      elif "@multi-source_" in name:
        [source_name, size] = name.split("@multi-source_")[1].split("@multi-size_")
        return "& \\multigate{" + size + "}{" + source_name + "} "
      elif "@multi-ghost_" in name:
        source_name = name.split("@multi-ghost_")[1]
        return "& \\ghost{" + source_name + "} "
      else:
        return "& \\gate{" + name + "}"
    else:
      if name == "I":
        return "--"
      elif "@multi-source_" in name:
        source_name = name.split("@multi-source_")[1].split("@multi-size_")[0]
        return source_name + "-"
      elif "@multi-ghost_" in name:
        return "|-"
      else:
        return name + "-"

  def flatten_name_layers(_name_layers):
    r""" Flatten the list _name_layers to make it more easily usable

    >>> circuit = [[('I',2),('H',1)],[('H',1),('H',1),('I',1)],[('I',1),('S',2)]]
    >>> flatten_name_layers(circuit)
    [['I', 'I', 'H'],
     ['H', 'H', 'I'],
     ['I', '@multi-source_S@multi-size_1', '@multi-ghost_S']]
    """
    _name_layers_flattened = []
    for layer in _name_layers:
      layer_flattened = []
      for gate in layer:
        if gate[1] == 1:
          layer_flattened += [gate[0]]
        elif gate[0] == "I":
          layer_flattened += ["I"]*gate[1]
        else:
          layer_flattened += ["@multi-source_" + gate[0] + "@multi-size_" + str(gate[1]-1)]
          layer_flattened += ["@multi-ghost_" + gate[0]]*(gate[1]-1)
      _name_layers_flattened.append(layer_flattened)
    return _name_layers_flattened

  flattened_name_layers = flatten_name_layers(name_layers)
  n_wires = len(flattened_name_layers[0])
  n_layers = len(flattened_name_layers)
  if to_latex:
    circuit_string = "\\begin{align*}\n  \\Qcircuit @C=1em @R=.7em {\n"
  else:
    circuit_string = ""
  for n_wire in range(n_wires):
    if to_latex:
      wire_string = "    "
    else: 
      wire_string = "-"
    for n_layer in range(n_layers):
      wire_string += gate_name_to_string(flattened_name_layers[n_layer][n_wire], to_latex)
    if to_latex:
      circuit_string += wire_string + "& \\qw \\\\\n"
    else:
      circuit_string += wire_string + "\n"
  if to_latex:
    circuit_string += "  }\n\\end{align*}" 
  return circuit_string


def layers_to_printable(layers):
  r"""layers is the same as matrix_layers but with matrix names embedded in them
  
  Example:  
    >>> I2 = matrix.identity(2)
    >>> I4 = matrix.identity(4)
    >>> H = matrix([[1,1],[1,-1]])/sqrt(2)
    >>> swap = matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    >>> layers = [[('I',I4),('H',H)],[('H',H),('H',H),('I',I2)],[('I',I2),('S',swap)]]
    [[('I',2),('H',1)],[('H',1),('H',1),('I',1)],[('I',1),('S',2)]]

  :param list[list[(str,matrix)]] layers: Circuit under the format described 
    above.
  :returns: list[List[(str,int)]] -- Circuit under the ready-to-print format 
    described above.
  """
  return [[(gate[0], log(gate[1].nrows())/log(2)) for gate in layer] for layer in layers]


def layers_to_matrix(layers):
  r"""layers is the same as matrix_layers but with matrix names embedded in them
  
  Example:
    >>> I2 = matrix.identity(2)
    >>> I4 = matrix.identity(4)
    >>> H = matrix([[1,1],[1,-1]])/sqrt(2)
    >>> swap = matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    >>> layers = [[('I',I4),('H',H)],[('H',H),('H',H),('I',I2)],[('I',I2),('S',swap)]]
    >>> layers_to_matrix(layers)
    [[I4,H],[H,H,I2],[I2,swap]]

  :param list[list[(str,matrix)]] layers: Circuit under the format described 
    above.
  :returns: list[List[(str,int)]] -- Circuit under the ready-to-run format 
    described above.
  """
  return [[gate[1] for gate in layer] for layer in layers]
