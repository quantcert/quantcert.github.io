#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This module builds the QFT in Qiskit and runs it
"""
import csv
import numpy as np

from qiskit import *
from qiskit.circuit.library import Permutation

from .coefficients_shapes import *
from . import evaluation
from . import run


def build_QFT_0_to_k(nWires, k, measure=False):
    r""" Builds the QFT on nWires wires up to the `k^{th}` state generates.
    
    Note that the whole QFT can be generated in Qiskit using the `QFT` method
    
    :param int nWires: number of wires
    :param int k: number of gates in the output
    :returns: QuantumCirctuit -- `k^{th}` first gates of the QFT circuit
    """
    partial_QFT = QuantumCircuit(nWires)
    gate_count = 0
    if gate_count == k:
        return partial_QFT
    for wire in range(nWires):
        partial_QFT.h(wire)
        gate_count += 1
        if gate_count == k:
            return partial_QFT
        for inlayer_nb in range(2, nWires-(wire-1)):
            partial_QFT.cu1(np.pi/(2**inlayer_nb), wire, inlayer_nb+(wire-1))
            gate_count += 1
            if gate_count == k:
                return partial_QFT
        partial_QFT.barrier()
    SWAP = Permutation(nWires, [nWires-1-wire for wire in range(nWires)])
    meas = QuantumCircuit(nWires, nWires)
    meas.barrier(range(nWires))
    meas.measure(range(nWires),range(nWires))
    result = partial_QFT + SWAP + meas if measure else partial_QFT + SWAP
    return result


def QFT_length(nWires):
    """ Returns the length (number of gates) of the QFT.
    
    :param int nWires: Number of qbits on the system.
    :returns: int -- Number of gates in the QFT.
    """
    return nWires + int(nWires*(nWires-1)/2) + 1


def all_QFT_circuits(nWires):
    """ Returns a list of quantum circuit of the QFT built up to `k` gates, with 
        `k` varying between 0 and the full length of the QFT.
    
    :param int nWires: Number of qbits on the system.
    :returns: list[QuantumCircuit] -- A lits the partial QFTs of length varying 
        between 0 and the full length of the QFT.
    """
    return [build_QFT_0_to_k(nWires, partial_QFT_len) 
            for partial_QFT_len in range(QFT_length(nWires)+1)]


def periodic_state(l,r,nWires):
    r""" Returns the periodic state `|\varphi^{l,r}>` of size `2^{nWires}`. We 
        have:
        
    `|\varphi^{l,r}> = \sum_{i=0}^{A-1}|l+ir>/\sqrt(A)` with
    `A = floor((2^{nWires}-l)/r)+1`

    In this definition, ``l`` is the shift of the state, and ``r`` is the period 
    of the state.

    Example:
        Since
        `|\varphi^{1,5}> = (|1>+|6>+|11>)/\sqrt(3)=(|0001>+|0110>+|1011>)/\sqrt(3)`,

        >>> periodic_state(1,5,4)
        (0, 1/3*sqrt(3), 0, 0, 0, 0, 1/3*sqrt(3), 0, 0, 0, 0, 1/3*sqrt(3), 0, 0, 0, 0)

    :param int l: The shift of the state.
    :param int r: The period of the state.
    :param int nWires: The size of the system (number of qubits).
    :returns: vector -- The state defined by ``l``, ``r`` and ``nWires`` 
        according to the definition given above.
    """
    N = 2**nWires
    sqrt_A = np.sqrt(np.ceil((N-l)/r))
    result = [0]*N
    for k in range(int(np.ceil((N-l)/r))):
        result[l+k*r] = 1/sqrt_A
    return result


def get_coef_from_optimization_file(filename, iteration, evaluation=False):
    """ The file fed in this function must be a csv file with one of columns 
        being named "iteration", an other one being named "coefficients" and if
        the `evaluation` parameter is set to `True` a column named
        "intricationValue". The "iteration" column must contain integers, the
        "coefficients" column must contain tuples of real numbers (the mermin
        coefficients) and the "intricationValue" column must contain real
        numbers.

    Example:
        >>> get_coef_from_optimization_file("../examples/QFT-optimized-coefficients/1-1-4.csv",2)
        ([(-0.197738971530022, -0.00983193840670331,  0.980205403028067), 
          (0.892812904656093,  -0.035586934469795,    0.449019695976237), 
          (-0.892282320991669, -0.00788653439204609, -0.451408974457756), 
          (0.982839628418978,   0.012341254672589,   -0.184048793102134)], 
         [(0.430894968126063,  -0.0211632981654613,   0.902153890006799), 
          (-0.984337747324624,  0.0105235779205133,   0.175978559772592), 
          (0.984221659519729,  -0.0118117927364157,  -0.176545196719091), 
          (0.883187500168151,   0.0083188846423974,   0.468946303647911)])

    :param str filename: Name of the CSV file containing the information about 
        the Mermin coefficients.
    :param int iteration: Designates the line from which the data must be
        retrieved.
    :param bool intricationValue: If `True`, the evaluation will be returned as 
        well as the coefficients.
    :returns: tuple[list[tuple[real]]], real(optional) -- The mermin 
        coefficients previously optimized in packed shape, eventually with the
        optimum computed with these coefficients.
    """
    with open(filename, newline='') as coef_file:
        reader = csv.DictReader(coef_file, delimiter=',')
        for row in reader:
            if row['iteration'] == str(iteration):
                if evaluation:
                    return eval(row['coefficients']), eval(row['intricationValue'])
                else:
                    return eval(row['coefficients'])
    raise LookupError("Iteration \"" + str(iteration) + "\" not found in file \""\
     + filename + "\"")


def _main(l, r, nWires, optimization_filepath, local=True, simulation=True, epsilon=0.1):
    """ Performs and prints the Mermin evaluation on each state generated by the 
        QFT ran on IBM's quantum processor.

    Example:
        >>> QFT._main(1,1,4,"../examples/QFT-optimized-coefficients/1-1-4.csv")                        
        0.99609375
        1.04541015625
        1.0166015625
        1.0234375
        1.0224609375
        0.734375
        1.02685546875
        1.0556640625
        1.05419921875
        1.0361328125
        1.0263671875
        0.96826171875

    :param int l: Shift of the periodic state given as input to the QFT.
    :param int r: Period of the periodic state given as input to the QFT.
    :param int nWires: Number of qubits of the periodic state given as input to 
        the QFT.
    :param Path optimization_filepath: Path of the file containing the result on 
        the optimization process. This file has to be compatible with the
        function  ``get_coef_from_optimization_file``.
    :param bool local: If ``True`` the local simulator is used, otherwise, a 
        call to the IBMQ API will be made.
    :param bool simulator: If ``True``, the job will be run on a simulator, be 
        it  local or on IBM's servers, otherwise, it will be run on the quantum 
        processors made available by IBM.
    :param real epsilon: For now, only used with simulation, this is the 
        precision parameter used as a proof of concept for RAC in python.
    :returns: None
    """
    if not local:
        run.load_IBMQ_account()
        print("Account loaded")

    s = periodic_state(l,r,nWires)
    preparation_circuit = QuantumCircuit(nWires)
    preparation_circuit.initialize(s,range(nWires))

    for k in range(QFT_length(nWires)+1):
        coeffs, theoreticalValue = get_coef_from_optimization_file(
            optimization_filepath, k, evaluation=True)
        QFT_partial = build_QFT_0_to_k(nWires,k)
        value = evaluation.evaluate_polynomial(nWires, 
            preparation_circuit + QFT_partial, coeffs, is_simulation=simulation, 
            local=local, monitor=not simulation)
        if simulation:
            assert abs(theoreticalValue-value) < epsilon
        print(value)

    # In case of problem, retrieve results using backend.retrieve_job(job_id)