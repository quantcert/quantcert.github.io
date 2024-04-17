#!/usr/bin/env python
# coding: utf-8


# NOTE: Code is written for Qiskit 0.46.0. 

import numpy as np
from collections import Counter
import itertools as it
import math
import csv

from qiskit import *
from qiskit_aer.noise import NoiseModel
from qiskit_aer import Aer 
from qiskit_ibm_provider import IBMProvider



# Function to give tensor product of 3 qubit operators
def tensor3(mat_1, mat_2, mat_3):
    return np.kron(mat_1, np.kron(mat_2, mat_3))

# Function to "multiply" string representation of 1-qubit operators, up to phase
def str_mult(op_par):
    single_op1, single_op2 = op_par
    if single_op1 == single_op2:
        return 'I'
    elif single_op1 == 'I':
        return single_op2
    elif single_op2 == 'I':
        return single_op1
    else:
        op_set = (single_op1, single_op2)
        if op_set in [('X','Y') , ('Y','X')]:
            return 'Z'
        elif op_set in [('X','Z') , ('Z','X')]:
            return 'Y'
        elif op_set in [('Z','Y') , ('Y','Z')]:
            return 'X'

# Function to give string product of 3-qubit operators        
def str_op_product(op1, op2):
    first_ops = (op1[0], op2[0])
    second_ops = (op1[1], op2[1])
    third_ops = (op1[2], op2[2])
    return str(str_mult(first_ops))+str(str_mult(second_ops))+str(str_mult(third_ops))

# Function to give product of operators across context in string form
def context_str_prod(context):
    return np.matmul(np.matmul(Ops[context[0]], Ops[context[1]]),Ops[context[2]])

# Function to give product of operators across context in ops form
def context_ops_prod(context):
    return np.matmul(np.matmul(context[0], context[1]),context[2])

# Function to check sign of overall context
def id_check(M):
    if (M.shape[0] == M.shape[1]) and (M == np.eye(M.shape[0])).all():
        return 1
    elif (M.shape[0] == M.shape[1]) and (M == -np.eye(M.shape[0])).all():
        return -1
    else:
        return 0
    
# Normalises a "unitary" matrix of the form A*A^+ = cI for some scalar c.
def normalise(mat):
    matrix = sp.Matrix(mat)
    norm = sp.sqrt(np.matmul(np.conj(matrix).T,matrix)[0,0])
    if norm == 0: 
       return matrix
    return matrix / norm

def str_to_array(string):
    arr = []
    for i in range(len(string)):
        arr.append(int(string[i]))
    return arr
    
    
# Measurement Functions

# Gives the parity of the line (1 for negative)
def S_factor(line):
    if line in neg_lines:
        return 1
    return 0

# Flips the player result if it's skew symmetric
def Y_change(op):
    return op.count('Y')%2

# Implements basis change and measurement gates for single qubit operator
def Op_measurement(op, qc, q, c):
    if op == 'X':
        qc.h(q)
        qc.measure(q,c)
    elif op == 'Z':
        qc.measure(q,c)
    elif op == 'Y':
        qc.sdg(q)
        qc.h(q)
        qc.measure(q,c)

# Implements basis change and delegation gates for single qubit operator
def Op_delegate(op, qc,q,qmeas):
    if op == 'X':
        qc.h(q)
        qc.cx(q,qmeas)
        qc.h(q)
    elif op == 'Z':
        qc.cx(q,qmeas)
    elif op == 'Y':
        qc.sdg(q)
        qc.h(q)
        qc.cx(q,qmeas)
        qc.h(q)
        qc.s(q)


# Finding a stabiliser state for 3 qubits

# Projects from an operator to a point in PG(5,2)
def project(point):
    string = ''
    for op in point:
        if op == 'X':
            string = string + '10'
        elif op == 'Z':
            string = string + '01'
        if op == 'Y':
            string = string + '11'
        elif op == 'I':
            string = string + '00'
    return string

def string_to_col(point):
    string = point + point
    binary_col = project(string)
    return str_to_array(binary_col)


# function to return stabiliser state for 3 qubits
def stabiliser_state(test_chi):
    
    # Forming tensor products of all Gamma generators
    O1 = np.kron(Ops[gamma_str[0]],Ops[gamma_str[0]])
    O2 = np.kron(Ops[gamma_str[1]],Ops[gamma_str[1]])
    O3 = np.kron(Ops[gamma_str[2]],Ops[gamma_str[2]])
    O4 = np.kron(Ops[gamma_str[3]],Ops[gamma_str[3]])
    O5 = np.kron(Ops[gamma_str[4]],Ops[gamma_str[4]])
    O6 = np.kron(Ops[gamma_str[5]],Ops[gamma_str[5]])
    O7 = np.kron(Ops[gamma_str[6]],Ops[gamma_str[6]])

    # New 6-qubit generators
    Ops_gen = [O1, O2, O3 ,O4, O5, O6]

    Chi_bin = format(test_chi, '064b') # test_chi input determines initial state Chi in binary form

    Chi = str_to_array(Chi_bin)
    Id = np.eye(64)
    Psi = Chi # start stabiliser state as initial state Chi
    Op_product = Id
    for op in Ops_gen: # iterate over product of identity+generators
        Op_product = np.matmul(0.5*(Id+op), Op_product)
        Psi = np.matmul(0.5*(Id+op), Psi).real
    
    # returning stabiliser state in array form with binary-encoded vectors (works for balanced state)
    psi_arr = []
    for term in Psi:
        if term != 0: 
            psi_arr.append(str(np.sign(term))+str(format(Psi.index(term), '06b')))
    
    return psi_arr

# The State: 
# -000111+001110+010101-011100+100011-101010-110001+111000


## Some functions for creating the game conditions

# Finds all intersecting lines with a given point
def intersecting_lines(point):
    intersecting_lines_arr = []
    for line in all_lines:
        if point in line:
            intersecting_lines_arr.append(line)
    return intersecting_lines_arr

# Returns the unique negative line in the canonical Eloily that a given point is on
def negative_line(point):
    for line in neg_lines:
        if point in line:
            return line



# Defining the Pauli operators

I = [[1, 0],
    [0, 1]]

X = [[0, 1],
    [1, 0]]

Y = [[0, -1.j],
    [1.j, 0]]

Z = [[1.0, 0],
    [0, -1.0]]


# Setting up dictionary of 3 qubit operators (strings <-> matrix reps)

# Creating the set of all three qubit operators

Ops_arr = [I,X,Y,Z]
Triple_ops = []
for i in Ops_arr:
    for j in Ops_arr:
        for k in Ops_arr:
            triple_op = tensor3(i,j,k)
            Triple_ops.append(triple_op)
Triple_ops.remove(Triple_ops[0])

Ops_str = ["I","X","Y","Z"]
Triple_ops_str = []
for i in Ops_str:
    for j in Ops_str:
        for k in Ops_str:
            triple_op = i+j+k
            Triple_ops_str.append(triple_op)
Triple_ops_str.remove('III')

Ops = {}
for i in range(len(Triple_ops)):
    Ops[Triple_ops_str[i]] = Triple_ops[i]

# Creating all positive points in GQ(2,4) (i.e. ignoring phase, anticommutation of gammas)

# Creating the generators of W(5,2) (Gamma_ij operators)

gamma_str = ['ZYI', 'YIX', 'XYI', 'IXY', 'YIZ', 'IZY', 'YYY']

# Creating the contexts of the canonical doily in terms of gamma indices i,j
# here gamma_ij := gamma_i * gamma_j, with gamma_i, gamma_j the i, jth entries in gamma_str
doily_indices = [[(1,2),(3,4),(5,6)],
                 [(1,2),(3,5),(4,6)],
                 [(1,2),(3,6),(4,5)],
                 [(1,3),(2,4),(5,6)],
                 [(1,3),(2,5),(4,6)],
                 [(1,3),(2,6),(4,5)],
                 [(1,4),(2,3),(5,6)],
                 [(1,4),(3,5),(2,6)],
                 [(1,4),(3,6),(2,5)],
                 [(1,5),(3,4),(2,6)],
                 [(1,5),(2,3),(4,6)],
                 [(1,5),(3,6),(2,4)],
                 [(1,6),(3,4),(2,5)],
                 [(1,6),(3,5),(2,4)],
                 [(1,6),(2,3),(4,5)]]

# recreating the above doily but mapping gamma indices to 3-qubit operators
doily_contexts_str = []
for indices in doily_indices:
    context_str = []
    for pair in range(3):
        context_str.append(str_op_product(str_op_product(gamma_str[indices[pair][0]-1], gamma_str[indices[pair][1]-1]),gamma_str[6]))
    doily_contexts_str.append(context_str)

# Adding the sign parity of each context
doily_contexts_with_signs = [[context, id_check(context_str_prod(context))] for context in doily_contexts_str]
# print("Doily contexts with signs:")
# print(doily_contexts_with_signs)
                            

# Creating the canonical double-six contexts also from the gamma generators

six_one = [gamma for gamma in gamma_str[0:6]]
six_two = [str_op_product(gamma, gamma_str[6]) for gamma in gamma_str[0:6]]
six_two_ops = [np.matmul(Ops[gamma], Ops[gamma_str[6]]*(-1.j)) for gamma in gamma_str[0:6]]


double_six_contexts_str = []
for i in range(6):
    for j in range(6):
        if i != j:
            context_str = [six_one[i], six_two[j], str_op_product(six_one[i], six_two[j])]
            double_six_contexts_str.append(context_str)
            
# Adding context sign parities
double_six_contexts_str_with_signs = [[context, id_check(context_str_prod(context))] for context in double_six_contexts_str]
# print("Double six contexts:")
# print((double_six_contexts_str_with_signs))


# Combining doily and double six to get eloily
all_contexts = doily_contexts_with_signs + double_six_contexts_str_with_signs

# Identifying all negative lines
neg_lines = []
for context in all_contexts:
    if context[1] == -1:
        neg_lines.append(context[0])
        

# Forming array of just the contexts (no parities)
all_lines = []
for context in all_contexts:
    all_lines.append(context[0])
points_on_all_lines = np.array(all_lines).flatten()
eloily_points = [item for item in Counter(points_on_all_lines).keys()]




# Creating the Eloily Game


# Implementing the game in a circuit 

# ------------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------------#


# YOUR IBM QUANTUM EXPERIENCE API TOKEN HERE:
IBMProvider.save_account('<API_TOKEN>', overwrite=True)


# ------------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------------#


# Load previously saved account credentials.
provider = IBMProvider()

# Select the backend and noise model
quantum_comp=provider.get_backend('ibm_brisbane')
noise_model = NoiseModel.from_backend(quantum_comp)

# Write column headers for the output date CSV
with open('results_file.csv', 'w', newline='') as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["Alice", "Bob", "Measures_Dict", "Bitstring_Answers", "Result"])
    

# Set whether doing simulation and if so whether using noise model
is_simulation = True
is_noisy = False
number_of_shots = 8192 # number of shots per game


# Allocating the qubits based on the CNOT connection topology

del_A, del_B = 28, 32 # Alice and Bob's delegation qubits
A1, B1 = 29,31 # Alice & Bobs first qubits A1, B1
A2, B2 = 35,36 
A3, B3 = 48,50
# If simulator doesn't accept such large qubit assignments change to any lower ones for simulation



# listing ways of choosing 2 lines from 5 (each point lies on 5 lines)
intersection_options = list(it.combinations(range(5),2))

Success_arr = []
# Iterating over all points in Eloily
for point in eloily_points:
    point_intersecting_lines = intersecting_lines(point) # all lines intersecting at that point
    
    # Iterating over all pairs of intersecting lines
    for a_ind, b_ind in intersection_options:
        Alice_line = point_intersecting_lines[a_ind]
        Bob_line = point_intersecting_lines[b_ind]
        
        # Line parities
        s_factor_a = S_factor(Alice_line)
        s_factor_b = S_factor(Bob_line)

        col = Alice_line.index(point) # defines "column" equivalent in Eloily for Alice's point
        row = Bob_line.index(point) # defines "row" equivalent for Bob's point      
        
        # Making the circuit
        q = QuantumRegister(max([del_A, del_B, A1, B1, A2, B2, A3, B3])+1) # ensuring we have enough qubits based on allocation choice
        c = ClassicalRegister(8)
        circuit = QuantumCircuit(q,c)

        # Applying gates to implement the stabiliser state:
        circuit.h (A1);
        circuit.h (A2);
        circuit.h (A3);
        circuit.cx (A1,B1);
        circuit.cx (A2,B2);
        circuit.cx (A3,B3);
        circuit.x (B1);
        circuit.x (B2);
        circuit.x (B3);
        circuit.z (A1);
        circuit.z (B2);
        circuit.z (B3);

        # Delegating and measuring Alice's qubits for first operator
        Op_delegate(Alice_line[0][0], circuit, A1, del_A)
        Op_delegate(Alice_line[0][1], circuit, A2, del_A)
        Op_delegate(Alice_line[0][2], circuit, A3, del_A)
        
        circuit.measure(del_A,c[0])

        # Delegating and measuring Bob's qubits for first operator
        Op_delegate(Bob_line[0][0], circuit, B1, del_B)
        Op_delegate(Bob_line[0][1], circuit, B2, del_B)
        Op_delegate(Bob_line[0][2], circuit, B3, del_B)
        
        circuit.measure(del_B,c[4])

        # Performing direct measurements on the original qubits
        Op_measurement(Alice_line[1][0], circuit, A1, c[1])
        Op_measurement(Alice_line[1][1], circuit, A2, c[2])
        Op_measurement(Alice_line[1][2], circuit, A3, c[3])
        
        Op_measurement(Bob_line[1][0], circuit, B1, c[5])
        Op_measurement(Bob_line[1][1], circuit, B2, c[6])
        Op_measurement(Bob_line[1][2], circuit, B3, c[7])

        # Setting up run on simulator or real QC 
        if is_simulation:
            local_simulator = Aer.get_backend('statevector_simulator')
#             local_simulator = AerSimulator()
            if is_noisy:
                result = execute(circuit, backend=local_simulator, noise_model = noise_model, shots=number_of_shots).result()
            else:
                result = execute(circuit, backend=local_simulator, shots=number_of_shots).result()
        else:
            result = execute(circuit, backend=quantum_comp, shots=number_of_shots).result()
        measures_dictionary = result.get_counts()

        # Getting results, noting IBM returns bitstrings backwards
        qubit_results = [key[::-1] for key in result.get_counts(circuit)]
        qubit_counts = list(result.get_counts(circuit).values())
        total_counts = sum(qubit_counts)

        # Configuring player results
        Alice_results = [qubits[0:4] for qubits in qubit_results]
        Bob_results = [qubits[4:8] for qubits in qubit_results]

        # Converts players's 4 bits to a 3-bit answer, one per qubit measurement
        Alice_ans = [str(int(qubits[0]))
                     +str((int(qubits[1])+int(qubits[2])+int(qubits[3]))%2)
                     +str((int(qubits[0])+int(qubits[1])+int(qubits[2])
                           +int(qubits[3])+s_factor_a)%2) 
                     for qubits in Alice_results]
        Bob_ans = [str(int(qubits[0]))
                     +str((int(qubits[1])+int(qubits[2])+int(qubits[3]))%2)
                     +str((int(qubits[0])+int(qubits[1])+int(qubits[2])
                           +int(qubits[3])+s_factor_b)%2) 
                   for qubits in Bob_results]

        # Checks to see for each shot whether they win the game
        checklist = [int(Alice_ans[i][col]) == int(Bob_ans[i][row]) for i in range(len(Alice_ans))]

        # Computes success rate for game
        true_counts = [checklist[i]*qubit_counts[i] for i in range(len(qubit_results))]
        success_rate = sum(true_counts)/float(total_counts)
        
        # Writes result to output file
        if not is_simulation:
            with open("results_file.csv",'a') as file:
                    writer = csv.writer(outcsv)
                    writer.writerow([Alice_line, Bob_line, measures_dictionary, Alice_ans+Bob_ans, success_rate])
        Success_arr.append(success_rate)
        

print("Success rate for each intersection point: ")
print("-"*30)
print("Global victory for Brisbane: ", np.mean(Success_arr))


