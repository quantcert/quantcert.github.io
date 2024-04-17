#!/usr/bin/env python
# coding: utf-8

# Copyright (C) 2024 Colm Kelleher.


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
            

# ----------------------------------------------------------------- #

# Below are some functions to check the success of a given game result

# ----------------------------------------------------------------- #


# Functions for converting results into bitstrings and check if game is won or not
def direct_measure_to_bit(bits):
	# converts triple of direct measured bits to one operator result
    return str((int(bits[0])+int(bits[1])+int(bits[2]))%2)

def bits_to_triplet(line, bits):
	# converts two bits to triplet answer using line parity
    s_factor = S_factor(line)
    return str(bits[0])+str(bits[1])+str((int(bits[0])+int(bits[1])+s_factor)%2)

def player_line_bits(line, bits):
	# returns triplet answer for bits
    reorg_bits = str(bits[0]) + str(bits[1])
    return bits_to_triplet(line, reorg_bits)

def coord(line, point):
	# gives coordinate of point in question on player line
    return line.index(point)

def coords(alice_line, bob_line, dylan_line, evan_line, point):
	# gives array of coordinates of point to check winning condition
    coords_arr = []
    for line in (alice_line, bob_line, dylan_line, evan_line):
        coords_arr.append(coord(line, point))
    return coords_arr

def intersection_check(alice_trip, bob_trip, dylan_trip, evan_trip, coords):
	# checks whether A*B*D*E = +1 condition, with A given by Alice measurement result on given point, etc.
    return not (int(alice_trip[coords[0]]) + int(bob_trip[coords[1]]) + int(dylan_trip[coords[2]]) + int(evan_trip[coords[3]]))%2



def results_check_single_del(point, alice_line, bob_line, daisy_line, evan_line, bit_result):
	# returns true if players have won the game or not
    coords_data = coords(alice_line, bob_line, daisy_line, evan_line, point)
    
    alice_bits = bit_result[0:3]
    bob_bits = bit_result[4:7]
    daisy_bits = bit_result[8:11]
    evan_bits = bit_result[12:15]
    
    alice_second_bit = direct_measure_to_bit(alice_bits)
    bob_second_bit = direct_measure_to_bit(bob_bits)
    daisy_second_bit = direct_measure_to_bit(daisy_bits)
    evan_second_bit = direct_measure_to_bit(evan_bits)
    
    alice_del_bit = bit_result[3]
    bob_del_bit = bit_result[7]
    daisy_del_bit = bit_result[11]
    evan_del_bit = bit_result[15]
    
    alice_trip = bits_to_triplet(alice_line, str(alice_del_bit)+str(alice_second_bit)) 
    bob_trip = bits_to_triplet(bob_line, str(bob_del_bit)+str(bob_second_bit))  
    daisy_trip = bits_to_triplet(daisy_line, str(daisy_del_bit)+str(daisy_second_bit))  
    evan_trip = bits_to_triplet(evan_line, str(evan_del_bit)+str(evan_second_bit))  
    
    return intersection_check(alice_trip, bob_trip, daisy_trip, evan_trip, coords_data)
    


# ----------------------------------------------------------------- #

# Below section implements Eloily as in 2-player game source code

# ----------------------------------------------------------------- #


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


# ----------------------------------------------------------------- #

# End of Eloily construction section

# ----------------------------------------------------------------- #



# Creating the 4-Player Eloily Game


# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #

# YOUR IBM QUANTUM EXPERIENCE API TOKEN HERE:
IBMProvider.save_account('<API_TOKEN>', overwrite=True)

# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #

# Load previously saved account credentials.
provider = IBMProvider()

# Select the backend and noise model
quantum_comp=provider.get_backend('ibm_brisbane')
noise_model = NoiseModel.from_backend(quantum_comp)

# Write column headers for the output date CSV
with open('4_player_results_file.csv', 'w', newline='') as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["Alice", "Bob", "Daisy", "Evan", "Measures_Dict", "Victory_checks", "Result"])
    
# Select whether to use a simulator and if so whether to include a noise model
is_simulation = True
is_noisy = False
number_of_shots = 8192# number of shots per game


choices = list(it.combinations(range(5),4)) # all possible 4-line choices at an intersection


# Select qubit assignments for players and their delegation qubits
alice = [0,1,2]
bob = [3,4,5]
daisy = [6,7,8]
evan = [9,10,11]

alice_del = [12]
bob_del = [13]
daisy_del = [14]
evan_del = [15]

Success_array = [] # store success rates

# loop over all points
for point in eloily_points:
    intersections = intersecting_lines(point)
    # loop over all 4-line choices per point
    for line_choices in choices:
        alice_line = intersections[line_choices[0]]
        bob_line = intersections[line_choices[1]]
        daisy_line = intersections[line_choices[2]]
        evan_line = intersections[line_choices[3]]

        # Making the circuit
        q = QuantumRegister(len(alice)+len(bob)+len(daisy)+len(evan)+len(alice_del)+len(bob_del)+len(daisy_del)+len(evan_del))
        c = ClassicalRegister(len(alice)+len(bob)+len(daisy)+len(evan)+len(alice_del)+len(bob_del)+len(daisy_del)+len(evan_del))
        circuit = QuantumCircuit(q,c)

        # Entangled state, based on stabiliser state:
        circuit.h (alice[0]);
        circuit.cx (alice[0],bob[0]);
        circuit.cx (alice[0],daisy[0]);
        circuit.cx (alice[0],evan[0]);

        circuit.h (alice[1]);
        circuit.cx (alice[1],bob[1]);
        circuit.cx (alice[1],daisy[1]);
        circuit.cx (alice[1],evan[1]);

        circuit.h (alice[2]);
        circuit.cx (alice[2],bob[2]);
        circuit.cx (alice[2],daisy[2]);
        circuit.cx (alice[2],evan[2]);


        # Delegating players first operator
        Op_delegate(alice_line[0][0], circuit, alice[0], alice_del[0])
        Op_delegate(alice_line[0][1], circuit, alice[1], alice_del[0])
        Op_delegate(alice_line[0][2], circuit, alice[2], alice_del[0])


        Op_delegate(bob_line[0][0], circuit, bob[0], bob_del[0])
        Op_delegate(bob_line[0][1], circuit, bob[1], bob_del[0])
        Op_delegate(bob_line[0][2], circuit, bob[2], bob_del[0])


        Op_delegate(daisy_line[0][0], circuit, daisy[0], daisy_del[0])
        Op_delegate(daisy_line[0][1], circuit, daisy[1], daisy_del[0])
        Op_delegate(daisy_line[0][2], circuit, daisy[2], daisy_del[0])


        Op_delegate(evan_line[0][0], circuit, evan[0], evan_del[0])
        Op_delegate(evan_line[0][1], circuit, evan[1], evan_del[0])
        Op_delegate(evan_line[0][2], circuit, evan[2], evan_del[0])

        
        # Measuring players's register qubits
        Op_measurement(alice_line[1][0], circuit, alice[0], 0)
        Op_measurement(alice_line[1][1], circuit, alice[1], 1)
        Op_measurement(alice_line[1][2], circuit, alice[2], 2)
        
        Op_measurement(bob_line[1][0], circuit, bob[0], 4)
        Op_measurement(bob_line[1][1], circuit, bob[1], 5)
        Op_measurement(bob_line[1][2], circuit, bob[2], 6)
        
        Op_measurement(daisy_line[1][0], circuit, daisy[0], 8)
        Op_measurement(daisy_line[1][1], circuit, daisy[1], 9)
        Op_measurement(daisy_line[1][2], circuit, daisy[2], 10)
        
        Op_measurement(evan_line[1][0], circuit, evan[0], 12)
        Op_measurement(evan_line[1][1], circuit, evan[1], 13)
        Op_measurement(evan_line[1][2], circuit, evan[2], 14)

        # Measuring players's delegation qubits
        circuit.measure(alice_del[0], 3)
        circuit.measure(bob_del[0], 7)
        circuit.measure(daisy_del[0], 11)
        circuit.measure(evan_del[0], 15)

        if is_simulation:
        	local_simulator = Aer.get_backend('qasm_simulator')
        	if is_noisy:
        		result = execute(circuit, backend=local_simulator, noise_model = noise_model, shots=number_of_shots).result()
        	else:
        		result = execute(circuit, backend=local_simulator, shots=number_of_shots).result()
        else:
            result = execute(circuit, backend=quantum_comp, shots=number_of_shots).result()
        
        # Getting results
        qubit_results = [key[::-1] for key in result.get_counts(circuit)]
        qubit_counts = list(result.get_counts(circuit).values())
        total_counts = sum(qubit_counts)

		# Checking to see if each game was victorious or not
        game_checks = [results_check_single_del(point, alice_line, bob_line, daisy_line, evan_line, bit_result) for bit_result in qubit_results]
        game_coords = []
        for i in range(len(game_checks)):
            if game_checks[i]:
                game_coords.append(i)

		# Compute success rate
        winning_count = sum([game_checks[i]*qubit_counts[i] for i in range(len(qubit_counts))])/number_of_shots
        Success_array.append(winning_count)
    
        if not is_simulation:
            with open('4_player_results_file.csv', 'w') as outfile:
                writer = csv.writer(outfile)
                for row in zip(point, alice_line, bob_line, daisy_line, evan_line, result.get_counts(), game_checks, winning_count):
                    writer.writerow(row)
        

victory = np.mean(Success_array)
print("Victory rate: ",victory)




