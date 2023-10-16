#!/usr/bin/env python
# coding: utf-8

# # This code plays the point-line Doily game using the Delegation method

# In[8]:


import numpy as np
import random
import itertools
from itertools import permutations
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel


# # Some functions

# In[12]:


# Provides the parity of a line for computation of third bit in context
def S_output(line):
    if line in neg_lines:
        return 1
    else:
        return 0

# Checks whether an operator is skew or symmetric
def Y_change(op):
    return op.count('Y')%2

# Performs a measurement on a qubit register based on associated operator in context
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

# Delegates the state of a qubit onto a delegation qubit based on associated operator
def Op_delegate_measurement(op, qc,q,qmeas):
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


# # Creating the Doily

# In[13]:


# Operators saved as strings
Op1 = 'IY'
Op2 = 'XY'
Op3 = 'XI'
Op4 = 'IX'
Op5 = 'XX'
Op6 = 'YY'
Op7 = 'ZZ'
Op8 = 'IZ'
Op9 = 'ZI'
Op10 = 'ZY'
Op11 = 'YZ'
Op12 = 'YX'
Op13 = 'ZX'
Op14 = 'YI'
Op15 = 'XZ'

Doily_points = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8, Op9, Op10, Op11, Op12, Op13, Op14, Op15]

line1 = [Op1, Op3, Op2]
line2 = [Op3, Op5, Op4]
line3 = [Op5, Op7, Op6] # negative
line4 = [Op7, Op9, Op8]
line5 = [Op9, Op1, Op10]
line6 = [Op6, Op14, Op1]
line7 = [Op2, Op12, Op7]
line8 = [Op8, Op15, Op3]
line9 = [Op4, Op13, Op9]
line10 = [Op10, Op11, Op5]
line11 = [Op11, Op2, Op13] # negative
line12 = [Op12, Op4, Op14]
line13 = [Op13, Op6, Op15]
line14 = [Op14, Op8, Op11]
line15 = [Op15, Op10, Op12] # negative

# Doily stored as an array of contexts
Doily = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15]

# The negative lines of this labelling
neg_lines = [line3, line11, line15]


# # Running the game on a quantum circuit

# In[20]:


is_simulation = True # Determines whether it's run on a simulator or real quantum computer
include_noise = False # Determines whether to include noise model

# Loading the IBM quantum experience account - include own credentials here
IBMQ.load_account()

# Determines project - include own info here
provider = IBMQ.get_provider(hub = 'hub', group = 'group', project = 'project')
quantum_comp=provider.get_backend('ibm_lagos')
noise_model = NoiseModel.from_backend(quantum_comp)

orders = list(permutations([0,1,2],3)) # permuting the operators along each line

orders_and_results = []
final_results = []

number_of_shots = 8192
simulator = Aer.get_backend('statevector_simulator')


# Iterating over all orders of the lines
for order in orders:
    doily_results = [] # storing results for a given order
    for doily_line in Doily: 
        Bob_line = [doily_line[i] for i in order]
        s_factor = S_output(doily_line)
        
        Alice_pos_results = [] # Storing the success rate for each of the 3 intersections of Bob's line
        for pos in range(3):
            Alice_point = Bob_line[pos]
            
            # Alice performs the skew-flip protocol using this variable
            Alice_y = Y_change(Alice_point)
        
            alice = [0,4] # Alice's qubit register, based on Lagos topoogy
            bob = [1,5] # Bob's qubit register, based on Lagos topoogy
            bob_del = 3
            
            # Creating the circuit and entangled state Psi
            circuit = QuantumCircuit(7,5)
            circuit.h(alice)
            circuit.cx(alice[0],bob[0])
            circuit.cx(alice[1],bob[1])

            # Measuring Alice's qubits on the first two registers 
            Op_measurement(Alice_point[0], circuit, alice[0], 0)
            Op_measurement(Alice_point[1], circuit, alice[1], 1)

            # Measuring Bob's first operator by delegating onto qubit bob_del
            Op_delegate_measurement(Bob_line[0][0], circuit, bob[0], bob_del)
            Op_delegate_measurement(Bob_line[0][1], circuit, bob[1], bob_del)
            circuit.measure(bob_del,2) # Measuring the delegation qubit
            
            # Measuring Bob's second operator directly on his register
            Op_measurement(Bob_line[1][0], circuit, bob[0], 3)
            Op_measurement(Bob_line[1][1], circuit, bob[1], 4)

            # Getting results
            if is_simulation:
                simulator = Aer.get_backend('qasm_simulator')
                if include_noise:
                    machine = "noisy simulator: "
                    result = execute(circuit,backend = simulator, 
                                     noise_model=noise_model, shots=number_of_shots).result()
                else:
                    machine = "noiseless simulator: "
                    result = execute(circuit,backend = simulator, shots=number_of_shots).result()
            else:
                machine = "IBM_Lagos"
                result = execute(circuit, backend=quantum_comp, shots=number_of_shots).result()
                
            measures_dictionary = result.get_counts()

            count = 0
            
            # Checking whether the players won, based on position of Alice's point
            for measure in measures_dictionary:
                if (pos==0) and (((measure[3]+measure[4]).count('1')+Alice_y)%2)==(((measure[2]).count('1'))%2):
                    count=count+measures_dictionary[measure]
                if (pos==1) and (((measure[3]+measure[4]).count('1')+Alice_y)%2)==(((measure[0]+measure[1]).count('1'))%2):
                    count=count+measures_dictionary[measure]
                if (pos==2) and (((measure[3]+measure[4]).count('1')+Alice_y)%2)==(((measure[0]+measure[1]+measure[2]).count('1')+s_factor)%2):
                    count=count+measures_dictionary[measure] 

            ordered_line_success = count/number_of_shots
            Alice_pos_results.append(ordered_line_success) # Storing successes along Bob's line
        
        avg_alice_order_success = np.mean(Alice_pos_results)
        doily_results.append(avg_alice_order_success)

    orders_and_results.append([order, doily_results]) # Storing all doily results and orderings
    final_results.append(np.mean(doily_results))

global_victory = max(final_results) # Taking the best performing order of the doily

print("Global victory for ", machine, global_victory)
        


# In[ ]:




