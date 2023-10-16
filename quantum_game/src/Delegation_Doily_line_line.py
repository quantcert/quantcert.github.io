#!/usr/bin/env python
# coding: utf-8

# # This code plays the line-line Doily game using the Delegation method

# In[1]:


import numpy as np
import random
import itertools as it
from itertools import permutations
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel


# In[2]:


# Creating 2-qubit operators
X=np.array([[0,1],[1,0]]) 
Z=np.array([[1,0],[0,-1]])
Y=np.array([[0,-1j],[1j,0]])
I = np.identity(2)


# # Creating dictionary for operators and their strings

# In[3]:


# Creating the set of all nontrivial 2-qubit operators in string form (i.e. the points in the Doily)
Ops_str = ["I","X","Y","Z"]
Doily_points = []
for i in Ops_str:
    for j in Ops_str:
        double_op_str = i+j
        Doily_points.append(double_op_str)
Doily_points.remove('II')




# # Some functions

# In[4]:


# Finding all lines for a given intersection point and all transform matrices for those lines
def intersecting_lines(Op):
    Lines = [line for line in Doily if Op in line]
    comb_list = list(it.combinations(Lines,2))
    return comb_list

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

# In[5]:


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

# In[7]:


is_simulation = True # Determines whether it's run on a simulator or real quantum computer
include_noise = False # Determines whether to include noise model

# Loading the IBM quantum experience account - include own credentials here
IBMQ.load_account()

# Determines project - include own info here
provider = IBMQ.get_provider(hub = 'hub', group = 'group', project = 'project')
quantum_comp=provider.get_backend('ibm_lagos')
noise_model = NoiseModel.from_backend(quantum_comp)

orders = list(permutations([0,1,2],3)) # permuting the operators along each line

orders_results = []
final_results = []

number_of_shots = 8192
simulator = Aer.get_backend('statevector_simulator')


# Iterating over all orders of the lines
for order in orders:
    particular_order_success = []    
    for point in Doily_points: # Cycling through all points of the doily
        point_intersecting_lines = intersecting_lines(point) # Getting all lines through that point
        for intersections in point_intersecting_lines:
            Alice_line = [intersections[0][i] for i in order]
            Bob_line = [intersections[1][i] for i in order]
            s_a = S_output(Alice_line) # lines parities for each player's lines
            s_b = S_output(Bob_line)
            Alice_y = [str(Alice_line[k]).count('Y')%2 for k in range(3)] # Alice performs skew-flip protocol
            row, col = Alice_line.index(point), Bob_line.index(point) # row, column equivalent for checking success
        
            alice = [5,1] # Alice's qubit register, based on Lagos topoogy
            bob = [3,2] # Bob's qubit register, based on Lagos topoogy

            d_alice = [6] # Delegation qubits
            d_bob = [0]
            
            # Creating the circuit and entangled state Psi
            circuit = QuantumCircuit(7,6)
            circuit.h(alice)
            circuit.cx(alice[0],bob[0])
            circuit.cx(alice[1],bob[1])

            # Delegating Alice's qubits to her delegation qubits 
            Op_delegate_measurement(Alice_line[0][0], circuit, alice[0], d_alice[0])
            Op_delegate_measurement(Alice_line[0][1], circuit, alice[1], d_alice[0])

            circuit.measure(d_alice[0],0)
            Op_measurement(Alice_line[1][0], circuit, alice[0], 1)
            Op_measurement(Alice_line[1][1], circuit, alice[1], 2)

            # Delegating Bob's qubits to his delegation qubits
            Op_delegate_measurement(Bob_line[0][0], circuit, bob[0], d_bob[0])
            Op_delegate_measurement(Bob_line[0][1], circuit, bob[1], d_bob[0])

            circuit.measure(d_bob[0],3)
            Op_measurement(Bob_line[1][0], circuit, bob[0], 4)
            Op_measurement(Bob_line[1][1], circuit, bob[1], 5)
            
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
            qubit_results = [key for key in measures_dictionary]
            qubit_counts = list(measures_dictionary.values())
            total_counts = sum(qubit_counts)

            # Getting player's answers as bit strings
            Alice_results = [qubits[3:6][::-1] for qubits in qubit_results]
            Bob_results = [qubits[0:3][::-1] for qubits in qubit_results]
            
            # Alice's answer before the skew-flip protocol
            Alice_ans_pre_skew_flip = [qubits[0]+str((int(qubits[1])+int(qubits[2]))%2)+str((int(qubits[0])+int(qubits[1])+int(qubits[2])+s_a)%2) for qubits in Alice_results]
            Alice_ans_y = [] # Storing Alice's answer after the skew-flip protocol
            for result in Alice_ans_pre_skew_flip:
                new_result = str()
                for i in range(3):
                    new_result_digit = str((int(result[i])+Alice_y[i])%2)
                    new_result += new_result_digit
                Alice_ans_y.append(new_result)
                
            # Bob's answer including third computed bit
            Bob_ans = [qubits[0]+str((int(qubits[1])+int(qubits[2]))%2)+str((int(qubits[0])+int(qubits[1])+int(qubits[2])+s_b)%2) for qubits in Bob_results]

            # Checking success for this line-line combo
            checklist = [Alice_ans_y[i][row] == Bob_ans[i][col] for i in range(len(Alice_ans_y))]
            true_counts = [checklist[i]*qubit_counts[i] for i in range(len(qubit_results))]
            success_rate = sum(true_counts)/float(total_counts)
            particular_order_success.append(success_rate)
    
    # Getting all success rates over all orders
    orders_results.append(sum(particular_order_success)/len(particular_order_success)) 
        
global_victory = max(orders_results) # Taking the best performing order of the doily

print("Global victory for ", machine, global_victory)
        


# In[ ]:




