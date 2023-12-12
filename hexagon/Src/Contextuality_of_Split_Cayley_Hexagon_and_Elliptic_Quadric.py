#!/usr/bin/env python
# coding: utf-8

# In[72]:


# Importing necessary libraries
import numpy as np
from qiskit import *
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit import QuantumCircuit, transpile
import sympy as sp


# In[54]:


# Defining the spin operator matrices

I = [[1, 0],
    [0, 1]]

X = [[0, 1],
    [1, 0]]

Y = [[0, -1.j],
    [1.j, 0]]

Z = [[1.0, 0],
    [0, -1.0]]


# # Creating the dictionary of all three qubit operators

# In[55]:


def tensor3(A,B,C): # tensor product of 3 operators
    return np.kron(A, np.kron(B,C))

Ops_arr = [I,X,Y,Z] # Matrix representation of operators
Triple_ops = []
for i in Ops_arr:
    for j in Ops_arr:
        for k in Ops_arr:
            triple_op = tensor3(i,j,k)
            Triple_ops.append(triple_op)
Triple_ops.remove(Triple_ops[0]) # Removing trivial operator III

Ops_str = ["I","X","Y","Z"] # String representation of operators
Triple_ops_str = []
for i in Ops_str:
    for j in Ops_str:
        for k in Ops_str:
            triple_op = i+j+k
            Triple_ops_str.append(triple_op)
Triple_ops_str.remove('III') # Removing trivial operator

# Creating a dictionary between representations
Ops = {}
for i in range(len(Triple_ops)):
    Ops[Triple_ops_str[i]] = Triple_ops[i]


# # Functions for finding 1-qubit operator products based on string representations

# In[56]:


def char_pos(indA, indB): # takes as input indices of 2 different nontrivial operators in [I,X,Z,Y] and returns index of product
    return (indA*indB + indA + indB - 3 + abs(indA-indB))%4

def char_phase(indA, indB): # gives the phase of the operator product given the inputs (e.g. (X, Y) -> iZ, gives i)
    return (np.sign(indB-indA)*np.power(-1,(indA*indB)%2))*(1.0j)*1

def char_product_phase(char_a, char_b): # returns string value for product of 1-pauli operators, as well as phase
    chars = ('I', 'X', 'Y', 'Z')
    indA, indB = chars.index(char_a), chars.index(char_b)
    if (indB - indA) == 0: # checks to see if the operators are the same
        return 'I', 1
    elif (indB*indA) == 0: # checks to see if one is the identity
        return chars[max(indA, indB)], 1
    else: # cyclic cases, e.g. (Y, X) -> -iZ
        return chars[char_pos(indA, indB)], char_phase(indA, indB)

def three_op_prod(op1, op2): # gives the 3-qubit pauli operator product of the two inputs, and the overall phase
    prod1, prod2, prod3 = char_product_phase(op1[0],op2[0]), char_product_phase(op1[1],op2[1]), char_product_phase(op1[2],op2[2])
    return prod1[0]+prod2[0]+prod3[0], prod1[1]*prod2[1]*prod3[1]

print(three_op_prod('XIX','YZI')[0])


# # Functions for creating PG(6,2), the parabolic quadric Q, and the Plucker equations

# In[57]:


def single_proj(op): # projects a 2-qubit operator into a vector over F_2
    return format(['I','X','Z','Y'].index(op), '02b')

def Proj(op_str): # projects an N-qubit operator into PG(2N-1,2)
    x = ''
    for op in op_str:
        x += single_proj(op)
    return [int(i) for i in x]

def vect_to_operator(bin_num): # Takes a vector in PG(5,2) and returns the corresponding 3-qubit operator
    op1_b = int(bin_num[0]+bin_num[3],2)
    op2_b = int(bin_num[1]+bin_num[4],2)
    op3_b = int(bin_num[2]+bin_num[5],2)
    ops = ['I','X','Z','Y']
    return ops[op1_b]+ops[op2_b]+ops[op3_b]

def p_q_comm(p,q): # The symplectic form giving commutation relations
    sum_part = 1 # Here we define it so <p,q> = 1 iff operators commute, for boolean checks
    N = int(len(p)/2)
    for i in range(N):
        sum_part += p[i]*q[N+i] + p[N+i]*q[i]
    return (sum_part)%2

def p_ij(x,y,i,j): # line definition used in Plucker equations
    return (x[i]*y[j] - x[j]*y[i])%2

def Extend(x): # Extends a vector in PG(5,2) to PG(6,2) via parabolic quadrant
    x.append((x[0]*x[3]+x[1]*x[4]+x[2]*x[5])%2)
    return x

def array_maker(bin_num): # Turns a binary string into an array for easier parsing
    return[int(i) for i in bin_num]

def Plucker_check_7(x,y): # Set of Plucker equations giving 63 lines in PG(6,2)
    
    # Plucker equations:
    P_1 = (p_ij(x,y,5,1) == p_ij(x,y,0,6))
    P_2 = (p_ij(x,y,0,2) == p_ij(x,y,6,1))
    P_3 = (p_ij(x,y,1,3) == p_ij(x,y,2,6))
    P_4 = (p_ij(x,y,2,4) == p_ij(x,y,6,3))
    P_5 = (p_ij(x,y,3,5) == p_ij(x,y,4,6))
    P_6 = (p_ij(x,y,4,0) == p_ij(x,y,6,5))
    P_7 = ((p_ij(x,y,0,3) + p_ij(x,y,1,4) + p_ij(x,y,2,5))%2 == 0)
    
    return P_1*P_2*P_3*P_4*P_5*P_6*P_7


def f4(x): # first mapping used for skew embedding
    return x[2]*x[4] + x[6]*x[3]

def f5(x): # second mapping used for skew embedding
    return x[3]*x[5] + x[6]*x[4]

def skew_map(x): # turns a point in the classical embedding into a point in the skew embedding
    return [x[0]+x[5] + f5(x), x[1] + x[2] + f4(x), x[2], x[3],x[4], x[5],x[6]]

Q_6 = [format(i,'06b') for i in range(1,64)] # all vectors in PG(5,2)


# # Building W(5,2)

# In[58]:


W_5_2 = []
for p_str in Q_6:
    p_lines = []
    for q_str in [string for string in Q_6 if string != p_str]:
        p_vec = Extend(array_maker(p_str)) # takes all binary strings in PG(5,2), turns them to vectors and extends to PG(6,2)
        q_vec = Extend(array_maker(q_str))
        if p_q_comm(p_vec[0:6],q_vec[0:6]): # checks only to see if corresponding operators commute
            op3_with_phase = three_op_prod(vect_to_operator(p_str), vect_to_operator(q_str)) # finds the third operator in the line
            p_lines.append([set([vect_to_operator(p_str), vect_to_operator(q_str), op3_with_phase[0]]), op3_with_phase[1]]) # last item gives sign for overall line
    [W_5_2.append(line) for line in p_lines if line not in W_5_2] # gets unique lines (i.e. ignores ordering)



# # Building the Classical Split-Cayley Hexagon

# In[59]:


SCH_Classical = []
for p_str in Q_6:
    p_lines = []
    for q_str in [string for string in Q_6 if string != p_str]:
        p_vec = Extend(array_maker(p_str))
        q_vec = Extend(array_maker(q_str))
        if Plucker_check_7(p_vec,q_vec) and p_q_comm(p_vec[0:6],q_vec[0:6]): # Same as for W(5,2) creation but also applies Plucker equations for only 63 lines
            op3_with_phase = three_op_prod(vect_to_operator(p_str), vect_to_operator(q_str))
            p_lines.append([set([vect_to_operator(p_str), vect_to_operator(q_str), op3_with_phase[0]]), op3_with_phase[1]]) 
    [SCH_Classical.append(line) for line in p_lines if line not in SCH_Classical] # gets unique lines (i.e. ignores ordering)



# # Building the Skew Split-Cayley Hexagon

# In[60]:


SCH_Skew = []
for p_str in Q_6:
    p_lines = []
    for q_str in [string for string in Q_6 if string != p_str]:
        p_vec = skew_map(Extend(array_maker(p_str))) # maps the points to the points in the skew embedding
        q_vec = skew_map(Extend(array_maker(q_str))) # maps the points to the points in the skew embedding
        if Plucker_check_7(p_vec,q_vec) and p_q_comm(p_vec[0:6],q_vec[0:6]):
            op3_with_phase = three_op_prod(vect_to_operator(p_str), vect_to_operator(q_str))
            p_lines.append([set([vect_to_operator(p_str), vect_to_operator(q_str), op3_with_phase[0]]), op3_with_phase[1]]) 
    [SCH_Skew.append(line) for line in p_lines if line not in SCH_Skew] # gets unique lines (i.e. ignores ordering)



# # Creating the Skew Complement

# In[61]:


Skew_comp = []
[Skew_comp.append(line) for line in W_5_2 if line not in SCH_Skew]
print(len(Skew_comp))


# # Building elliptic quadric Q^(-)(5,2)

# In[62]:


# Some necessary functions using the gamma matrix method

# function to give product of operators across a context in string form
def context_str_prod(context):
    return np.matmul(np.matmul(Ops[context[0]], Ops[context[1]]), Ops[context[2]])

# function to check sign of overall context, based on product of 3 operator matrices
def id_check(M):
    if (M.shape[0] == M.shape[1]) and (M == np.eye(M.shape[0])).all():
        return 1
    elif (M.shape[0] == M.shape[1]) and (M == -np.eye(M.shape[0])).all():
        return -1
    else:
        return 0


# In[63]:


# Creating the generators (Gamma_ij operators)

gamma_str = ['ZYI', 'YIX', 'XYI', 'IXY', 'YIZ', 'IZY', 'YYY']

# Creating the contexts of the doily
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
doily_contexts_str = []
for indices in doily_indices:
    context_str = []
    for pair in range(3): # for operator products we don't include the overall phase in the three_op_prod output
        context_str.append(three_op_prod(three_op_prod(gamma_str[indices[pair][0]-1], gamma_str[indices[pair][1]-1])[0],gamma_str[6])[0])
    doily_contexts_str.append(context_str)

# adds line parities to contexts and uses context sets to ignore ordering of operators
doily_contexts_with_signs = [[set(context), id_check(context_str_prod(context))] for context in doily_contexts_str]


# Creating the double-six contexts using gamma generators
six_one = [gamma for gamma in gamma_str[0:6]] # first set of double six
six_two = [three_op_prod(gamma, gamma_str[6])[0] for gamma in gamma_str[0:6]] # second set
six_two_ops = [np.matmul(Ops[gamma], Ops[gamma_str[6]]*(-1.j)) for gamma in gamma_str[0:6]] 

# Creates contexts as operator strings
double_six_contexts_str = []
for i in range(6):
    for j in range(6):
        if i != j:
            context_str = [six_one[i], six_two[j], three_op_prod(six_one[i], six_two[j])[0]]
            double_six_contexts_str.append(context_str)
            
# Adds line parities and un-orders contexts
double_six_contexts_str_with_signs = [[set(context), id_check(context_str_prod(context))] for context in double_six_contexts_str]

# Creates the elliptic quadric as the union of the doily and double-six
Elliptic_quadric = doily_contexts_with_signs + double_six_contexts_str_with_signs

# Identifies the negative lines in the elliptic quadric
neg_lines = []
for context in Elliptic_quadric:
    if context[1] == -1:
        neg_lines.append(context[0])
       


# # Creating the Complements

# In[64]:


Skew_comp = [] # complement of the skew-embedded split Cayley hexagon
[Skew_comp.append(line) for line in W_5_2 if line not in SCH_Skew]
print(len(Skew_comp))

Classical_comp = [] # complement of the classically-embedded split Cayley hexagon
[Classical_comp.append(line) for line in W_5_2 if line not in SCH_Classical]
print(len(Classical_comp))



# # Saving Arrangements and degrees

# In[84]:


Arrangements = [W_5_2, SCH_Classical, Classical_comp, SCH_Skew, Skew_comp, Elliptic_quadric]

# associated degree of contextuality for each arrangement
degrees = [63, 0, 0, 0, 24, 9]


# # Testing Contextuality on a Quantum Computer simulator

# In[66]:


# Checks operator, changes basis using relevant gates, makes measurement
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

# Checks operator, changes basis using relevant gates, delegates onto delegation qubit, changes basis back to computational
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
    


# # Creating a noise model

# In[80]:


IBMQ.save_account('<INSERT IBMQ API TOKEN>') # Use user credentials here
provider = IBMQ.get_provider(hub = 'hub', group = 'group', project = 'project') # insert instance parameters here

backend_choice = 'ibm_brisbane' # Use backend of choice from IBM Quantum Experience
quantum_comp=provider.get_backend(backend_choice) 
noise_model = NoiseModel.from_backend(quantum_comp)


# # Running over all lines in an arrangement on a quantum backend or simulator

# In[81]:


# change this to change the arrangement based on Arrangements array above.
Arrangement_choice = 4 # 4 for Skew Hexagon comeplement, 5 for elliptic quadric


Arrangement = Arrangements[Arrangement_choice] # Sets the arrangement based on choice

is_simulation = True # Set to "False" to run on real backend specified in quantum_comp
is_noisy = True # If running a simulator, can engage noise model based on backend

number_of_shots = 8192
simulator = Aer.get_backend('statevector_simulator')

qubits = [11,30,31] # set of register qubits used for context, check backend topology on IBM
d_qubits = [12,17,29] # set of register qubits used for delegations, check backend topology

contextuality_counts_real = [] # counts measured positive and negative lines, and the line parities for comparison

for line in Arrangement:
   context = line[0]
   parity = line[1]
   # creates circuit with necessary number of qubits based on register assignment, and 3 classical bits from delegation measurements
   circuit = QuantumCircuit(max(max(qubits), max(d_qubits))+1,3) 

   # Delegating each operator onto a new qubit
   op_num = 0
   for operator in context:
       for i in range(3):
           Op_delegate_measurement(operator[i], circuit, qubits[i], d_qubits[op_num])
       op_num += 1

   # Measuring on delegation qubits
   for j in range(3):
       circuit.measure(d_qubits[j],j)

   if is_simulation:
       local_simulator = Aer.get_backend('qasm_simulator')
       if is_noisy:
           result = execute(circuit, backend=local_simulator, noise_model = noise_model, shots=number_of_shots).result()
       else:
           result = execute(circuit, backend=local_simulator, shots=number_of_shots).result()
       backend_output = 'simulator'
   else:
       result = execute(circuit, backend=quantum_comp, shots=number_of_shots).result()
       backend_output = quantum_comp
   measures_dictionary = result.get_counts()

   pos_count = 0
   neg_count = 0
   for measure in measures_dictionary:
       if sum([int(bit) for bit in measure])%2 == 0:
           pos_count=pos_count+measures_dictionary[measure]
       else:
           neg_count=neg_count+measures_dictionary[measure]
   
   contextuality_counts_real.append([pos_count, neg_count, parity])

   


# In[83]:


# Checking Results

deg = degrees[Arrangement_choice] # degree of contextuality of the chosen arrangement

total_lines = len(contextuality_counts_real)

# Counts positive lines
parity_check = sum([1 for line in contextuality_counts_real if int(line[2]) == 1]) 

# Computes the metric \chi from the Cabello inequalities
contextuality_check_real = sum([int(line[2])*(line[0] - line[1])/number_of_shots for line in contextuality_counts_real])

print("Total # Lines: ",total_lines)
print("Total # Positive Lines: ",parity_check)
print("Total # Computationally Satisfied Lines: ",contextuality_check_real)
print("")

# Outputs the cabello limits depending on quantum mechanics or NCHV models
print("Cabello limit assuming QM: ", total_lines)
print("Cabello limit assuming NCHV: ", total_lines - 2*deg)

# Compares to computed result for \chi - if higher than NCHV limit then no NCHV models satisfy conditions
print("Result (", backend_output , "): ",contextuality_check_real)


# In[ ]:




