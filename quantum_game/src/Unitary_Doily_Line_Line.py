#!/usr/bin/env python
# coding: utf-8

# # The following code covers the Doily Line-Line game using the Unitary method

# In[9]:


import numpy as np
import sympy as sp
from itertools import combinations
import itertools as it
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit import QuantumCircuit, transpile
from qiskit.tools.monitor import job_monitor
from qiskit.providers.aer.noise import NoiseModel

# Some functions to create the unitary basis transformation
# Normalises a "unitary" matrix of the form A*A^+ = cI for some scalar c.
def normalise(mat):
    matrix = sp.Matrix(mat)
    norm = sp.sqrt(np.matmul(np.conj(matrix).T,matrix)[0,0])
    if norm == 0: 
       return matrix
    return matrix / norm



# This function gives the span of the matrix, i.e. two 4-dim vectors, associated to the given eigenvalue
def eigenvect_span(mat, eig_val_sample):
    
    # The eig_val_sample input is the value of the eigenvalue for which we're getting the span, i.e. +1 or -1
    eigvals = np.linalg.eig(mat)[0]
    eigvals = [round(float(eig), 1) for eig in eigvals] # need this line as complex spans were giving floating point errors
    eigvects = np.linalg.eig(mat)[1]
    
    # below gives a 4x4 matrix where only the selected span vectors are nonzero
    all_vectors = np.where(np.array(eigvals) == eig_val_sample, np.array(eigvects), 0)
    
    # removes zero vectors to just give the span
    eig_span = all_vectors[:, ~np.all(all_vectors == 0, axis = 0)]
    
    return eig_span.T # Outputs a 4x2 matrix of two vectors


    
# This function finds the common vector between two eigenvector spans of the two operators in question
def common_span_vector(arr_1, arr_2, eig_val_sample_1, eig_val_sample_2):
    
    # input the +1 or -1 eigenvalue for each of matrix 1 and 2. Computes the span for each matrix for its chosen eigenvalue
    vector_1 = eigenvect_span(arr_1, eig_val_sample_1)
    vector_2 = eigenvect_span(arr_2, eig_val_sample_2)
    
    # concatenates the two vector spaces (one negated) to find the nullspace - the one-vector space of common eigenvectors
    concat = np.concatenate((vector_1, -vector_2))
    concat = sp.Matrix(concat).T
    concat = concat.evalf(10) # removes floating point errors via rounding
    common_vec = concat.nullspace() # should return a single vector from matrix 1 - the common vector to both spans
    a = common_vec[0][0] # now have values for the unknowns a and b that are the coefficients in the span of the first matrix
    b = common_vec[0][1]
    w1 = vector_1[0]
    w2 = vector_1[1]
    if a != 0:# create the vector that's in the span of both matrices
        return (w1*a + w2*b)/a 
    elif b != 0:
        return (w1*a + w2*b)/b
    else:
        return (w1*a + w2*b)



# This function outputs the unitary transformation matrix
def Unitary_Transform(arr_1, arr_2):
    # eigenvalue combinations, given in a chosen order. Left one from arr_1 and right from arr_2
    eigenvalue_combos = [[1,1],[1,-1],[-1,1],[-1,-1]]
    
    A = []
    
    # the below goes through all combinations of eigenvalue choices for the spans
    for [i,j] in eigenvalue_combos:
        vector = sp.Matrix(common_span_vector(arr_1, arr_2, i, j)) # gets the vector common to the eigen spans of each array and each eigenvalue combination
        A.append(vector)
    
    # conjugates the matrix A for the reverse unitary operation, and normalises it. Note the appending above already gives a transposed matrix of vectors
    return normalise(np.reshape(np.conj(A), (4,4)))

# ---------------------------------------
# Some functions to determing intersecting lines in the Doily geometry
# Finding all lines for a given intersection point and all transform matrices for those lines
def intersecting_lines(Op, perm):
    Lines = [line for line in Doily_lines if Op in line]
    Reordered_lines = []
    for line in Lines:
        Reordered_lines.append([line[j] for j in perm])
    comb_list = list(it.combinations(Reordered_lines,2))
    return comb_list


# Below gives the unitary transformation for each of the combination of lines going through a given point
def point_transformations(Op, player, perm):
    intersecting_pairs_of_lines = intersecting_lines(Op, perm)
    num_combs = len(intersecting_pairs_of_lines)

    alice_transforms = []
    bob_transforms = []
    # Below just selects first two operators for each of A and B's lines
    for i in range(num_combs):
        alice_transform = Unitary_Transform(Ops[intersecting_pairs_of_lines[i][0][0]], Ops[intersecting_pairs_of_lines[i][0][1]]).evalf(11)
        bob_transform = Unitary_Transform(Ops[intersecting_pairs_of_lines[i][1][0]], Ops[intersecting_pairs_of_lines[i][1][1]]).evalf(11)
        alice_transforms.append(alice_transform)
        bob_transforms.append(bob_transform)
    if player[0] == 'A':
        return alice_transforms
    elif player[0] == 'B':
        return bob_transforms


# In[4]:


# Defining the various operators

I = [[1, 0],
    [0, 1]]

X = [[0, 1],
    [1, 0]]

Y = [[0, -1.j],
    [1.j, 0]]

Z = [[1.0, 0],
    [0, -1.0]]

IX = np.kron(I, X)
XI = np.kron(X, I)
XX = np.kron(X, X)

IY = np.kron(I, Y)
YI = np.kron(Y, I)
YY = np.kron(Y, Y)

IZ = np.kron(I, Z)
ZI = np.kron(Z, I)
ZZ = np.kron(Z, Z)

XY = np.kron(X, Y)
YX = np.kron(Y, X)

XZ = np.kron(X, Z)
ZX = np.kron(Z, X)

YZ = np.kron(Y, Z)
ZY = np.kron(Z, Y)



# # Making a dictionary to go from text to operators

# In[5]:


Ops_arr = [I,X,Y,Z]
Double_ops = []
for i in Ops_arr:
    for j in Ops_arr:
        double_op = np.kron(i,j)
        Double_ops.append(double_op)
Double_ops.remove(Double_ops[0])
# print(Double_ops_w_neg)

Ops_str = ["I","X","Y","Z"]
Double_ops_str = []
for i in Ops_str:
    for j in Ops_str:
        double_op = i+j
        Double_ops_str.append(double_op)
Double_ops_str.remove('II')
# print(Double_op_str_w_neg)

Ops = {}
for i in range(len(Double_ops)):
    Ops[Double_ops_str[i]] = Double_ops[i]


# # Creating the Doily

# In[6]:


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
Doily_lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15]

# The negative lines of this labelling
negative_lines = [line3, line11, line15]


# In[ ]:





# In[ ]:





# In[11]:


# Running the results on a quantum circuit

is_simulation = True # Determines whether it's run on a simulator or real quantum computer
include_noise = False # Determines whether to include noise model

# Loading the IBM quantum experience account - include own credentials here
IBMQ.load_account()

# Determines project - include own info here
provider = IBMQ.get_provider(hub = 'hub', group = 'group', project = 'project')
quantum_comp=provider.get_backend('ibm_lagos')
noise_model = NoiseModel.from_backend(quantum_comp)

# Setting permutations for orderings
perms = list(it.permutations([0,1,2], 3))

doily_successes = [] # array storing success rates for each permutation
for permutation in perms:
    Success_matrix = [] # Stores success rates for a given permutation
    for point in Doily_points:
        S_vect_point = [] # stores success rate for a given set of intersecting lines
        alice_transforms = point_transformations(point, 'Alice', permutation)
        bob_transforms = point_transformations(point, 'Bob', permutation)
        point_intersecting_lines = intersecting_lines(point, permutation)
        for n in range(len(point_intersecting_lines)):
             # Reordering all lines given a permutation
            Alice_line = point_intersecting_lines[n][0]
            Bob_line = point_intersecting_lines[n][1]

            col = Alice_line.index(point) # "column" index
            row = Bob_line.index(point) # "row" index

            # Creating unitary transformation matrices for each pair of intersecting lines at the point
            Unitary_A = alice_transforms[n].evalf(8)
            Unitary_B = bob_transforms[n].evalf(8)

            # Alice Y value if she has any Y's in her line
            Alice_y = [str(Alice_line[k]).count('Y')%2 for k in range(3)]

            # Additional S Factor (line parity) if either's line is negative, to determine third measurement
            if point_intersecting_lines[n][0] in negative_lines:
                s_factor_a = 1
            else:
                s_factor_a = 0
            if point_intersecting_lines[n][1] in negative_lines:
                s_factor_b = 1
            else:
                s_factor_b = 0

            # Creating circuit
            q = QuantumRegister(4)
            c = ClassicalRegister(4)
            circuit = QuantumCircuit(q,c)
            
            # Creating the shared state Psi
            circuit.h (q[0]);
            circuit.h (q[2]);
            circuit.cx (q[0],q[1]);
            circuit.cx (q[2],q[3]);
            
            # Setting the unitary transformations
            circuit.unitary(Unitary_A,[0,2])
            circuit.unitary(Unitary_B,[1,3])
            
            # Measuring onto output bits
            circuit.measure (q[0],c[0]);
            circuit.measure (q[2],c[1]);
            circuit.measure (q[1],c[2]);
            circuit.measure (q[3],c[3]);

            # Getting results
            if is_simulation:
                simulator = Aer.get_backend('qasm_simulator')
                if include_noise:
                    machine = "noisy simulator: "
                    result = execute(circuit,backend = simulator, 
                                     noise_model=noise_model).result()
                else:
                    machine = "noiseless simulator: "
                    result = execute(circuit,backend = simulator).result()
            else:
                machine = quantum_comp
                result = execute(circuit, backend=quantum_comp, shots=8192).result()
                
            qubit_results = [key for key in result.get_counts(circuit)]
            qubit_counts = list(result.get_counts(circuit).values())
            total_counts = sum(qubit_counts)

            Alice_results = [qubits[2:4] for qubits in qubit_results]
            Bob_results = [qubits[0:2] for qubits in qubit_results]
            
            # Alice's answers before skew-flip protocol, including third bit
            Alice_ans_pre_flip = [qubits+str((int(qubits[0])+int(qubits[1])+s_factor_a)%2) for qubits in Alice_results]
            Alice_ans_y = [] # Storing Alice's answers post skew-flip protocol
            for result in Alice_ans_pre_flip:
                new_result = str()
                for i in range(3):
                    new_result_digit = str((int(result[i])+Alice_y[i])%2)
                    new_result += new_result_digit
                Alice_ans_y.append(new_result)
                
            # Bob's answers, including third bit
            Bob_ans = [qubits+str((int(qubits[0])+int(qubits[1])+s_factor_b)%2) for qubits in Bob_results]
            
            # Determines whether they won or lost the game
            checklist = [Alice_ans_y[i][col] == Bob_ans[i][row] for i in range(len(Alice_ans_y))]

            # Determines how many successful plays
            true_counts = [checklist[i]*qubit_counts[i] for i in range(len(qubit_results))]
            success_rate = sum(true_counts)/float(total_counts)
            S_vect_point.append(round(success_rate, 1))

        Success_matrix.append(S_vect_point)

    doily_successes.append(np.mean(Success_matrix))


global_victory = max(doily_successes) # Using best permutation result

# Printing global result
print("Global victory for ", machine, global_victory)


# In[ ]:




