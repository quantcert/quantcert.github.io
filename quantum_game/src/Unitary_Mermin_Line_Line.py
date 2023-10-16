#!/usr/bin/env python
# coding: utf-8

# # The following code covers the Mermin Line-Line game using the Unitary method

# In[1]:


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



# In[2]:


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



# # Creating the 10 distinct Mermin Grids (up to line permutations)

# In[3]:


# Note here we include some sign changes to ensure all columns are negative contexts and all rows are positive
Mermin_grids_all = [np.array([[IX, XI, XX],
                      [ZI, IZ, ZZ],
                      [-ZX, -XZ, YY]]),
                   np.array([[XX, -YZ, -ZY],
                      [YY, YI, IY],
                      [ZZ, IZ, ZI]]),
                   np.array([[ZY, YZ, XX],
                      [XZ, ZX, YY],
                      [YX, XY, ZZ]]),
                    np.array([[XI, XZ, IZ],
                      [XY, YX, ZZ],
                      [-IY, ZY, -ZI]]),
                   np.array([[XX, YZ, ZY],
                      [IX, ZX, ZI],
                      [-XI, XY, -IY]]),
                   np.array([[IY, YI, YY],
                      [XY, YX, ZZ],
                      [-XI, -IX, XX]]),
                   np.array([[ZX, IX, ZI],
                      [XY, YX, ZZ],
                      [YZ, -YI, -IZ]]),
                   np.array([[YX, IX, YI],
                      [ZY, ZI, IY],
                      [XZ, -ZX, -YY]]),
                   np.array([[YZ, IZ, YI],
                      [XY, XI, IY],
                      [ZX, -XZ, -YY]]),
                   np.array([[XZ, IZ, XI],
                      [ZY, YZ, XX],
                      [YX, -YI, -IX]])]


# In[4]:


# Running the game on IBM quantum experience


# In[10]:


is_simulation = True # Determines whether it's run on a simulator or real quantum computer
include_noise = False # Determines whether to include noise model

# Loading the IBM quantum experience account - include own credentials here
IBMQ.load_account()

# Determines project - include own info here
provider = IBMQ.get_provider(hub = 'hub', group = 'group', project = 'project')
quantum_comp=provider.get_backend('ibm_lagos')
noise_model = NoiseModel.from_backend(quantum_comp)

mermin_successes = []

# Cycling through each of the 10 grids
for Mermin_grid_init in Mermin_grids_all:

    # Creating the arrays for the reordered grids, permutations defining the reorderings, and the success rates of each grid
    M_rearranged = []
    grid_success = []

    # Creating the permutations that will reorder each grid
    comb_list = list(it.permutations([0,1,2],3))
    for switch_col in comb_list:
        for switch_row in comb_list:
            Mermin_grid_reord = Mermin_grid_init[:,switch_col]
            Mermin_grid_reord = Mermin_grid_reord[[switch_row]]
            M_rearranged.append(Mermin_grid_reord)
            switch_perm.append([switch_row,switch_col])
    
    arrangement = 0
    
    # Running the code on IBMQ for each reordered grid
    for Mermin_grid in M_rearranged:
        
        # Creating arrays of the unitary matrices for each rown and column for Alice and Bob
        row_transforms = []
        col_transforms = []

        for i in range(3):
            # takes the first two operators in each context and makes a unitary transformation for that context
            row_transform = Unitary_Transform(Mermin_grid[0,i,0], Mermin_grid[0,i,1]).evalf(11)
            col_transform = Unitary_Transform(Mermin_grid[0,0,i], Mermin_grid[0,1,i]).evalf(11)
            row_transforms.append(row_transform)
            col_transforms.append(col_transform)

        Success_matrix = [] # Stores the success rates for a given grid arrangement

        # Running through each column and row to apply the unitary transforms in a circuit
        for col in range(3):
            S_vect = [] # Success rate for a given column
            for row in range(3):

                # Calling unitary transformation matrices for the row / column to 8 significant figures
                Unitary_A = row_transforms[row].evalf(8)
                Unitary_B = col_transforms[col].evalf(8)

                # Creating circuit
                q = QuantumRegister(4)
                c = ClassicalRegister(4)
                circuit = QuantumCircuit(q,c)
                
                # Creating the shared entangled state Psi
                circuit.h (q[0]);
                circuit.h (q[1]);
                circuit.cx (q[0],q[2]);
                circuit.cx (q[1],q[3]);
                
                # Implementing the unitary transformations
                circuit.unitary(Unitary_A,[0,1])
                circuit.unitary(Unitary_B,[2,3])
                
                # Measuring onto the classical register
                circuit.measure (q[0],c[0]);
                circuit.measure (q[1],c[1]);
                circuit.measure (q[2],c[2]);
                circuit.measure (q[3],c[3]);
                
                # Executing job
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
                    machine = "IBM_Lagos"
                    result = execute(circuit, backend=quantum_comp, shots=8192).result()

                # Getting results
                qubit_results = [key for key in result.get_counts(circuit)]
                qubit_counts = list(result.get_counts(circuit).values())
                total_counts = sum(qubit_counts)

                # Computing Alice and Bob's answers
                Alice_results = [qubits[2:4] for qubits in qubit_results]
                Bob_results = [qubits[0:2] for qubits in qubit_results]
                
                # Alice's answer before performing the skew-flip protocol
                Alice_ans_pre_skew_flip = [qubits+str((int(qubits[0])+int(qubits[1]))%2) for qubits in Alice_results]

                # Performing the skew-flip protocol
                Alice_y = [int(not np.imag(Mermin_grid[0,row,i]).any()==np.zeros((4,4)).any()) for i in range(3)]
                Alice_ans_y = []
                # Flips Alice's bit for each Y
                for result in Alice_ans_pre_skew_flip:
                    new_result = str()
                    for i in range(3):
                        new_result_digit = str((int(result[i])+Alice_y[i])%2)
                        new_result += new_result_digit
                    Alice_ans_y.append(new_result)
                    
                # Bob's answer including third computed bit from line parity
                Bob_ans = [qubits+str((int(qubits[0])+int(qubits[1])+1)%2) for qubits in Bob_results]

                # Check to see if players have won the game
                checklist = [Alice_ans_y[i][col] == Bob_ans[i][row] for i in range(len(Alice_ans_y))]

                # Computes the success percentage across each shot
                true_counts = [checklist[i]*qubit_counts[i] for i in range(len(qubit_results))]
                success_rate = sum(true_counts)/float(total_counts)
                S_vect.append(round(success_rate, 5))
                
            # Success matrix - gives % successful runs at each intersection of the Mermin grid
            Success_matrix.append(S_vect)
            
        # Gives an array of all success matrices
        grid_success.append(np.mean(np.array(Success_matrix).T))
 
    mermin_successes.append(max(grid_success)) # taking the best result out of all orderings for a given grid


# Pulling the best grid and ordering as the success matrix
global_victory = max(mermin_successes)

# Printing global result
print("Global victory for ", machine, global_victory)


# In[ ]:




