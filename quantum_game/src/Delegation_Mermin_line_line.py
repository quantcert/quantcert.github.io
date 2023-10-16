#!/usr/bin/env python
# coding: utf-8

# # This code provides the results for playing the Mermin line-line game using the delegation method

# In[4]:


import numpy as np
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel
import random
import itertools as it
from ast import literal_eval

orders = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
# Orders above correspond to orders on the doily data set


# ## Creating the Doily

# In[5]:


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

line1 = [Op1, Op3, Op2]
line2 = [Op3, Op5, Op4]
line3 = [Op5, Op7, Op6] # negative
line4 = [Op7, Op9, Op8]
line5 = [Op9, Op1, Op10]
line6 = [Op6, Op14, Op1]
line7 = [Op2, Op12, Op7]
line8 = [Op8, Op15, Op3]
line9 = [Op4,Op13,Op9]
line10 = [Op10, Op11, Op5]
line11 = [Op11, Op2, Op13] # negative
line12 = [Op12, Op4, Op14]
line13 = [Op13, Op6, Op15]
line14 = [Op14, Op8, Op11]
line15 = [Op15, Op10, Op12] # negative

# Storing the doily as an array of lines
Doily = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15]

neg_lines = [line3, line11, line15] # the negative lines


# ## Creating the Mermin Grids

# In[6]:


Mermin_grids_all = [np.array([['IX', 'XI', 'XX'],
                      ['ZI', 'IZ', 'ZZ'],
                      ['ZX', 'XZ', 'YY']]),
                   np.array([['XX', 'YZ', 'ZY'],
                      ['YY', 'YI', 'IY'],
                      ['ZZ', 'IZ', 'ZI']]),
                   np.array([['ZY', 'YZ', 'XX'],
                      ['XZ', 'ZX', 'YY'],
                      ['YX', 'XY', 'ZZ']]),
                    np.array([['XI', 'XZ', 'IZ'],
                      ['XY', 'YX', 'ZZ'],
                      ['IY', 'ZY', 'ZI']]),
                   np.array([['XX', 'YZ', 'ZY'],
                      ['IX', 'ZX', 'ZI'],
                      ['XI', 'XY', 'IY']]),
                   np.array([['IY', 'YI', 'YY'],
                      ['XY', 'YX', 'ZZ'],
                      ['XI', 'IX', 'XX']]),
                   np.array([['ZX', 'IX', 'ZI'],
                      ['XY', 'YX', 'ZZ'],
                      ['YZ', 'YI', 'IZ']]),
                   np.array([['YX', 'IX', 'YI'],
                      ['ZY', 'ZI', 'IY'],
                      ['XZ', 'ZX', 'YY']]),
                   np.array([['YZ', 'IZ', 'YI'],
                      ['XY', 'XI', 'IY'],
                      ['ZX', 'XZ', 'YY']]),
                   np.array([['XZ', 'IZ', 'XI'],
                      ['ZY', 'YZ', 'XX'],
                      ['YX', 'YI', 'IX']])]

neg_lines_unordered = [['XX', 'ZZ', 'YY'], ['YZ', 'XY', 'ZX'], ['XZ', 'ZY', 'YX']]
neg_lines_perms = [list(it.permutations(line, 3)) for line in neg_lines_unordered]
neg_lines = np.reshape(neg_lines_perms, (18,3))
neg_set = [set(line) for line in neg_lines_unordered] # storing negative lines as unordered sets


# In[7]:


# Formatting the grids to show unordered rows and columns, for comparison with results
Mermin_grids_rows_cols = []
for Grid in Mermin_grids_all:
    grid_rows = [set(row) for row in Grid]
    grid_cols = [set(col) for col in Grid.T]
    Mermin_grids_rows_cols.append([grid_rows,grid_cols])


# In[8]:


# Pulling the data from the Doily line-line game
# Each line is the ordering, and the results for each context of the doily

doily_line_line_data = []
with open("Doily_line_line_results_with_lines.txt",'r') as file:
    for line in file:
        order = literal_eval(line[0:9])
        Alice_line = literal_eval(line[9:27])
        Bob_line = literal_eval(line[27:45])
        result = literal_eval(line[45:])
        doily_line_line_data.append([order,Alice_line,Bob_line,result])



# # Pulling results out for Mermin Grids from Doily data

# In[14]:


All_Mermin_averages = []
for Grid in Mermin_grids_rows_cols: # Looping over all mermin grids
    grid_results = []
    for order in orders: # Looping over all orders
        order_results = []
        order_data = [line for line in doily_line_line_data if line[0] == order] # limiting data to that of that order
        for result in order_data: # Checking each individual ordering
            Alice = set(result[1])
            Bob = set(result[2])
            if (Alice in Grid[0] and Bob in Grid[1]) or (Alice in Grid[1] and Bob in Grid[0]): # checking if the lines exist in the grid
                order_results.append(result[3])
        grid_results.append(np.mean(order_results)) # appending the best order for that grid
    All_Mermin_averages.append(max(grid_results))
print(All_Mermin_averages)


# In[15]:


# Identifying out the best grid
best_grid_avg = max(All_Mermin_averages)
best_index = All_Mermin_averages.index(best_grid_avg)
print('Best results are for Grid ',best_index)

print("Grid:")
print(Mermin_grids_all[best_index])
print("Average: ",best_grid_avg)


# In[ ]:




