#!/usr/bin/env python
# coding: utf-8

# # This code provides the results for playing the Mermin point-line game using the delegation method

# In[2]:


import numpy as np
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel
import random
import itertools as it

orders = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
# Orders above correspond to orders on the doily data set


# ## Creating the Doily

# In[27]:


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

# In[29]:


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


# In[22]:


# Pulling the data from the Doily point-line game
# Each line is the ordering, and the results for each context of the doily

doily_point_line_data = []
with open("Doily_point_line_results.txt",'r') as file:
    for line in file:
        data = line.split()
        results =[data[0][0:3],float(data[0][4:-2])]
        for datum in data[1:]:
            results.append(float(datum[:-2]))
        doily_point_line_data.append(results)

# Forming a dictionary based on the orders and the results per point of the doily
dict_pl = {}
for i in range(6):
    dict_pl[doily_point_line_data[i][0]] = doily_point_line_data[i][1:]


# In[23]:


keys = list(dict_pl.keys()) # orders of the doily

# Forming an array of all results over all orders
Results_array = [dict_pl[key] for key in keys]

# Getting best result per context over the orders
doily_max_results = [max(line) for line in np.transpose(Results_array)]

# Attaching relevant line to each result
doily_max_results_w_lines = []
for i in range(len(Doily)):
    doily_max_results_w_lines.append([set(Doily[i]), doily_max_results[i]])


# In[31]:


# Formatting the grids to show unordered rows and columns, for comparison with results
Mermin_grids_rows_cols = []
for Grid in Mermin_grids_all:
    grid_rows = [set(row) for row in Grid]
    grid_cols = [set(col) for col in Grid.T]
    Mermin_grids_rows_cols.append([grid_rows,grid_cols])


# ## Pulling results out for Mermin Grids from Doily data

# In[33]:


All_Mermin_averages = [] # Storing the average result for each Mermin grid
for Grid in Mermin_grids_rows_cols: # Grids of unordered lines
    grid_results = []
    for result in doily_max_results_w_lines: # Using the doily results array with best result per context
        Bob_line = result[0] # Setting Bob's line
        if (Bob_line in Grid[0]) or (Bob_line in Grid[1]):
            grid_results.append(result[1])
    grid_average = np.mean(grid_results)
    All_Mermin_averages.append(grid_average)


# In[34]:


# Identifying the best grid

best_grid_avg = max(All_Mermin_averages)
best_index = All_Mermin_averages.index(best_grid_avg)
print('Best results are for Grid ',best_index)
Grid = Mermin_grids_all[best_index]
print("Grid:")
print(Grid)
print("Average: ",best_grid_avg)


# In[ ]:




