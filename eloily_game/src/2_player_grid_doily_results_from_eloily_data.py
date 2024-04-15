#!/usr/bin/env python
# coding: utf-8


import numpy as np
import re
from ast import literal_eval
import csv


# # Extracting data from full Eloily results datafile



# Pull the full dataset of all Eloily 2-player games:
with open('eloily_full_2_player_game_results.csv', newline='') as eloily_data:
    csvread = csv.reader(eloily_data)
    batch_data = list(csvread)

alice_lines_str = [row[0] for row in batch_data] # array storing Alice's lines as strings
alice_lines_str = alice_lines_str[1:] # deletes header entry
alice_lines = [literal_eval(line) for line in alice_lines_str] # casts from string to array

bob_lines_str = [row[1] for row in batch_data] # array storing Bob's lines as strings
bob_lines_str = bob_lines_str[1:]
bob_lines = [literal_eval(line) for line in bob_lines_str] 

game_results_str = [row[4] for row in batch_data] # array storing game victory results in i'th position for game involving Alice's and Bob's i'th lines
game_results_str = game_results_str[1:]
game_results = [literal_eval(result) for result in game_results_str] 

print(alice_lines[0:3])
print(bob_lines[0:3])
print(game_results[0:3])
print(np.mean(game_results))
print(len(alice_lines))



# 
# # Extracting Grid results from Eloily data



# Functions to define spin operator multiplication for string representations
def operator_multiplication(op_1, op_2):
    op_set = [op_1, op_2]
    if op_1 == op_2:
        return 'I'
    elif op_set in [['X','Y'], ['Y','X']] or op_set in [['Z','I'], ['I','Z']] :
        return 'Z'
    elif op_set in [['X','Z'], ['Z','X']] or op_set in [['Y','I'], ['I','Y']]:
        return 'Y'
    else:
        return 'X'

def point_multiplication(point_1, point_2):
    return_str = ''
    for i in range(len(point_1)):
        return_str += operator_multiplication(point_1[i], point_2[i])
    return return_str



# Defining array of all intersections of lines
intersections = [[set(alice_lines[i]), set(bob_lines[i])] for i in range(len(alice_lines))]



# Construct an array of all grids in W(5,2), but not unique grids
sample_size = len(intersections)
all_grids = []
all_grids_points = []
for i in range(sample_size):
    intersection_1 = (intersections[i])
    for j in range(i+1,sample_size,1):
        intersection_2 = (intersections[j])
        if not any(line in intersection_1 for line in intersection_2):
            line_1_a, line_1_b = intersection_1[0], intersection_1[1]
            line_2_a, line_2_b = intersection_2[0], intersection_2[1]
            point_1 = line_1_a&line_1_b
            point_2 = line_2_a&line_2_b
            if len((line_1_a-point_1)&line_2_a)==1 and len(line_1_b&(line_2_b-point_2))==1:
                midpoint_a = list(line_1_a-point_1-(line_1_a&line_2_a))[0]
                midpoint_b = list(line_2_b-point_2-(line_1_b&line_2_b))[0]
                midpoint_c = list(line_2_a-point_2-(line_2_a&line_1_a))[0]
                midpoint_d = list(line_1_b-point_1-(line_2_b&line_1_b))[0]
                if point_multiplication(midpoint_a, midpoint_b) == point_multiplication(midpoint_c, midpoint_d):
                    newline_1 = set([midpoint_a, midpoint_b, point_multiplication(midpoint_a, midpoint_b)])
                    newline_2 = set([midpoint_c, midpoint_d, point_multiplication(midpoint_c, midpoint_d)])
                    grid = [line_1_a, line_1_b, line_2_a, line_2_b, newline_1, newline_2]
                    grid_points = []
                    all_grids.append(grid)
                    for line in grid:
                        for point in line:
                            grid_points.append(point)
                    grid_points = set(grid_points)
                    all_grids_points.append(grid_points)
    

# Checks if the grids are unique
index_arr = []
for grid_points in all_grids_points:
    index_arr.append(all_grids_points.index(grid_points))

# Makes array of unique grids, but with some zero entries
unique_grids_arr = list(np.zeros(max(index_arr)+1))
for index in index_arr:
    unique_grids_arr[index] = all_grids[index]
    
# Removes the zero entries
unique_nonzero_grids = []
for grid_val in unique_grids_arr:
    if grid_val != 0:
        unique_nonzero_grids.append(grid_val)
        
# Getting average result per grid
all_grid_results = []
all_grid_mean_results = []
for grid in unique_nonzero_grids:
    grid_results = []
    for line_1 in grid:
        for line_2 in grid:
            if line_1 != line_2:
                for intersection in intersections:
                    if line_1 in intersection and line_2 in intersection:
                        grid_results.append(game_results[intersections.index(intersection)])
    all_grid_results.append(set(grid_results))
    all_grid_mean_results.append(np.mean(grid_results))
    
# Finding the best grid
max_grid_result = max(all_grid_mean_results)
max_index = all_grid_mean_results.index(max_grid_result)
max_grid = unique_nonzero_grids[max_index]
print("Best grid: ",max_grid)
print("Best grid results: ",max_grid_result)


# # Extracting Doily results


# Function to convert single spin operators to vector components in PG(5,2)
def single_op_to_vector(op):
    if op == 'I':
        return '00'
    elif op == 'X':
        return '10'
    elif op == 'Z':
        return '01'
    else:
        return '11'

# Function to convert vector components in PG(5,2) to single spin operators
def vector_point_to_op(vect_point):
    if vect_point == '00':
        return 'I'
    elif vect_point == '10':
        return 'X'
    elif vect_point == '01':
        return 'Z'
    else:
        return 'Y'

# Converts full vector in PG(5,2) to 3-qubit spin operator string
def vector_to_operator(vect):
    vect_1, vect_2, vect_3 = vect[0:2], vect[2:4], vect[4:6]
    return vector_point_to_op(vect_1)+vector_point_to_op(vect_2)+vector_point_to_op(vect_3)

# Converts 3-qubit spin operator string to full vector in PG(5,2) 
def operator_to_vector(operator):
    vect = ''
    for op in operator:
        vect += single_op_to_vector(op)[0]
    for op in operator:
        vect += single_op_to_vector(op)[1]
    return vect

# Quadric in symplectic polar space W(5,2)
def Q_0(p):
    return (int(p[0])*int(p[3]) + int(p[1])*int(p[4]) + int(p[2])*int(p[5]))%2

# Inner product of two vectors in PG(5,2)
def inner_prod(p,q):
    return (int(p[0])*int(q[3]) + int(p[1])*int(q[4]) + int(p[2])*int(q[5]) + int(p[3])*int(q[0]) + int(p[4])*int(q[1]) + int(p[5])*int(q[2]))%2

# Symplectic form defining commutation relations between vectors
def Q_p(p,q):
    return (Q_0(p) + inner_prod(p,q))%2

# Provides all points that commute with operator given by vector p
def perp(p_op, ops):
    p = operator_to_vector(p_op)
    perp_set = []
    for op in ops:
        q = operator_to_vector(op)
        if not inner_prod(p,q):
            perp_set.append(op)
    return perp_set


# ## Making the set W(5,2)


# Gives all 3-qubit operators in string format
Ops = ['I', 'X', 'Y', 'Z']
Ops_str = []
for i in Ops:
    for j in Ops:
        for k in Ops:
            Ops_str.append(i+j+k)
Ops_str.remove('III') 

GQ_24_ops = [] # Set of points in GQ(2,4) (canonical labelling)
Non_elliptic_ops = [] # Set of points outside GQ(2,4)
for op in Ops_str:
    if op.count('I')==1:
        GQ_24_ops.append(op)
    else:
        Non_elliptic_ops.append(op)



doily_sets = [] # array storing all doilies in GQ(2,4)

# Doily can be defined as perp set of elements in W(5,2) outside of GQ(2,4)
for non_elliptic_op in Non_elliptic_ops:
    perp_set = perp(non_elliptic_op, GQ_24_ops)
    doily_sets.append(perp_set) 
doily_identity = perp('III',GQ_24_ops)


# ## Finding best doily results


doily_averages = []
doily_all_intersection_sets = []
for doily in doily_sets:
    doily_intersection_lines = []
    doily_intersections = []
    for intersection in intersections:
        flat_intersection = [op for op in list(intersection[0])] + [op for op in list(intersection[1])]
        if set(flat_intersection).issubset(set(doily)):
            doily_intersection_lines.append(intersection)
            doily_intersections.append(game_results[intersections.index(intersection)])
    doily_averages.append(np.mean(doily_intersections))
    doily_all_intersection_sets.append(doily_intersection_lines)
max_doily_result = max(doily_averages)
max_index = doily_averages.index(max_doily_result)
max_doily = doily_sets[max_index]
print("Best doily: ",max_doily)
print("Best doily results: ",max_doily_result)




