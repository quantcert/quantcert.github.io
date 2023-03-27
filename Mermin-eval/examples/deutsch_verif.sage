#!/usr/local/bin/python
# coding: latin-1

# to use this file, run in the same folder `sage deutsch_verif.sage`. It will 
# output in the terminal either Deutsch's algorithm is correct or not.

load("deutsch.sage")

const_f_dic = [
	{0:0,1:0},
	{0:1,1:1}
]
var_f_dic = [
	{0:0,1:1},
	{0:1,1:0}
]
is_correct = True
for f in var_f_dic:
	is_var = deutsch(f)
	if not is_var:
		is_correct = False
for f in const_f_dic:
	is_var = deutsch(f)
	if is_var:
		is_correct = False

print("Deutsch's algorithm is correct: " + str(is_correct))
