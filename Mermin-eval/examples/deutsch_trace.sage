#!/usr/local/bin/python
# coding: latin-1

# to use this file, run in the same folder `sage deutsch_trace.sage` It will 
# output a tex file that can be an input for another tex document such as 
# `sage_inclue.tex`. 

load("deutsch.sage")

f_dic = {
	0:1,
	1:0
}

deutsch(f_dic, True, True, "sage_output.sage.tex")
