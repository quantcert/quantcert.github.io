import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_eval.run_circuit import *

H = matrix(field, [[1, 1], [1, -1]])/sqrt(2) # can only be used in symbolic ring
# H = matrix(field, [[1, 1], [1, -1]])/sqrt2 # fix for UniversalCyclotomicField
I2 = Matrix.identity(2)
Uf = Matrix.identity(4)


def U(f):
	""" Creates the oracle needed in the Deutsch algorithm
	args:
		f (dict): a dictionary mapping {0, 1} to itself
	returns:
		matrix: the matrix U such as U|x, y> = |x, f(x) XOR y> for all x, y in {0, 1}
	"""
	f_tot = {}
	for i in range(0, 4):
		# we operate bitwise, first bit is i%2, second is i//2,  
		# i is considered as the vector
		# Uf*(i1,i2) = (i1, i1 XOR i2)
		f_tot[i] = (i//2)*2 + (i%2)^^(f[i//2]) 
	result = Matrix(field, 4, 4)
	for i in f_tot:
		result[i,f_tot[i]]=1
	return result


def deutsch(f, output=False, output_file=False, file_path=None):
	""" Determines if a boolean function is constant or not
	args:
		f (dict): a dictionary mapping {0, 1} to itself
		output (bool): if true, outputs are enabled
		output_to_file (bool): if true, trace is returned in file, else it is printed
		file_path (str): path to the file to output the commands to
	returns:
		bool: False if f is constant, True otherwise
	"""
	if file_path:
		file = open(file_path, "w")
	else:
		if output and output_file:
			raise ValueError("output_file set to True but no valid path to file")
		file = None
	layers = [[H,H], [U(f)], [H,I2]]
	V0 = vector(field, [0,1,0,0])
	Vf = run(layers, V0, output, output_file, file)[0][-1]
	m = measure(Vf, 1)
	if m[1]:
		return True
	else:
		return False