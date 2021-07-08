import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_eval.run_circuit import *

X = matrix([[0, 1], [1, 0]])
H = matrix(field, [[1, 1], [1, -1]])/sqrt(2) # can only be used in symbolic ring
# H = matrix(field, [[1, 1], [1, -1]])/sqrt2 # fix for UniversalCyclotomicField
I2 = Matrix.identity(2)
Uf = Matrix.identity(4)


def U(f):
	""" Creates the oracle needed in the Deutsch-Jozsa algorithm
	args:
		f (dict): a dictionary mapping [0..(2**n)-1] to {0, 1}
	returns:
		matrix: the matrix U such as U|x, y> = |x, f(x) XOR y> for all x, y in 
			[0..(2**n)-1]*{0, 1}
	"""
	result = matrix(2*len(f))
	for i in range(0, len(f)):
		if f[i] == 0:
			result[2*i:2*i+2,2*i:2*i+2] = matrix.identity(2)
		else:
			result[2*i:2*i+2,2*i:2*i+2] = X
	return result


def deutsch_jozsa(f, output=False, output_file=False, file_path=None,epsilon=10^(-5)):
	""" Determines if a boolean function is constant or balanced
	args:
		f (dict): a dictionary mapping [0..(2**n)-1] to {0, 1}, f should either be 
			constant or balanced, in any other case, the result displayed won't carry any
			meaning, if f has a size other than a power of 2, an error will be raised
		output (bool): if true, outputs are enabled
		output_to_file (bool): if true, trace is returned in file, else it is printed
		file_path (str): path to the file to output the commands to
	returns:
		bool: True if f is constant, False if f is balanced
	"""
	if file_path:
		file = open(file_path, "w")
	else:
		if output and output_file:
			raise ValueError("output_file set to True but no valid path to file")
		file = None
	nb_qubits = log(len(f))/log(2)
	layers = [[H]*(nb_qubits+1), [U(f)], [H]*nb_qubits+[I2]]
	# v0 = [0]*nb_qubits+[1]
	# v0_kets = [ [1,0] if elt == 0 else [0,1] for elt in v0 ]
	# built_V0 = matrix([1])
	# for ket in v0_kets:
	# 	built_V0 = kronecker(built_V0,ket)
	V0 = vector(field, [0,1]+[0]*(2*len(f)-2))
	# print built_V0
	# print V0
	Vf = run(layers, V0, output, output_file, file)[0][-1]
	m = 1
	for i in range(1,nb_qubits+1):
		m *= measure(Vf, i)[0]
	if abs(1-m) < epsilon:
		return True
	elif abs(m) < epsilon:
		return False
	else: 
		return None 

n = 4
balanced_fs = []
for subset in Subsets(range(2^n),2^(n-1)):
	f = {}
	for i in range(2^n):
		f[i] = 1 if i in subset else 0
	balanced_fs += [f]

constant_fs = []
for const in [0,1]:
	f = {}
	for i in range(2^n):
		f[i] = const
	constant_fs += [f]

correct = True
for f in balanced_fs:
	correct = correct and not deutsch_jozsa(f)
for f in constant_fs:
	correct = correct and deutsch_jozsa(f)

print("Deutsch-Jozsa's algorithm is correct: " + str(correct))