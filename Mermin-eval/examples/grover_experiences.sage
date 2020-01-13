import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_eval.grover import *

for n in range(4,5):
  v = target_state_ket_string_to_vector("0"*n)
  grover(v, verbose=True, file_name="grover-results/"+str(n)+"qbit.csv", use_precomputed=True, artificial=True)