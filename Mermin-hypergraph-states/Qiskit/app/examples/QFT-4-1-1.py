import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_on_qiskit import QFT

QFT._main(1,1,4,"QFT-optimized-coefficients/1-1-4.csv")