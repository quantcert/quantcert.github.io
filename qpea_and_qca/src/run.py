from qpea import *

# For test running
test = True

if test:
    nb_qubits = 3
    case = "Z"
    n_second_register = 1
    run(nb_qubits, n_second_register, case)
else:
    nb_qubits = 3
    n_second_register = 1
    while nb_qubits < 6:
        for case in ["Z", "G"]:
            run(nb_qubits, n_second_register, case)
        clear = True  # To separate the different figures
        nb_qubits += 1
