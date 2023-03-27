import os
import sys
from pathlib import Path
import numpy as np

sys.path.append(os.getcwd() + "/../")

from mermin_on_qiskit import coefficients_shapes as cs
from mermin_on_qiskit import evaluation
from mermin_on_qiskit import hypergraphstates as hg
from mermin_on_qiskit.hypergraphstates_optimization import hypergraphstates as hgo
from mermin_on_qiskit.hypergraphstates_optimization import mermin_polynomials as mpo

def create_folder_if_needed(_filename):
    dirname = os.path.dirname(_filename)
    if not os.path.exists(dirname) and dirname != '':
        os.makedirs(dirname)


# ### var initializations
# random.seed(1)
n = 3
hyperedges = [[0, 1, 2]]
# temportary fix
run_id = 1  # Used to register run data

chemin_fichier = Path("__logs__/")
graph_states_path = Path("hypergraphstate/")
mu_path = Path("mu/")
probability_path = Path("probability/")
filename = str(n) + "-" + str(run_id) + ".txt"

if __name__ == "__main__":
    # ### Calcul du graphstate
    state_vector_init = hgo.state_vector_initialisation(n)
    state_vector = hgo.hyperedges_computation(n, state_vector_init, hyperedges)

    graph_filename = chemin_fichier / graph_states_path / filename
    create_folder_if_needed(graph_filename)

    with open(graph_filename, "w") as mon_vecteur:
        mon_vecteur.write("The state vector is :")
        mon_vecteur.write("\n" + str(state_vector))

    # ### mu optimization
    a_a_p_coeffs = mpo.xbest_calculation(
        n, False, vector=state_vector, file_path=graph_filename, alpha=5, 
        alpha_minimum=0.00001, c_maximum=40, saving_file=True)
    a_a_p_coeffs_packed = cs.coefficients_format_mixed_to_packed(a_a_p_coeffs)
    # at this point, the operators are stored in reverse order, this is fixed by 
    # the following loop. TODO : find a better fix ?
    a_a_p_coeffs_packed = [operator[::-1] for operator in a_a_p_coeffs_packed]

    # ### circuit preparation
    circuit = hg.circuit_initialisation(n, hyperedges)
    hg.edges_layout(n, hyperedges, circuit)

    # ### Mermin evaluation
    mermin_eval = evaluation.evaluate_polynomial(
        n, circuit, a_a_p_coeffs_packed, shots=1024, is_simulation=True)

    print(mermin_eval)

    prob_filename = chemin_fichier / probability_path / filename
    create_folder_if_needed(prob_filename)

    with open(prob_filename, "w") as vector_probability:
        vector_probability.write("\n" + str(a_a_p_coeffs))
        vector_probability.write("\n")
        vector_probability.write("The probability is :")
        vector_probability.write("\n")
        vector_probability.write("\n" + str(mermin_eval))
