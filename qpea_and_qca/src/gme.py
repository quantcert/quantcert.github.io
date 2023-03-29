import numpy as np
import random as rd


ket_0 = np.array([1, 0])
ket_1 = np.array([0, 1])


def first_coefficients_generation(n):
    """ Generates the very first coefficients for the calculation of GME by 
    quantum walk algorithm.

    Example :
        >>> first_coefficients_generation(3)
        [[-0.11275809+0.60216392j -0.47908013+0.62862266j]
         [-0.366532  -0.25643112j  0.67716151-0.58425138j]
         [ 0.74835654-0.53252389j -0.3953533 -0.00874982j]]

    :parameter int n: The number of qubits.
    :returns: np.array[list[complex]] -- A table of n lists containing 3 coefficients.
    """
    first_coefficient_list = []

    for counter in range(n):
        alpha_init = complex(rd.uniform(-1, 1), rd.uniform(-1, 1))
        beta_init = complex(rd.uniform(-1, 1), rd.uniform(-1, 1))
        norm = np.sqrt(abs(alpha_init) ** 2 + abs(beta_init) ** 2)
        first_coefficient_list.append([alpha_init / norm, beta_init / norm])

    first_coefficients = np.array(first_coefficient_list)
    return first_coefficients


def new_coefficients_generation(n, old_coefficients, step):
    """ Generates randomly new parameters for GME maximization by quantum walk algorithm.


    :parameter int n: The number of qubits
    :parameter np.array[list[complex]] old_coefficients: The table of the coefficients that didn't maximize the gme
    :parameter float step: The value of the step (used in the random walk algorithm)
    :returns: np.array[list[complex]] -- Table of new coefficients
    """
    new_coefficients_list = []
    for counter in range(n):
        alpha_init = old_coefficients[counter][0] + complex(rd.uniform(-1, 1) * step, rd.uniform(-1, 1) * step)
        beta_init = old_coefficients[counter][1] + complex(rd.uniform(-1, 1) * step, rd.uniform(-1, 1) * step)
        norm = np.sqrt(abs(alpha_init) ** 2 + abs(beta_init) ** 2)
        new_coefficients_list.append([alpha_init / norm, beta_init / norm])
    new_coefficients = np.array(new_coefficients_list)
    return new_coefficients


def vector_calculation(n, coeffs, first_base=ket_0, second_base=ket_1):
    """ Calculation of the vector on which the GME will be calculated.
    This vector depends on the basis which can be |0> or |00>. The default value is the first mentioned.

    :parameter int n: The number of qubits
    :parameter np.array[list[complex]] coeffs: The table of the coefficients
    :parameter list[int] first_base: The vector corresponding to |0> or |00>
    :parameter list[int] second_base: The vector corresponding to |1> or |11>
    :returns: np.array[list[complex]] -- Vector of length 2^n
    """
    vector = 1
    for i in range(n):
        vector_part = (coeffs[i][0] * first_base) + (coeffs[i][1] * second_base)
        vector = np.kron(vector, vector_part)
    return vector


def gme_calculation(vector, psi):
    """ Calculates "GME", the result of the calculation of the vector with the generated one, formed to be "separable".

    :parameter np.array[complex] vector: The separable, ie non entangled vector.
    :parameter np.array[complex] psi: The state vector from which the gme is computed.
    :returns: float -- The value of the calculation.

    """
    v_dagger = np.conj(vector).T
    gme = abs(v_dagger.dot(psi)) ** 2
    return gme


def gme_best_calculation(n, psi, step, alpha_min, c_max):
    """ Maximizes GME. The algorithm used here is called the Random walk method.
        We randomly generate first parameters which are used to calculate GME.
        The first value of GME is called first_gme. Then, we calculate new
        parameters based on the previous ones and a variable called the descent
        step. With the new parameters, we calculate a new value of GME. If this
        value is better than the previous one, we keep it and continue the
        researches until a counter is at its maximum value. The goal here is to
        take a big circle of research scope and to reduce it (by decreasing the
        decent step) more and more until the maximum value of GME is found.

    :parameter int n: The number of qubits.
    :parameter list[complex] psi: The vector for the calculation of GME.
    :parameter int step: The initial value of the descent step.
    :parameter float alpha_min: The minimum value of the descent step (which is the length of the radius).
    :parameter int c_max: The maximum value of the counter.
    :returns: float -- The value that is, expectedly, Mu, the optimal maximal value of the gme.
    """
    c = 0
    step_it = step
    alpha_min_it = alpha_min
    c_max_it = c_max

    x0_tab = first_coefficients_generation(n)
    first_vector = vector_calculation(n, x0_tab)
    first_gme = gme_calculation(first_vector, psi)

    # Initial parameters
    xbest_tab = x0_tab
    vector_best = first_vector
    gme_best = first_gme

    # Former iterations
    while step_it > alpha_min_it:
        while c < c_max_it:
            nouv_coef_mk = new_coefficients_generation(n, xbest_tab, step_it)
            new_vector = vector_calculation(n, nouv_coef_mk)
            gme_new_value = gme_calculation(new_vector, psi)
            if gme_new_value > gme_best:
                xbest_tab = nouv_coef_mk
                vector_best = new_vector
                gme_best = gme_new_value
                c = 1
            else:
                c += 1
        c = 1
        step_it = step_it / 2

    gme = 1 - gme_best
    return gme
