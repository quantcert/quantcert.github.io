import numpy as np
import random as rd


# The Pauli matrices
matrix_x = np.array([[0, 1], [1, 0]])
matrix_y = np.array([[0, complex(0, -1)], [complex(0, 1), 0]])
matrix_z = np.array([[1, 0], [0, -1]])


def first_coefficients_generation(n):
    """ Generates the very first coefficients for the calculation of MU by random walk algorithm.

    Example :
        >>> first_coefficients_generation(3)
        [[-0.71636549  0.16623039 -0.67763407]
         [-0.36958952  0.54390424  0.75337359]
         [ 0.48738401 -0.78929548  0.37345613]
         [-0.51019129 -0.53759251  0.6713413 ]
         [ 0.89590429 -0.08368889  0.43629311]
         [ 0.42877383  0.76027758  0.48798669]]

    :parameter int n: The number of qubits.
    :return: np.array(list(complex)) -- A table of 2 * n lists containing three (03) coefficients each.

    """
    first_coefficients_list = []
    maximum = 2 * n
    for counter in range(maximum):
        alpha_init = rd.uniform(-1, 1)
        beta_init = rd.uniform(-1, 1)
        gamma_init = rd.uniform(-1, 1)
        norm = np.sqrt(alpha_init ** 2 + beta_init ** 2 + gamma_init ** 2)
        first_coefficients_list.append([alpha_init / norm, beta_init / norm, gamma_init / norm])
    first_coefficients_array = np.array(first_coefficients_list)
    return first_coefficients_array


def new_coefficients_generation(n, old_array, step):
    """ Generates randomly new parameters for MU maximization by random walk algorithm.


    :parameter int n: The number of qubits
    :parameter np.array(list(float)) old_array: The table of the coefficients that didn't maximize Mu
    :parameter float step: The value of the step (used in the random walk algorithm)
    :return: np.array(list(complex)) -- Table of new coefficients

    """
    maximum = 2 * n
    new_coefficients_list = []
    for counter in range(maximum):
        alpha_init = old_array[counter][0] + (rd.uniform(-1, 1) * step)
        beta_init = old_array[counter][1] + (rd.uniform(-1, 1) * step)
        gamma_init = old_array[counter][2] + (rd.uniform(-1, 1) * step)
        norm = np.sqrt(alpha_init ** 2 + beta_init ** 2 + gamma_init ** 2)
        new_coefficients_list.append([alpha_init / norm, beta_init / norm, gamma_init / norm])
    new_coefficients_array = np.array(new_coefficients_list)
    return new_coefficients_array


def a_matrix(n, t):
    """ Computes the matrix "a" of the Mermin polynomial based on three Pauli's matrices.

    Example :
        >>> a_matrix(3, [[-0.7,0.2,-0.7],[-0.36,0.5,0.7],[0.5,-0.8,0.4],[-0.5,-0.5,0.7], [0.9,-0.1,0.4],[0.4,0.7,0.5]])
        [[ 0.43+0.j    0.9 +0.08j]
         [ 0.9 -0.08j -0.43+0.j  ]]

    :parameter int n: The number of qubits.
    :parameter: list(int) t : A list of 2 * n lists containing three (03) floats.
    :return: np.array(complex) -- The matrix a of dimension 2*2.

    """
    a = t[2 * n - 2][0] * matrix_x + t[2 * n - 2][1] * matrix_y + t[2 * n - 2][2] * matrix_z
    return a


def a_prime_matrix(n, t):
    """ Computes the matrix "a'" of the Mermin polynomial ased on three Pauli's matrices.

    Example :
    >>> a_prime_matrix(3, [[-0.7,0.2,-0.7],[-0.4,0.5,0.7],[0.5,-0.8,0.4],[-0.5,-0.5,0.7],[0.9,-0.1,0.4],[0.4,0.7,0.5]])
        [[ 0.48+0.j    0.43-0.76j]
         [ 0.43+0.76j -0.48+0.j  ]]

    :parameter int n: The number of qubits.
    :parameter: list(int) t : A list of 2 * n lists containing three (03) floats.
    :return: np.array(complex) -- The matrix a' of dimension 2*2.

    """
    a_prime = t[2 * n - 1][0] * matrix_x + t[2 * n - 1][1] * matrix_y + t[2 * n - 1][2] * matrix_z
    return a_prime


def mermin(n, t):
    """ Calculates the Mermin polynomial "mn" by recurrence.

    Example :
        >>> mermin(2, [[0.62, -0.035, -0.78],[-0.22, 0.12, -0.97],[0.45, 0.68, -0.57],[0.2, 0.16, -0.97]])
        [[ 0.4066 +0.j       -0.37475+0.5798j   -0.5214 -0.05095j   0.1575 -0.206825j]
         [-0.37475-0.5798j   -0.4066 +0.j        0.1905 +0.199575j   0.5214 +0.05095j]
         [-0.5214 +0.05095j   0.1905 -0.199575j -0.4066 +0.j   0.37475-0.5798j]
         [ 0.1575 +0.206825j  0.5214 -0.05095j   0.37475+0.5798j   0.4066 +0.j]]

    :parameter int n: The number of qubits.
    :parameter: list(int) t : A list of 2 * n lists containing three (03) floats.
    :return: np.array(complex) -- The polynomial "mn" of dimension 2^n * 2^n.

    """
    if n == 1:
        mn = a_matrix(1, t)
        mn_prime = a_prime_matrix(1, t)
    else:
        mn = 0.5 * (np.kron(mermin(n - 1, t),
                            (a_matrix(n, t) + a_prime_matrix(n, t)))
                    + np.kron(mermin_prime(n - 1, t),
                              (a_matrix(n, t) - a_prime_matrix(n, t))))
    return mn


def mermin_prime(n, t):
    """ Calculates the Mermin polynomial "mn'" by recurrence from matrices "a" and "a'"

    Example :
        >>> mermin_prime(2, [[0.62, -0.035, -0.78],[-0.22, 0.12, -0.97],[0.45, 0.68, -0.57],[0.2, 0.16, -0.97]])
        [[ 0.9029+0.j       -0.21775+0.2046j    0.0454+0.0854j  -0.2085+0.210225j]
         [-0.21775-0.2046j   -0.9029+0.j       -0.0895-0.296975j  -0.0454-0.0854j]
         [ 0.0454-0.0854j   -0.0895+0.296975j -0.9029+0.j   0.21775-0.2046j]
         [-0.2085-0.210225j -0.0454+0.0854j    0.21775+0.2046j   0.9029+0.j]]

    :parameter int n: The number of qubits.
    :parameter: list(int) t : A list of 2 * n lists containing three (03) floats.
    :return: np.array(complex) -- The polynomial "mn'" of dimension 2^n * 2^n.

    """
    if n == 1:
        mn = a_matrix(1, t)
        mn_prime = a_prime_matrix(1, t)
    else:
        mn_prime = 0.5 * (np.kron(mermin_prime(n - 1, t),
                                  (a_prime_matrix(n, t) + a_matrix(n, t))) +
                          np.kron(mermin(n - 1, t),
                                  (a_prime_matrix(n, t) - a_matrix(n, t))))
    return mn_prime


def mu_calculation(mn, vector):
    """ Calculates "MU", the result of the calculation of the vector with the mermin polynomial.

    :parameter np.array(complex) mn: The Mermin polynomial "mn".
    :parameter list(int) vector: The state vector.
    :return: float -- The value of the calculation.

    """
    v_dagger = np.conj(vector).T  # The transpose conjugate of the vector
    mu = abs(vector.dot(mn.dot(v_dagger)))
    return mu


def xbest_calculation(n, step, step_min, c_max, vector):
    """ Maximizes Mu. The algorithm used here is called the Random walk method.
        The principle is simple. We randomly generate first parameters which are
        used to calculate Mu. The first value of Mu is called Mu0. Then, we
        calculate new parameters based on the previous ones and a variable
        called the descent step. With the new parameters, we calculate a new
        value of Mu. If this value is better than the previous one, we keep it
        and continue the researches until a counter is at its maximum value. The
        goal here is to take a big circle of research scope and to reduce it (by
        decreasing the decent step) more and more until the maximum value of Mu
        is found.

    :parameter int n: The number of qubits.
    :parameter int step: The initial value of the descent step.
    :parameter float step_min: The minimum value of the descent step (which is the length of the radius).
    :parameter int c_max: The maximum value of the counter.
    :parameter list(int) vector: The vector for the calculation of Mu.
    :return: float -- The value that is, expectedly, the optimal maximal value of Mu.

    """
    c = 0

    x0_tab = first_coefficients_generation(n)
    m0 = mermin(n, x0_tab)
    first_mu = mu_calculation(m0, vector)

    # Initialisation
    xbest_tab = x0_tab
    mu_best = first_mu

    while step > step_min:
        while c < c_max:
            tableau_nouv_coef_mk = new_coefficients_generation(n, xbest_tab, step)
            mk = mermin(n, tableau_nouv_coef_mk)
            mu_new_value = mu_calculation(mk, vector)
            if mu_new_value > mu_best:
                xbest_tab = tableau_nouv_coef_mk
                mu_best = mu_new_value
                c = 1
            else:
                c += 1
        c = 1
        step = step / 2

    return mu_best
