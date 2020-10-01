import numpy as np
import random as rd

# The Pauli matrices
matrix_x = np.array([[0, 1], [1, 0]])
matrix_y = np.array([[0, complex(0, -1)], [complex(0, 1), 0]])
matrix_z = np.array([[1, 0], [0, -1]])


def first_coefficients_generation(n):
    """ Generates the very first coefficients for the calculation of MU.

    Example :
        >>> first_coefficients_generation(2)
        [[ 0.24006446 -0.97020025 0.03287126]
         [0.72092088 -0.59054414  0.36267162]
         [-0.76022821 0.64723032 -0.056089]
         [-0.0278048  0.48298397 -0.87518763]]

    :param int n: The number of qubits.
    :returns: np.array(list(float)) -- The table of list of the coefficient 
        taken randomly.
    """
    first_coefficient_list = []  # Empty list
    maximum = 2 * n  # Maximum of iterations

    # Filling the list
    for counter in range(maximum):
        # Production of parameters
        alpha_init = rd.uniform(-1, 1)
        beta_init = rd.uniform(-1, 1)
        gamma_init = rd.uniform(-1, 1)
        norm = np.sqrt(alpha_init ** 2 + beta_init ** 2 + gamma_init ** 2)
        first_coefficient_list.append([alpha_init / norm, beta_init / norm, gamma_init / norm])

    first_coefficients = np.array(first_coefficient_list)
    return first_coefficients


def a_matrix(n, t):
    """ Constitutes the matrix a of the Mermin polynomial.

    :param int n: The number of qubits.
    :param np.array(list(float)) t: The table of coefficients.
    :returns: np.array(complex) -- The matrice a of mn.
    """
    a = t[2 * n - 2][0] * matrix_x + \
        t[2 * n - 2][1] * matrix_y + \
        t[2 * n - 2][2] * matrix_z  # a = x*X + y*Y + z*Z
    return a


def a_prime_matrix(n, t):
    """ Constitutes the matrix a' of the Mermin polynomial.

    :param int n: The number of qubits.
    :param np.array(list(float)) t: The table of coefficients.
    :returns: np.array(complex) -- The matrice a' of mn.
    """
    # a' = x'*X + y'*Y + z'*Z
    a_prime = t[2 * n - 1][0] * matrix_x + t[2 * n - 1][1] * matrix_y + t[2 * n - 1][2] * matrix_z
    return a_prime


def mermin(n, t):
    """ Calculates the Mermin polynomial `mn`.

    :param int n: The number of qubits.
    :param np.array(list(float)) t: The table of coefficients.
    :returns: np.array(complex) -- The Mermin polynomial `mn`.
    """
    if n == 1:
        # M1 is the first Mermin polynomial
        mn = a_matrix(1, t)
    else:
        # Formula of recurrence
        mn = 0.5 * (np.kron(mermin(n - 1, t),
                            (a_matrix(n, t) + a_prime_matrix(n, t)))
                    + np.kron(mermin_prime(n - 1, t),
                              (a_matrix(n, t) - a_prime_matrix(n, t))))
    return mn


def mermin_prime(n, t):
    """ Calculates the Mermin polynomial `mn'`.

    :param int n: The number of qubits
    :param np.array(list(float)) t: The table of coefficients
    :returns: np.array(complex) -- The Mermin polynomial `mn'`
    """
    if n == 1:
        mn_prime = a_prime_matrix(1, t)
    else:
        # Recurrence formula
        mn_prime = 0.5 * (np.kron(mermin_prime(n - 1, t),
                                  (a_prime_matrix(n, t) + a_matrix(n, t))) +
                          np.kron(mermin(n - 1, t),
                                  (a_prime_matrix(n, t) - a_matrix(n, t))))
    return mn_prime


def mu_calculation(mn, mn_prime, vector, type_of_mu):
    """ Calculates MU, the value of the calculation of the vector with the 
        mermin polynomial.

    :param np.array(complex) mn: The Mermin polynomial `mn`.
    :param np.array(complex) mn_prime: The Mermin polynomial `mn'`.
    :param list(int) vector: The state vector.
    :param type_of_mu bool: If False, the classical calculation will be made. 
        If not, another method is used.
    :returns: float -- The value of the calculation.
    """
    v_dagger = np.conj(vector).T
    if type_of_mu:
        # <V|Mn|V>² + <V|M'n|V>²
        mu = vector.dot(mn.dot(v_dagger)) ** 2 + vector.dot(mn_prime.dot(v_dagger)) ** 2  
    else:  
        # |<V|Mn|V>|
        mu = abs(vector.dot(mn.dot(v_dagger)))
    return abs(mu)


def new_coefficients_generation(n, old_coefficients, alpha):
    """ Random generation of new parameters for MU maximization

    :param int n: The number of qubits
    :param np.array(list(float)) old_coefficients: The table of the coefficients 
        that didn't maximize Mu
    :param int alpha: The value of the descent step (used in the random walk 
        method)
    :returns:  np.array(list(float)) -- Table of new coefficients
    """
    maximum = 2 * n
    new_coefficients_list = []
    for counter in range(maximum):
        # Parameters generation
        alpha_init = old_coefficients[counter][0] + (rd.uniform(-1, 1) * alpha)
        beta_init = old_coefficients[counter][1] + (rd.uniform(-1, 1) * alpha)
        gamma_init = old_coefficients[counter][2] + (rd.uniform(-1, 1) * alpha)
        norm = np.sqrt(alpha_init ** 2 + beta_init ** 2 + gamma_init ** 2)
        new_coefficients_list.append(
            [alpha_init / norm, beta_init / norm, gamma_init / norm])
    new_coefficients = np.array(new_coefficients_list)
    return new_coefficients


def mu_file_saving(n, file_path, first_mu, first_parameters, maximal_mu, maximisation_parameters, type_of_mu):
    """ Creates a file to write the various calculation parameters

    :param int n: The number of qubits.
    :param string file_path: The path where the file is to be saved
    :param float first_mu: The value of the calculation of the vector with the 
        mermin polynomial
    :param np.array(list(float)) first_parameters: The first coefficients of Mu 
        calculation
    :param float maximal_mu: The value of the maximal Mu calculated

    """
    with open(file_path, "w") as my_file:
        if type_of_mu == 0:
            my_file.write("MU formula of calculation : |<V|Mn|V>|")
        else:
            my_file.write("MU formula of calculation : <V|Mn|V>² + <V|M'n|V>²")
        my_file.write("\n")
        my_file.write("\n" + "Number of qubits :")
        my_file.write("\n" + str(n))
        my_file.write("\n")
        my_file.write("\n" + "First value of MU : ")
        my_file.write("\n" + str(first_mu))
        my_file.write("\n" + "First parameters : ")
        my_file.write("\n" + str(first_parameters))
        my_file.write("\n")
        my_file.write("\n" + "Maximal value of MU  : ")
        my_file.write("\n" + str(maximal_mu))
        my_file.write("\n" + "Maximal MU parameters : ")
        my_file.write("\n" + str(maximisation_parameters))
        my_file.write("\n")


def xbest_calculation(n, type_of_mu, alpha, alpha_minimum, c_maximum, vector, file_path, saving_file=True):
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

    :param int n: The number of qubits.
    :param type_of_mu bool: If False is specified, the classical calculation 
        will be made. If not, another method is used.
    :param int alpha: The value of the descent step.
    :param int alpha_minimum: The minimum value of the descent step (which is 
        the length of the radius).
    :param int c_maximum: The maximum value of the counter.
    :param list(int) vector: The vector for the calculation of Mu.
    :param string file_path: The path where the file is to be saved.
    :param boolean saving_file: If set to True, a file will be created / 
        overloaded with the information about the calculation of Mu. If not, 
        only the calculations are made.
    :returns: np.array(list(float)) -- The array that contains the parameters 
        that maximizes Mu.
    """
    c = 0

    # First iteration
    x0_tab = first_coefficients_generation(n)
    m0 = mermin(n, x0_tab)
    mnp0 = mermin_prime(n, x0_tab)
    first_mu = mu_calculation(m0, mnp0, vector, type_of_mu) 

    # Initial parameters
    xbest_tab = x0_tab
    mu_best = first_mu

    # Former iterations
    while alpha > alpha_minimum:
        while c < c_maximum:
            # Generation of the new parameters
            nouv_coef_mk = new_coefficients_generation(n, xbest_tab, alpha)
            mk = mermin(n, nouv_coef_mk)
            mkp = mermin_prime(n, nouv_coef_mk)
            mu_new_value = mu_calculation(mk, mkp, vector, type_of_mu)
            if mu_new_value > mu_best:
                xbest_tab = nouv_coef_mk
                mu_best = mu_new_value
                c = 1
            else:
                c += 1
        c = 1
        alpha = alpha / 2
    if saving_file:
        mu_file_saving(n, file_path, first_mu, x0_tab, mu_best, xbest_tab, type_of_mu)
    return xbest_tab
