import numpy as np


def convert_in_binary(number_to_convert, number_of_bits=0):
    """ Converts an int into a string containing its bits

    Example :
        >>> convert_in_binary(5,3)
        101
        >>> convert_in_binary(5,5)
        00101

    :param int number_to_convert: The number that is going to be converted.
    :param int number_of_bits: The number of bits required.
    :returns: str -- The converted number.
    """
    return format(number_to_convert, "b").zfill(number_of_bits)


def mermin_coeffs_to_U3_coeffs(x, y, z):
    """ Generates the coefficients of the U3 gate from mermin coefficients:
        x*X + y*Y + z*Z = U3(theta, phi, -phi-pi)

    :param float x: The coefficient alpha for the matrix X.
    :param float y: The coefficient beta for the matrix Y.
    :param float z: The coefficient gamma for the matrix Z.
    :returns: (float, float) -- The two angles of U3 gate.
    """
    theta = np.arccos(z)
    phi = 0 if np.sin(theta) == 0 else np.arccos(x/np.sin(theta))
    if y/np.sin(theta) < 0:
        phi = - phi
    return theta, phi


def U3_gates_placement(n, n_measure, a_a_p_coeffs, circuit):
    """ Places the U3 gates according to the mermin_IBM monomial

    :param int n: the size of the register to be evaluated
    :param int n_measure: The measure to be performed. Dictates whether a_i or
        a'_i is used on each wire
    :param list[list[real]] a_a_p_coeffs: Contains the list of coefficients for 
        a_i and a'_i in the packed shape
    :param QuantumCircuit circuit: The circuit on which the measures are
        appended
    :returns: None
    """
    mermin_monomial_description = convert_in_binary(n_measure, n)
    a_coeffs, a_p_coeffs = a_a_p_coeffs
    for k in range(n):
        if mermin_monomial_description[k] == '0':
            theta, phi = mermin_coeffs_to_U3_coeffs(*a_coeffs[k])
        else:
            theta, phi = mermin_coeffs_to_U3_coeffs(*a_p_coeffs[k])
        circuit.u3(theta, np.pi, -phi-np.pi, k)
    circuit.measure(range(n), range(n))
