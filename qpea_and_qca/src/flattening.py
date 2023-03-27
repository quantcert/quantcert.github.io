import numpy as np
from itertools import combinations
from math import comb, floor


def maximum_of_flattening(n):
    """ Calculates the maximum of flattening types that can be done. If the reuslt of the calculation is
    a float, then, this number turns out to be the largest integer less or equal to our number.

    Examples:
        >>> maximum_of_flattening(4)
        2
        >>> maximum_of_flattening(7)
        3

    :param int n: The number of qubits.
    :returns: int -- The maximum number of flattening types.
    """
    return floor(n / 2)


def maximum_of_matrices(n, type_of_flattening):
    """ Calculates the maximum of matrices for a type of flattening. If the type of flattening is a divider of the
    number of qubits, then, the first result is divided by 2.

    Examples:
        >>> maximum_of_matrices(2, 1)
        1.0
        >>> maximum_of_matrices(4, 2)
        3.0

    :param int n: The number of qubits.
    :param int type_of_flattening : The type of flattening.
    :returns: float -- The maximum number of matrices.
    """
    max_matrices = comb(n, type_of_flattening)
    if n % type_of_flattening == 0:
        max_matrices = max_matrices / 2
    return max_matrices


def convert_in_binary(number, bits=0):
    """ Converts an int into a string containing the number in bits.

    Examples:
        >>> convert_in_binary(5,3)
        '101'
        >>> convert_in_binary(5,5)
        '00101'

    :param int number: the number that is going to be converted.
    :param int bits: the number of bits required.
    :returns: String -- The converted number.
    """
    return format(number, "b").zfill(bits)


def putting_in_list(number):
    """ Puts every figure of the number into a list type.

    Example:
        >>> putting_in_list("000")
        [0, 0, 0]

    :param string number: The number to split in a string format.
    :returns: list(int) -- The list with every digit of the number in a case.
    """
    number_in_caracters = str(number)
    number_in_list = []
    for caracter in number_in_caracters:
        number_in_list.append(int(caracter))
    return number_in_list


def list_of_combinaisons(n, type_of_flattening):
    """ Realizes the possible combinations for a type of flattening.

    Examples:
        >>> list_of_combinaisons(5, 1)
        [(1,), (2,), (3,), (4,), (5,)]
        >>> list_of_combinaisons(4, 2)
        [(1, 2), (1, 3), (2, 3)]

    :param int n: The number of qubits.
    :param int type_of_flattening : The current type of flattening.
    :returns: list[list(int)] -- The list of every possible combination.
    """
    if type_of_flattening == 1:
        n = n + 1
    l = list(range(1, n))
    combinaisons_list = []
    for i in combinations(l, type_of_flattening):
        combinaisons_list.append(i)
    return combinaisons_list


def list_of_combinaisons_for_invariant(n):
    """ Realizes the possible combinations for a type of flattening.

    Example:
        >>> list_of_combinaisons(4)
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    :param int n: The number of qubits.
    :returns: list[list(int)] -- The list of every possible combination.
    """
    L = list(range(n))
    combinaisons_list = []
    for i in combinations(L, 2):
        combinaisons_list.append(i)
    return combinaisons_list


def flattening_and_rank_calculation(n, T):
    """ Construct the matrices of every flattening type and calculates their rank.

    Example:
        >>> flattening_and_rank_calculation(3, [1, 0, 0, 0, 0, 0, 0, 1])
        ([2, 2, 2], [[array([[1, 0, 0, 0], [0, 0, 0, 1]]), array([[1, 0, 0, 0], [0, 0, 0, 1]]), array([[1, 0, 0, 0], [0, 0, 0, 1]])]])

    :param int n: The number of qubits.
    :param list(int) T : A list filled with the coefficients corresponding to the vector.
    :returns: list(list(int)) -- A list of all the matrices after all flatennings.
    """
    q = maximum_of_flattening(n)
    all_coef_matrices = []
    for m in range(1, q + 1):
        list_combinaisons = list_of_combinaisons(n, m)
        coef_matrices = []
        for i in range(len(list_combinaisons)):
            l = list_combinaisons[i]
            M = [[]]
            for h in range(2 ** m - 1):
                M.append([])
            for j in range(2 ** n):
                bin_j = convert_in_binary(j, n)
                v = putting_in_list(bin_j)
                w = ""
                for k in range(m):
                    w += str(v[l[k] - 1])
                M[int(w, 2)].append(T[j])
            M = np.array(M)
            coef_matrices.append(M)
        all_coef_matrices.append(coef_matrices)
    return all_coef_matrices
