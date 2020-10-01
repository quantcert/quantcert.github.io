#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""There are three format used for the algorithms: 

1.  A flat list of coefficients, organized as such: 
    `[x1,y1,z1, x2,y2, ..., xn,yn,zn, x'1,y'1,z'1, x'2,y'2, ..., x'n,y'n,z'n]`.
    This format is called the *unpacked coefficients* and is used for the QFT 
    optimization.
2.  A list of coefficients grouped by families of operators:
    `[[[x1,y1,z1], [x2,y2, ..., [xn,yn,zn]], [[x'1,y'1,z'1], [x'2,y'2, ..., [x'n,y'n,z'n]]]`
    in other words, you have the whole `a` family and the the whole `a'` family,
    and in a family, you have `a1`, `a2`, and so on... Each `a` is described by
    it's three coefficients: `x`, `y` and `z`.
    This format is called *packed coefficients* and is used to easily manipulate
    coefficients.
3.  A list of coefficients grouped by operator:
    `[[x1,y1,z1], [x'1,y'1,z'1], [x2,y2,z2], [x'2,y'2,z'2], ...]`
    in other words, the list is formed as such: `[a1, a'1, a2, a'2, ...]`
    This format was previously used for evaluation in Qiskit, allowing for a
    simpler data flow. But it has the inconvenient of being less true to the 
    maths behind all this so it has been dropped. 
    This format is called *mixed*

With the functions of this module, one may switch between *1.* and *2.* and 
between *2.* and *3.*, allowing tho switch freely between the three formats.
"""


def coefficients_format_unpacked_to_packed(_a_a_prime_coeffs):
    r"""
    Packs a list of elements in two lists of lists of three elements

    Example:
        >>> coefficients_format_unpacked_to_packed([1,2,3,4,5,6,7,8,9,10,11,12])
        ([[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]])

    :param list[any] _a_a_prime_coeffs: List of elements (unpacked coefficients).
    :returns: tuple[list[list[any]]] -- Lists of lists of elements as described 
        above (packed coefficients). 
    """
    _a_a_prime_coeffs_packed = [[_a_a_prime_coeffs[3*_i], _a_a_prime_coeffs[3*_i+1], 
        _a_a_prime_coeffs[3*_i+2]] for _i in range(int(len(_a_a_prime_coeffs)/3))]
    _a_coeffs = _a_a_prime_coeffs_packed[:int(len(_a_a_prime_coeffs_packed)/2)]
    _a_prime_coeffs = _a_a_prime_coeffs_packed[int(len(_a_a_prime_coeffs_packed)/2):]
    return _a_coeffs, _a_prime_coeffs


def coefficients_format_packed_to_unpacked(_a_coeffs, _a_prime_coeffs):
    r"""
    Unpacks two lists of lists of three elements to one list of elements

    Example:
        >>> coefficients_format_packed_to_unpacked([[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]])
        [1,2,3,4,5,6,7,8,9,10,11,12]

    :param list[list[any]] _a_coeffs, _a_prime_coeffs: Lists of lists of 
        elements as described above (packed coefficients). 
    :returns: list[any] -- List of elements (unpacked coefficients).
    """
    _a_a_prime_coeffs = [k for element in _a_coeffs+_a_prime_coeffs for k in element]
    return _a_a_prime_coeffs
    

def coefficients_format_mixed_to_packed(_a_a_prime_coeffs):
    r"""Format the coefficients in the shape previously used for evaluation

    Example:
        >>> coefficients_format_mixed_to_packed([[1, 2, 3], [7, 8, 9], [4, 5, 6], [10, 11, 12]])                   
        ([[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]])

    :param list[list[any]] _a_a_prime_coeffs: List of lists of elements as 
        described above (mixed coefficients). 
    :returns: tuple(list[list[any]]) -- Tuple of list of list of elements (packed
        coefficients).
    """
    n2 = int(len(_a_a_prime_coeffs)/2)
    _a_coeffs = [_a_a_prime_coeffs[2*i] for i in range(n2)]
    _a_p_coeffs = [_a_a_prime_coeffs[2*i+1] for i in range(n2)]
    return _a_coeffs, _a_p_coeffs


def coefficients_format_packed_to_mixed(_a_coeffs, _a_prime_coeffs):
    r"""Format the coefficients in the shape now used for evaluation

    Example:
        >>> coefficients_format_packed_to_mixed([[1,2,3],[4,5,6]], [[7,8,9],[10,11,12]])                           
        [[1, 2, 3], [7, 8, 9], [4, 5, 6], [10, 11, 12]]

    :param list[list[any]] _a_coeffs, _a_prime_coeffs: Lists of lists of 
        elements as described above (packed coefficients). 
    :returns: list[any] -- List of list of elements (mixed coefficients).
    """
    _a_a_prime_coeffs = []
    for index in range(len(_a_coeffs)):
        _a_a_prime_coeffs.append(_a_coeffs[index])
        _a_a_prime_coeffs.append(_a_prime_coeffs[index])
    return _a_a_prime_coeffs