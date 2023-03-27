from rank import *
from math import sqrt


def printing_cms(n, vector):
    """Calculates the cumulation of the coefficient matrices values which is an invariant. In fact, at first, from the
    input vector, all the flattenings are made. Then, a calculation is done from each of them. In the end, the
    invariant is just the sum of all the results.

    Example:
        >>> printing_cms(3, [1 / sqrt(2), 0, 0, 0, 0, 0, 0, 1 / sqrt(2)])
        1

    :param int n: The number of qubits.
    :param list(int) vector : A list filled with the coefficients corresponding to the vector.
    :returns: int invariant_total : The invariant.
    """
    all_cms = flattening_and_rank_calculation(n, vector)
    number_of_flattening_types = maximum_of_flattening(n)
    invariant_total = 0
    un_nombre = 0

    for type_of_flattening in range(1, number_of_flattening_types + 1):
        invariant_col = 0
        number_of_matrices_per_flattening_type = len(list_of_combinaisons(n, type_of_flattening))
        for i in range(number_of_matrices_per_flattening_type):
            one_cm = all_cms[type_of_flattening - 1][i]
            one_cm = one_cm.T
            columns_length = one_cm.shape[0]
            list_of_columns_combinations = list_of_combinaisons_for_invariant(columns_length)
            for k in range(len(list_of_columns_combinations)):
                colonne_i = one_cm[list_of_columns_combinations[k][0]]
                colonne_j = one_cm[list_of_columns_combinations[k][1]]
                invariant_col_inter = ((np.linalg.norm(colonne_i) ** 2 * np.linalg.norm(colonne_j) ** 2) -
                                       abs((colonne_i.dot(np.conj(colonne_j.T))) ** 2))
                invariant_col = invariant_col + invariant_col_inter
            if invariant_col < 0:
                invariant_col *= -1
            invariant_cm = sqrt(4 * invariant_col)
            invariant_col = 0
            invariant_total = invariant_total + invariant_cm
        un_nombre = un_nombre + number_of_matrices_per_flattening_type

    invariant_total = invariant_total / un_nombre
    return invariant_total
