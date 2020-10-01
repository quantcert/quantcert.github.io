from .mermin_polynomials import *


def convert_in_binary(number_to_convert, number_of_bits=0):
    """ Converts an int into a string containing its bits.

    Example:
        >>> convert_in_binary(5,3)
         101
         and
        >>> convert_in_binary(5,3)
        00101

    :param int number_to_convert: the number that is going to be converted.
    :param int number_of_bits: the number of bits required.
    :returns: String number : The converted number.
    """
    return format(number_to_convert, "b").zfill(number_of_bits)


def putting_in_list(number):
    """ Puts every figure of the number in a list.

    Example:
        >>> putting_in_list(000)
         [0, 0, 0]

    :param int number: The number to split.
    :returns: list(int) -- The list with every digit of the number in a case.
    """
    number_in_caracters = str(number)
    number_in_list = []
    for caracter in number_in_caracters:
        number_in_list.append(int(caracter))
    return number_in_list


def states_formation(n):
    """ Calculates every state for n qubits.

    Example:
        >>> states_formation(3)
         [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

    :param int n: The number of qubits.
    :returns: list(list(int)) -- A list of all the possible states of the qubits 
        which are also contained in a list.
    """
    list_of_states = []
    maximum = 2 ** n  # Maximum
    for counter in range(maximum):
        list_of_states.append(putting_in_list(convert_in_binary(counter, n)))  # Putting numbers in a tab
    return list_of_states


def state_vector_initialisation(n):
    """ Initializes the state vector; which is to create an array with the size 
        of 2 to the power of n. Every number in the array is equal to 1 over the
        square root of 2 to the power of n.

    Example:
        >>> state_vector_initialisation(3)
        [0.35355339 0.35355339 0.35355339 0.35355339 0.35355339 0.35355339
         0.35355339 0.35355339]

    :param int n: The number of qubits.
    :returns: list(int) -- The initialized state vector.
    """
    state_vector_init = np.zeros(2 ** n)  # Array full of zeros
    maximum = 2 ** n  # Maximum
    for counter in range(maximum):
        state_vector_init[counter] = (1 / (2 ** n) ** 0.5)
    return state_vector_init


def hyperedges_computation(n, state_vector, hyperedges):
    """ Puts the phases in the right places. Wherever there is an edge or an 
        hyperedge between some vertices, a minus sign is put where those 
        vertices are all in state 1.

    Example:
        >>> hyperedges_computation(2, [0.5, 0.5, 0.5, 0.5], [[0,1]])
        [0.5, 0.5, 0.5, -0.5]

    :param int n: The number of qubits.
    :param list(int) state_vector: The initialized state vector.
    :param list[list[int] hyperedges: a list containing the lists of the 
        vertices which are linked by an hyperedge.
    :returns: list(int) -- The correct state vector.
    """
    max_hyperedges = len(hyperedges)
    max_state_vector = 2 ** n
    # Writing every state in a binary form
    state_vector_en_binaire = states_formation(n)  
    for counter_hyperedges in range(max_hyperedges):
        for counter_state_vector in range(max_state_vector):
            state = state_vector_en_binaire[counter_state_vector]
            # Throws True if the state matches an edge and False if not
            verified_state = corresponding_state_determination(
                n, state, hyperedges[counter_hyperedges])
            if verified_state:
            # If True is returned, apply a phase to the corresponding element
                state_vector[counter_state_vector] *= -1
    return state_vector


def corresponding_state_determination(n, state, hyperedges):
    """ Determines the states corresponding to an hyperedge. That is the states 
        where the vertices linked by the hyperedge are equal to 1.

    :param int n: The number of qubits.
    :param list(int) state: The state of the vector.
    :param list(list(int)) hyperedges: a list containing the lists of the 
        vertices which are linked by an hyperedge.
    :param: boolean -- True if the state is in an hyperedge.
    """
    maximum = len(hyperedges)
    verified_state = 0
    for counter in range(maximum):
        x = hyperedges[counter]
        # n - x - 1 :returns the index where the digit (either 0 or 1) of the 
        # state, must be checked
        if state[n - x - 1] == 1:
            verified_state += 1
    # If the sum is greater than the number of hyperedges, it means that the 
    # state checks the condition of the hyperedge
    return verified_state >= maximum
