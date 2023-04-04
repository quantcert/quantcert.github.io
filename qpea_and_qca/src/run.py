from qpea import *


def tests(n, test=True):
    """ Runs the tests. A default setting is set for simple tests. But, all the calculations may be made by putting the
     variable "test" to false.

    :param int n: The number of qubits.
    :param bool test : If activated, simples tests will be run.

    """
    if test:
        n = 3
        case = "Z"
        n_second_register = 1
        run(n, n_second_register, case, test=True)
    else:
        n = 3
        n_second_register = 1
        while n < 6:
            for case in ["Z", "G"]:
                run(n, n_second_register, case, test=False)
            clear = True  # To separate the different figures
            n += 1


# For test running
tests(3)
