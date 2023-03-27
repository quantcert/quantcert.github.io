from builtins import complex
from math import e, pi
from datetime import datetime
import matplotlib.pyplot as plt

from coefficient_matrices import *
from mermin_polynomials import *
from gme import *

# Variables for the random walk algorithm
type_of_mu = 0
alpha = 5
alpha_min = 0.00001
c_max = 200


def cas_matrice_z(n, angles, basis):
    """ Generates the vectors depending on the basis and also the angles for the quantum phase estimation algorithm.

    :param int n: The number of qubits.
    :param np.array(float) angles : The table of the angle values.
    :param int basis: The value of the number of qubits on the second register.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithm will be.

    """
    if basis == 1:
        psi = psi_calculation_qpea_base_0(n, angles)
    else:
        psi = psi_calculation_qpea_base_00(n, angles)
    return psi


def cas_matrice_g(n, angles, basis):
    """ Generates the vectors depending on the basis and also the angles for the quantum counting algorithm.

    :param int n: The number of qubits.
    :param np.array(float) angles : The table of the angle values.
    :param int basis: The value of the number of qubits on the second register.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithm will be.

    """
    if basis == 1:
        psi = psi_calculation_qca_base_0(n, angles)
    else:
        psi = psi_calculation_qpea_base_00(n, angles)
    return psi


def angles_determination(min, max, step):
    """ Generates a table of different values that corresponds to an angle which will be involved in the algorithms
    calculations and variations.

    :param int min: The minimal value of the angle.
    :param int max: The maximal value of the angle.
    :param float step: The step of variation of the angle.
    :returns: np.array(float) -- The table of float values of the angle.

    """
    angles = np.arange(min, max, step)
    return angles


def psi_calculation_qca_base_0(n, thetas):
    """ Computes the vector for each angle for the quantum counting algorithm when there are only one qubit on the
    second register. This means that the total number of qubits of all the system is n + 1.

    :param int n: The number of qubits.
    :param np.array(float) thetas : The table of float values of the angle.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithm will be.
    performed.

    """
    norm = 1 / sqrt(2 ** n)
    max = 2 ** (n - 1)
    liste = []
    for theta in thetas:
        psi_liste = []
        for x in range(max):
            zero = -2 * 1j * pi * theta * (x + 0.5)
            psi_liste.append((e ** zero) * norm)
            un = 2 * 1j * pi * theta * (x - 0.5)
            psi_liste.append((e ** un) * norm)
        liste.append(psi_liste)
    psi = np.array(liste, dtype=object)
    return psi


def psi_calculation_qca_base_00(n, thetas):
    """ Computes the vector for each angle for the quantum counting algorithm when there are two qubits on the
    second register. This means that the total number of qubits of all the system is n + 2.

    :param int n: The number of qubits.
    :param np.array(float) thetas : The table of float values of the angle.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithm will be
    performed.

    """
    norm = 1 / sqrt(2 ** (n + 2))
    max = 2 ** (n - 2)
    liste = []
    for theta in thetas:
        psi_liste = []
        for x in range(max):
            zero = -2 * 1j * pi * theta * (x + 0.5)
            psi_liste.append((e ** zero) * norm)
            psi_liste.append((e ** zero) * norm)
            un = 2 * 1j * pi * theta * (x - 0.5)
            psi_liste.append((e ** un) * norm)
            psi_liste.append((e ** un) * norm)
        liste.append(psi_liste)
    psi = np.array(liste, dtype=object)
    return psi


def psi_calculation_qpea_base_0(n, phis):
    """ Computes the vector for each angle for the quantum phase estimation algorithm when there are only one qubit on
    the second register. This means that the total number of qubits of all the system is n + 1.

    :param int n: The number of qubits.
    :param np.array(float) phis : The table of float values of the angle.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithm will be
    performed.

    """
    norm = 1 / sqrt(2 ** n)
    max = 2 ** (n - 1)
    liste = []
    for phi in phis:
        psi_liste = []
        for x in range(max):
            psi_liste.append(norm)
            psi_liste.append(norm * (e ** (complex(0, 2 * pi * x * phi))))
        liste.append(psi_liste)
    psi = np.array(liste, dtype=object)
    return psi


def psi_calculation_qpea_base_00(n, phis):
    """ Computes the vector for each angle for the quantum phase estimation algorithm when there are two qubits on
    the second register. This means that the total number of qubits of all the system is n + 2.

    :param int n: The number of qubits.
    :param np.array(float) phis : The table of float values of the angle.
    :returns: np.array(complex) -- A table of vectors corresponding to each angle value on which the algorithms will be
    performed.

    """
    norm = 1 / sqrt(2 ** n)
    max = 2 ** (n - 2)
    liste = []
    for phi in phis:
        psi_liste = []
        for x in range(max):
            psi_liste.append(norm)
            psi_liste.append(norm * ((-1) ** x))
            psi_liste.append(norm * (e ** (-2j * pi * x * phi)))
            psi_liste.append(norm * (e ** (2j * pi * x * phi)))
        liste.append(psi_liste)
    psi = np.array(liste, dtype=object)
    return psi


def calcul_invariant_cm(n, psi):
    """ Calculates the value of the invariant based on the coefficient matrices.

    :param int n: The number of qubits.
    :param np.array(float) psi : The table of vectors for the algorithm.
    :returns: np.array(complex) -- The list of the coefficients calculated.

    """
    liste_cms = []
    for i in range(len(psi)):
        cm = printing_cms(n, psi[i])
        liste_cms.append(cm)
    return liste_cms


def calcul_invariant_mu(n, angles, psi):
    """ Calculates the value of the invariant based on Mermin polynomials.

    :param int n: The number of qubits.
    :param np.array(float) angles : The table of the angle values.
    :param np.array(float) psi : The table of vectors for the algorithm.
    :returns: np.array(complex) -- The list of the coefficients calculated.

    """
    liste_mus = []

    for i in range(len(angles)):
        liste_mus.append(xbest_calculation(n, alpha, alpha_min, c_max, psi[i]))
    return liste_mus


def calcul_invariant_gme(n, angles, psi):
    """ Calculates the value of the invariant based on Mermin polynomials.

    :param int n: The number of qubits.
    :param np.array(float) angles : The table of the angle values.
    :param np.array(float) psi : The table of vectors for the algorithm.
    :returns: np.array(complex) -- The list of the coefficients calculated.

    """
    liste_gmes = []

    for i in range(len(angles)):
        liste_gmes.append(gme_best_calculation(n, psi[i], alpha, alpha_min, c_max))
    return liste_gmes


def run(n_first_register, n_second_register, case, cm=True, mu=True, gm=True):
    """ Generates a table of different values that corresponds to an angle which will be involved in the algorithms
    calculations and variations.

    :param int n_first_register: The number of qubits on the first register of the system.
    :param int n_second_register: The value of the number of qubits on the second register.
    :param str case: The value that determines the type of algorithm : "G" for qca or "Z" for qpea.
    :param boolean cm : This determines if the invariant based on the coefficient matrices should be calculated or not
    :param boolean mu : This determines if the invariant based on the Mermin polynials should be calculated or not
    :param boolean gm : This determines if the geometric measure invariant should be calculated or not
    :returns: np.array(float) -- The table of float values of the angle.

    """
    # Initialisations
    cms = []
    mus = []
    gmes = []
    psis = []
    figure_name = ""
    figure_title = ""

    # Deletes the old records
    plt.clf()

    # Angles
    angles = angles_determination(min=0, max=1, step=0.01)

    # Total number of qubits definition
    n_total = n_first_register + n_second_register

    # Basis choice
    if case == "Z":
        figure_name = "QPEA"
        psis = cas_matrice_z(n_total, angles, n_second_register)
        figure_title = "Quantum Phase Estimation Algorithm "
    if case == "G":
        figure_name = "QCA"
        psis = cas_matrice_g(n_total, angles, n_second_register)
        figure_title = "Quantum Counting Algorithm "

    partial_name = str(n_first_register) + "+" + str(n_second_register) + "qubits"

    # Figure title
    figure_title += "for " + partial_name
    plt.title(figure_title)

    if cm:
        figure_name += "_CM"
        if not cms:
            cms = calcul_invariant_cm(n_total, psis)
    if mu:
        figure_name += "_MU"
        mus = [0.9999999999067279, 1.0103062222152035, 1.0405665920265939, 1.0889422751646143, 1.1527524278805945, 1.2288389766268188, 1.3138405208539408, 1.4043483587175283, 1.496980728967107, 1.5884231296220788, 1.6754690852864125, 1.755078083330678, 1.8244644421456928, 1.8811919708574691, 1.9233067103801844, 1.9494805496877794, 1.9591481218698696, 1.9526719570212117, 1.9314873490675082, 1.8982328255914078, 1.856830989996974, 1.812463816884072, 1.7713548382279636, 1.7402410919708633, 1.7255335895156148, 1.7320508073458023, 1.7620354373100198, 1.8145278515514063, 1.885676776176094, 1.9694898609990459, 2.0590476587362216, 2.1473898028424565, 2.2281717566662085, 2.2960387077701476, 2.346805761501569, 2.3775198481459032, 2.386450618004572, 2.3730514404152783, 2.337857957727014, 2.282412727350107, 2.209140544981537, 2.1212166619215393, 2.0224569280724434, 1.9171923595909324, 1.8101627493296928, 1.7064220309937372, 1.6112106582550543, 1.5298234744484558, 1.4672424510969109, 1.4277324369209297, 1.4142135620999805, 1.4277324347673055, 1.4672425897845138, 1.529823120689342, 1.6112107722854376, 1.706397114563703, 1.810163018834823, 1.917192903264069, 2.0224570878922603, 2.121216645317109, 2.2091406434618444, 2.2824129499826364, 2.337857920775528, 2.3730517264330504, 2.386451442652399, 2.3775198912375686, 2.346805770061271, 2.2960387074996693, 2.2281717723242105, 2.147389829660271, 2.059047723972027, 1.9694897162155698, 1.8856766120231359, 1.8145398129064096, 1.7620349975820364, 1.732050807057174, 1.7255337868883218, 1.7402460220754987, 1.7713547492766084, 1.8124639943559362, 1.856831077120687, 1.898232855205713, 1.9314873101625183, 1.952671957350256, 1.9591481147363243, 1.9494805983340204, 1.9233059764470963, 1.881191270166835, 1.824464380198904, 1.7550778407042589, 1.6754690439110307, 1.5884231310548558, 1.4969807310340382, 1.404348356864725, 1.313840520917688, 1.2288389760180813, 1.1527524298488945, 1.0889422804377134, 1.0405665867274447, 1.0103062903380824]
        if not mus:
            mus = calcul_invariant_mu(n_total, angles, psis)
    if gm:
        figure_name += "_GME"
        gmes = [9.988188054421698e-11, 0.005170899886881686, 0.020556444565457466, 0.045779310777429694, 0.080224273599112, 0.12305814295608786, 0.17325662410511855, 0.22963710323288644, 0.2908961490792845, 0.35565034746109636, 0.42247898202151435, 0.47934455650004193, 0.49831527703978074, 0.49870809915894976, 0.49088593550877246, 0.48000292329213623, 0.4691605712365421, 0.4604855812965428, 0.4555075951999775, 0.4552237160102125, 0.4599552617985945, 0.4690252430098939, 0.48040828118946455, 0.49095458362373856, 0.49782072437542035, 0.5000000001063967, 0.4983002159372162, 0.494401447354348, 0.4900606178792827, 0.48663489763556866, 0.4848846956219197, 0.48500431552267653, 0.48676182332598195, 0.48964925076788546, 0.4930139555600691, 0.49618187595248864, 0.4985806392388925, 0.4998423471713588, 0.4998466790699849, 0.4986876518187885, 0.4965932524111609, 0.49384423816450496, 0.4907185856135208, 0.4874634278661204, 0.4842855895710807, 0.48135176643242483, 0.4787925688308179, 0.47682868046764715, 0.4790143143573391, 0.4868346734349447, 0.50000000003228, 0.48683467347205034, 0.4790143144678166, 0.4768286811591098, 0.4787925688776139, 0.4813517664087923, 0.48428558957908063, 0.4874634278304283, 0.4907185856598357, 0.49384423815238965, 0.4965932523141716, 0.49868765183897823, 0.4998466791367343, 0.4998423472561506, 0.49858063922618834, 0.49618187593169905, 0.49301395553496363, 0.48964925079317767, 0.4867618233612544, 0.4850043155206073, 0.4848846955959414, 0.4866348975706738, 0.49006061784854515, 0.49440144733215374, 0.49830021591848384, 0.5000000000229394, 0.49782072430287116, 0.4909545835843562, 0.4804082811491477, 0.46902524303397386, 0.4599552618131959, 0.45522371591874444, 0.4555075951963953, 0.4604855812303088, 0.46916057126002564, 0.4800029233215385, 0.49088593547903514, 0.49870809910605485, 0.49831527707243817, 0.47934455649636887, 0.42247898198326606, 0.3556503474739049, 0.2908961490297365, 0.2296371032404957, 0.17325662406164344, 0.12305814306230756, 0.08022427346340089, 0.045779310579230015, 0.020556444624852954, 0.005170900003777734]
        if not gmes:
            gmes = calcul_invariant_gme(n_total, angles, psis)

    if cms:
        plt.plot(angles, cms, label="CM")
    if mus:
        plt.plot(angles, mus, label="MU")
    if gmes:
        plt.plot(angles, gmes, label="GME")
    plt.legend()
    plt.xlabel("ANGLE")

    # Figure name and saving
    figure_name += "_" + partial_name + "_" + str(datetime.today().strftime('%d-%m-%Y-%H-%M'))
    plt.savefig(figure_name)

    # Shows the image
    plt.show()
