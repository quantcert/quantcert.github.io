## New and improved bounds on the contextuality degree of multi-qubit configurations <br> _and_ <br> Revealing contextuality of quantum configurations with a SAT solver

Qontextium
==========

Qontextium is a program which estimates the contextuality degree 
of quantum configurations. It is developed in the C language and is 
using the OpenMP parallelization library.

Copyright (C) 2024 Axel Muller and Alain Giorgetti.

Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France

Contact: axel.muller AT femto-st.fr

Copyright
=========

This program is distributed under the GNU GPL 2. See the enclosed file [LICENSE](./LICENSE).

Installation
============

Linux is required to run this program.
The program is located in this [folder](https://github.com/quantcert/quantcert.github.io/tree/master/contextualityDegree).
`flex` and `bison` need to be installed in your machine for the program to run.
These tools are commonly available and easy to install on most Linux distributions. 
If you don't have them already, you can typically install them using your distribution's package manager. 
For example, on Debian-based systems: 

    sudo apt-get install flex bison

The compilation is automatic when running Qontextium for the first time.

Execution
=========

Here are the possible options to execute the program:

    ./qontextium QUBITS_NUMBER CONFIGURATION [OPTIONS] [--solver SOLVER]

QUBITS_NUMBER is the number of qubits in the configuration.

CONFIGURATION [OPTIONS] is the configuration to be generated. It can be one of the following:

--import file_name: import a configuration from a file

--elliptics [--complement]: generate elliptic configurations (or their complement or all of them)

--hyperbolics [--complement]: generate hyperbolic configurations (or their complement or all of them)

--perpsets [--complement]: generate perpsets configurations (or their complement or all of them)

--subspaces k: generates a configuration with all subspaces of dimension k as contexts

--affine: generate affine configurations

SOLVER

--solver sat: use a SAT solver to estimate the contextuality degree

--solver retrieve: checks a solution from a solution code

For example, to compute the contextuality degree of totally isotropic subspaces of dimension 1 
(lines) for 2 qubits, run this command:

    ./qontextium 2 --subspaces 1

#### Graph Visualization

There are multiple filters you can use to visualize the contextual graphs :

NOTHING
- This filter always returns false, meaning no lines will pass through this filter. It effectively filters out all lines.


ALL
- This filter always returns true, meaning all lines will pass through this filter.


POINT DEGREE
- This filter checks if any of the points in the geometry (for the given line i) match a specific point degree specified by the user.


OBSERVABLE VALUE
- This filter checks if any of the points in the geometry for the line i match a specific observable value given by the user


LINE DEGREE
- This filter checks if the sorted degrees of the points in the geometry for each line exactly match the degrees specified by the user


SYMMETRIC
- This filter checks if all points in the geometry for line i are symmetric. An observable is symmetric if the number of Y's is even.

## References

|                         |                                                    |
|-------------------------|----------------------------------------------------|
|<a id="MSGDH23"/>[MSGDH23]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and de Boutray, Henri and Holweck, Frédéric. *Revealing contextuality of quantum configurations with a SAT solver*. Journées nationales du GDR GPL (Génie de la Programmation et du Logiciel), CNRS, groupe de travail LVP (Langages et Vérification de Programmes), 5-8 juin 2023, Rennes, France. Session posters et démos, 6 juin 2023. [Poster](23poster.pdf){:target="_blank"}|
|<a id="MSGDH24"/>[MSGDH24]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and de Boutray, Henri and Holweck, Frédéric. *New and improved bounds on the contextuality degree of multi-qubit configurations*. Mathematical Structures in Computer Science, 2024. [Article](https://doi.org/10.1017/S0960129524000057){:target="_blank"}|
