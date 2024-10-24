# Qontextium

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

--import assignment [FILE]: imports a configuration from a file (see ./misc/qa_grid.txt for an example) and estimates its contextuality degree

--import hypergram [FILE1] [FILE2]: imports a hypergram from two files, checks assignability. When assignable, estimates the contextuality degree. The first file describes the hypergraph. The second file describes a Gram matrix on the same vertices (see ./misc/grid_hypergraph.txt and ./misc/grid_gram.txt for an example)

--elliptics [--complement]: generates elliptic configurations (or their complement or all of them)

--hyperbolics [--complement]: generates hyperbolic configurations (or their complement or all of them)

--perpsets [--complement]: generates perpsets configurations (or their complement or all of them)

--subspaces k: generates a configuration with all subspaces of dimension k as contexts

--affine: generates affine configurations

SOLVER

--solver sat: uses a SAT solver to estimate the contextuality degree (eventually gives the optimum)

--solver heuristic: uses the heuristic method presented in [MSGHK24](#MSGHK24) to estimate the contextuality degree (faster, but no guarantee of finding the optimum)

--solver retrieve: checks a solution from a solution code

For example, to compute the contextuality degree of totally isotropic subspaces of dimension 1 (lines) for 2 qubits, run this command:

    ./qontextium 2 --subspaces 1

Example files are provided in the `misc` folder, for example you can import a Peres-Mermin magic square with this command:

    ./qontextium 2 --import assignment ./misc/qa_grid.txt

#### Graph Visualization

There are multiple filters you can use to visualize the contextual graphs :

NOTHING
- This filter always returns false, meaning no lines will pass through this filter. It effectively filters out all lines.


ALL
- This filter always returns true, meaning all lines will pass through this filter.


POINT DEGREE
- This filter checks if any of the points in the geometry (for a given line) match a specific point degree specified by the user.


OBSERVABLE VALUE
- This filter checks if any of the points in the geometry for a line match a specific observable value given by the user


LINE DEGREE
- This filter checks if the sorted degrees of the points in the geometry for each line exactly match the degrees specified by the user


SYMMETRIC
- This filter checks if all points in the geometry for a line are symmetric. An observable is symmetric if the number of Y's is even.

## References

|                         |                                                    |
|-------------------------|----------------------------------------------------|
|<a id="MSGDH23"/>[MSGDH23]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and de Boutray, Henri and Holweck, Frédéric. *Revealing contextuality of quantum configurations with a SAT solver*. Journées nationales du GDR GPL (Génie de la Programmation et du Logiciel), CNRS, groupe de travail LVP (Langages et Vérification de Programmes), 5-8 juin 2023, Rennes, France. Session posters et démos, 6 juin 2023. [Poster](23poster.pdf){:target="_blank"}|
|<a id="MSGDH24"/>[MSGDH24]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and de Boutray, Henri and Holweck, Frédéric. *New and improved bounds on the contextuality degree of multi-qubit configurations*. Mathematical Structures in Computer Science, 2024. [Article](https://doi.org/10.1017/S0960129524000057){:target="_blank"}|
|<a id="MSGHK24"/>[MSGHK24]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and Holweck, Frédéric and Kelleher, Colm. *A new heuristic approach for contextuality degree estimates and its four- to six-qubit portrayals*. [Article](https://doi.org/10.48550/arXiv.2407.02928){:target="_blank"}|
|<a id="MG24"/>[MG24]|Muller, Axel and Giorgetti, Alain. *An abstract structure determines the contextuality degree of observable-based Kochen-Specker proofs*. [Article](https://arxiv.org/html/2410.14463v1){:target="_blank"}|
