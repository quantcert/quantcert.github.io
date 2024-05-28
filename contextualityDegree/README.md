## New and improved bounds on the contextuality degree of multi-qubit configurations <br> _and_ <br> Revealing contextuality of quantum configurations with a SAT solver

Qontextium
==========

Qontextium is a program which estimates the contextuality degree 
of quantum configurations. It is developped in the C language and is 
using the OpenMP parallelization library.

Copyright (C) 2024 Axel Muller and Alain Giorgetti.

Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France

Contact: axel.muller AT femto-st.fr

Copyright
=========

This program is distributed under the GNU LGPL 2.1 See the enclosed file [LICENSE](./LICENSE).

Installation
============

The program is located in this [folder](https://github.com/quantcert/quantcert.github.io/tree/master/contextualityDegree).
The programs flex and bison need to be installed in your machine for the program to run.
The installation is automatic when running Qontextium for the first time.

Execution
=========


    ./qontextium QUBITS_NUMBER CONFIGURATION [OPTIONS] [--solver SOLVER]

QUBITS_NUMBER                   : the number of qubits in the configuration

CONFIGURATION                   : the configuration to be generated. It can be one of the following:

    --import file_name          : import a configuration from a file
    --elliptics [--complement]  : generate elliptic configurations
    --hyperbolics [--complement]: generate hyperbolic configurations
    --perpsets [--complement]   : generate perpsets configurations
    --subspaces k               : generates a configuration with subspaces of dimension k as contexts
    --affine                    : generate affine configurations

SOLVER

    --solver sat                : use a SAT solver to estimate the contextuality degree
    --solver retrieve           : checks a solution from a solution code

For example to compute the contextuality degree of totally isotrpic subspaces of dimension 1 
(lines) for 2 qubits, you need to run this command :

    ./qontextium 2 --subspaces 1


## References

|                         |                                                    |
|-------------------------|----------------------------------------------------|
|<a id="dHG21"/>[dHG+21]|Muller, Axel and Saniga, Metod and Giorgetti, Alain and de Boutray, Henri and Holweck, Frédéric. *Revealing contextuality of quantum configurations with a SAT solver*. Journées nationales du GDR GPL (Génie de la Programmation et du Logiciel), CNRS, groupe de travail LVP (Langages et Vérification de Programmes), 5-8 juin 2023, Rennes, France. Session posters et démos, 6 juin 2023. [Poster](23poster.pdf){:target="_blank"}|
