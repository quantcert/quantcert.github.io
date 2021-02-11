# Testing Quantum Contextuality of the Symplectic Polar Space on a Noisy Intermediate Scale Quantum Computer

Copyright (C) 2021 Frédéric Holweck.

Contact: frederic[at]utbm.fr

## Description

We tested macroscopic contextual inequalities with the online quantum computers of the [IBM Quantum Experience](https://quantum-computing.ibm.com/). More precisely
we provide codes to generate the points (multiqubit Pauli matrices) and lines (triple of mutually commuting operators whose product is +I, -I) of 
some specific finite geometries and we evaluate for these geometries some contextual inequalities studied by [Adan Cabello](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.82.032110).
The computations show that for the symplectic polar space of rank 2 and order 2 - a point-line configuration known as the doily - and for the symplectic polar space of rank 3 and order 2, the inequalities are strongly violated by the results provided by the IBM Quantum Computers. For the symplectic polar space of rank 3 and order 2, our experiment involves the measurement of 315 different 3-qubits contexts. 
Codes to test the inequalities on hyperbolic quadrics (subgeometries corresponding to the Mermin-Pere squares in the rank 2 case) are also provided.
 

The codes are written using python and 
the python library [Qiskit](https://www.qiskit.org/) and are available, as well as examples of results, in the following folder [Sources and Results](https://github.com/quantcert/quantcert.github.io/tree/master/Testing_contextuality/Src).

More details  can be found in the article 
[[H21]](https://arxiv.org/abs/2101.03812).

## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

<a id="H21"/>[H21]  Frédéric Holweck **Testing Quantum Contextuality of the Symplectic Polar Space on a Noisy Intermediate Scale Quantum Computer**  [arXiv](https://arxiv.org/abs/2101.03812)
