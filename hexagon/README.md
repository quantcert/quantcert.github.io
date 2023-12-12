# Testing Quantum Contextuality of the Skew Hexagon Complement and 3-qubit Elliptic Quadric on a Noisy Intermediate Scale Quantum Computer

Copyright (C) 2023 Colm Kelleher.

Contact: colm.kelleher[at]utbm.fr

## Description

We tested the Cabello contextuality inequalities of the following two geometries: the complement of the Skew split Cayley Hexagon and the Elliptic Quadric both contained in symplectic polar space W(5,2). For this we used the online quantum computers provided by [IBM Quantum Experience](https://quantum-computing.ibm.com/).

We generated the set of all 3-qubit operators, and their commutation relations (W(5,2)). From this we constructed some subgeometries of commutation relations, including the classically embedded split Cayley Hexagon, the skew embedded Split Cayley Hexagon, both of their complements in W(5,2), and the elliptic quadric GQ(2,4). These geometries are constructed as lines (lists) of mutually commuting 3 qubit operators (stored as strings, e.g. 'XIY') along with the parity condition on the lines (+ or - I). For example, [['ZZI', 'YYI', 'XXI'], -1] indicates the line formed from the mutually commuting operators ZZI, YYI and XXI, whose product is -I, the identity. 

We then tested the Cabello inequality (https://journals.aps.org/pra/abstract/10.1103/PhysRevA.82.032110) for two of the geometries that are known to be contextual - the Skew hexagon complement and the elliptic quadric. 

This involved implementing a circuit on a quantum computer using the qiskit library, making measurements of each operator along each context, and summing the expectation measurement values along positive lines minus those for negative lines. Cabello shows that in a non-contextual hidden variable model NCHV, this is bounded by the number of overall lines minus twice the degree of contextuality. We show that the measured value for this expression is higher than this bound, indicating no NCHV models.
 

The codes are written using python and 
the python library [Qiskit](https://www.qiskit.org/) and are available, 
as well as examples of results, in the following folder 
[Src](https://github.com/quantcert/quantcert.github.io/tree/master/hexagon/Src).

## Copyright




