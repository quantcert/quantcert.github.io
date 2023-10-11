# Implementing 2-qubit pseudo-telepathy games on noisy intermediate scale quantum computers

Copyright (C) 2023 Colm Kelleher

Contact: colm.kelleher[at]utbm.fr

## Description

We contextual quantum games with the online quantum computers of the [IBM Quantum Experience](https://quantum-computing.ibm.com/). More precisely
we provide codes to generate the points (multiqubit Pauli matrices) and lines (triple of mutually commuting operators whose product is +I, -I) of 
some specific finite geometries and using these geometries implemented some Mermin and Mermin-like games on a noisy intermediate-scale quantum computer.
We tested out two different geometries - the Mermin grid and Doily - of symplectic polar space of rank 2. For each game we played it using two methods - a unitary basis transformation method and a delegation qubit method. The latter allowed for a secondary scenario, the point-line scenario, for the mermin and doily games. For each scenario and game, the classical bounds allowed by non-contextual hidden variable models were not beaten, although were beaten using the noisy and noiseless simulators.
Codes to play the games are provided.
 

The codes are written using python and 
the python library [Qiskit](https://www.qiskit.org/) and are available, 
as well as examples of results, in the following folder 
[Src](https://github.com/quantcert/quantcert.github.io/tree/master/quantum_game/src).



## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

<a id="K23"/>[K23]  Colm Kelleher **Implementing 2-qubit pseudo-telepathy games on noisy intermediate scale quantum computers**
