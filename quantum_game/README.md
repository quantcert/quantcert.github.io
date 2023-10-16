# Implementing 2-qubit Pseudo-Telepathy Games on Noisy Intermediate Scale Quantum Computers

Copyright (C) 2023 Colm Kelleher.

Contact: colmkelleher[at]utbm.fr

## Description

We implemented the well-known Mermin game on a 7- and 16-qubit noisy intermediate scale quantum computer, via IBM's Quantum Experience (https://quantum-computing.ibm.com/). The objective was to test non-contextual hidden variable models via the success rate of Alice and Bob playing the game. 
In addition, variations of the game were played to reduce the effect of noise on the system. Variations include changing the game from the Mermin to the Doily game via the geometrical structure of the classical strategy, and changing the game from a line-line to a point-line scenario.
Two different implementation methods were used for both doily and Mermin games. The unitary method utilises basis transformations to provide outputs associated with sequential measurements on a context. The delegation method utilises auxiliary qubits for repeated measurements of the same shared state.
Both noiseless and noisy simulators were used, and the Lagos machine.
 

The codes are written using python and 
the python library [Qiskit](https://www.qiskit.org/) and are available, 
as well as examples of results, in the following folder 
[Src](https://github.com/quantcert/quantcert.github.io/tree/master/quantum_game/src).

More details  can be found in the article 
[[K23]](https://arxiv.org/abs/2310.07441).

## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

<a id="K23"/>[K23]  Colm Kelleher, Mohammad Roomy, Frédéric Holweck **Implementing 2-qubit Pseudo-Telepathy Games on Noisy Intermediate Scale Quantum Computers**  [arXiv:2310.07441](https://arxiv.org/abs/2310.07441)
