# Exploiting Finite Geometries for Better Quantum Advantages in Mermin-Like Games

Copyright (C) 2024 Colm Kelleher.

Contact: colmkelleher[at]utbm.fr

## Description

Quantum games embody non-intuitive consequences of quantum phenomena, such as entanglement and contextuality. The Mermin-Peres game is a simple example, demonstrating how two players can utilise shared quantum information to win a no - communication game with certainty, where classical players cannot. In this paper we look at the geometric structure behind such classical strategies, and borrow ideas from the geometry of symplectic polar spaces to maximise this quantum advantage. We introduce a new game called the Eloily game with a quantum-classical success gap of $0.2\overline{6}$, larger than that of the Mermin-Peres and doily games. We simulate this game in the IBM Quantum Experience and obtain a success rate of $1$, beating the classical bound of $0.7\overline{3}$ demonstrating the efficiency of the quantum strategy. 
 

The codes are written using python and 
the python library [Qiskit](https://www.qiskit.org/) and are available, 
as well as examples of results, in the following folder 
[Src](https://github.com/quantcert/quantcert.github.io/tree/master/eloily_game/src).

> File "Eloily_2_player_game_brisbane_source_code.py" provides code for running 2-player Eloily game for simulator, noisy simulator, and quantum backend.

> File "Eloily_4_player_game_brisbane_source_code.py" provides code for running 4-player Eloily game for simulator, noisy simulator, and quantum backend.
Both of the above run on the "ibm_brisbane" backend, but that can be changed in-code.

> File "eloily_full_2_player_game_results.csv" gives all results for running the 2-player game.

> File "2_player_grid_doily_results_from_eloily_data.py" gives the results for all Mermin and Doily games, which are subgeometries contained in the eloily.  It gathers these results from the results of the eloily game. It provides the specific grid and doily that give the best respective results.

More details  can be found in the article 
[[K24]](https://arxiv.org/abs/2403.09512).

## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

<a id="K24"/>[K24]  Colm Kelleher, Frédéric Holweck, Péter Lévay **Exploiting Finite Geometries for Better Quantum Advantages in Mermin-Like Games**  [arXiv:2403.09512](https://arxiv.org/abs/2403.09512)
