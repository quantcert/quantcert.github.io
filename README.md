# Projects related to the quantum certification

This Git repository lists projects related to our work around quantum
certification.

## Method of constructing a block design by Key and Moori

_Paper abstract:_<br/> 
In the field of experimental mathematics, the programs used to obtain various
results rarely follow good software engineering practices, making these results
difficult to reproduce and evaluate. This article presents a formalization of a
method for of combinatorial structures (called block systems) and a validation
of their properties.

More information related to this paper can be found in the
[Designs](Designs) part of this GitHub repository.

## GT IQ 2019's days

The annual meeting of GT IQ is organised in Besançon for 2019, more details can
be found on [GT-IQ'19 (fr)](https://quantcert.github.io/GT-IQ'19)

## Mermin Polynomials for Entanglement Evaluation in Grover’s algorithm and Quantum Fourier Transform

_Paper abstract:_<br/> 
The entanglement of a quantum system can be valuated
using Mermin polynomials. This gives us a means to study entanglement evolution
during the execution of quantum al- gorithms. We first consider Grover’s quantum
search algorithm, noticing that states during the algorithm are maximally
entangled in the direction of a single constant state, which allows us to search
for a single optimal Mermin operator and use it to evaluate entanglement through
the whole execution of Grover’s algorithm. Then the Quantum Fourier Transform is
also studied with Mermin polynomials. A different optimal Mermin operator is
searched at each execution step, since in this case there is no single direction
of evolution. The results for the Quantum Fourier Transform are compared to
results from a previous study of entanglement with Cayley hyperdeterminant. All
our computations can be replayed thanks to a structured and documented
open-source code that we provide.

More information related to this paper can be found in the
[Mermin-eval](Mermin-eval) part of this GitHub repository.

## Entanglement and non-locality of four-qubits connected hypergraph states

_Paper abstract:_<br/> 
We study entanglement and non-locality of connected four-qubit hypergraph
states. One obtains the SLOCC classification from the known LU-orbits. We then
consider Mermin’s polynomials and show that all four-qubit hypergaph states
exhibit non-local behavior. "Finally, we implement some of the corresponding
inequalities on the IBM Quantum Experience.

More information related to this paper can be found in the
[Mermin-hypergraph-states](Mermin-hypergraph-states) part of this GitHub
repository.

## Contextuality degree of quadrics in multi-qubit symplectic polar spaces

_Paper abstract:_<br/> 
Quantum contextuality takes an important place amongst the concepts of quantum
computing that bring an advantage over its classical counterpart. For a large
class of contextuality proofs, aka. observable-based proofs of the
Kochen-Specker Theorem, we formulate the contextuality property as the absence
of solutions to a linear system and define for a contextual configuration its
degree of contextuality. Then we explain why subgeometries of binary symplectic
polar spaces are candidates for contextuality proofs. We report the results of a
software that generates these subgeometries, decides their contextuality and
computes their contextuality degree for some small symplectic polar spaces. We
show that quadrics in the symplectic polar space Wn are contextual for n=3,4,5.
The proofs we consider involve more contexts and observables than the smallest
known proofs. This intermediate size property of those proofs is interesting for
experimental tests, but could also be interesting in quantum game theory.

More information related to this paper can be found in the
[Magma-contextuality](Magma-contextuality) part of this GitHub repository.

## Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank

_Paper abstract:_<br/> 
We study certain physically-relevant subgeometries of binary symplectic polar
spaces of small rank N , W(2N-1,2), when the points of these spaces are
parametrized by canonical N -fold products of Pauli matrices and the associated
identity matrix (i.e., N-qubit observables). Key characteristics of a subspace
of such W(2N-1,2) are: the number of its negative lines, distribution of types
of observables, character of the geometric hyperplane the subspace shares with
the distinguished (non-singular) quadric of W(2N-1,2) and the structure of its
Veldkamp space. W(3,2) features three negative lines of the same type and its
W(1,2)’s are of five different types. W(5,2) is endowed with 90 negative lines
of two types and its W(3,2)’s split into 13 types. 279 out of 480 W(3,2)’s with
three negative lines are composite, i.e. they all originate from the two-qubit
W(3,2) by selecting in the latter a geometric hyperplane and formally adding to
each two-qubit observable, at the same position, the identity matrix if an
observable lies on the hyperplane and the same Pauli matrix for any other
observable. Further, given a W(3,2) and any of its geometric hyperplanes, there
are other three W(3,2)’s possessing the same hyperplane. There is also a
particular type of W(3,2)’s, a representative of which features a point each
line through which is negative. Finally, W(7,2) is found to possess 1908
negative lines of five types and its W(5,2)’s fall into as many as 29 types.
1524 out of 1560 W(5,2)’s with 90 negative lines originate from the three-qubit
W(5,2). Remarkably, the difference in the number of negative lines for any two
distinct types of four-qubit W(5,2)’s is a multiple of four.

More information related to these papers can be found in the
[Magma-contextuality](Magma-contextuality) part of this GitHub repository.

## Testing quantum contextuality of binary symplectic polar spaces on a Noisy Intermediate Scale Quantum Computer

_Paper abstract:_<br/> 
The development of Noisy Intermediate Scale Quantum Computers (NISQC) provides 
for the Quantum Information community new tools to perform quantum experiences 
from an individual laptop. It facilitates interdisciplinary research in the sense 
that theoretical descriptions of properties of quantum physics can be translated 
to experiments easily implementable on a NISCQ. In this note I test large 
state-independent inequalities for quantum contextuality on finite geometric 
structures encoding the commutation relations of the generalized N-qubit Pauli 
group. The bounds predicted by Non-Contextual Hidden Variables theories are 
strongly violated in all conducted experiences.

More information related to this paper can be found in the
[Testing_contextuality](Testing_contextuality) part of this GitHub repository.

## Computational studies of entanglement and quantum contextuality properties towards their formal verification

_Thesis abstract:_<br/> 
Although current quantum computers are limited to the use of a few quantum bits,
the foundations of quantum programing have been growing over the last 20 years.
These foundations have been theorized as early as in the 80’s but the complexity
of their implementation caused these leads to be out of reach until very
recently. In this context, the objective of this thesis is to contribute to the
adaptation of the methods of formal specification and deductive verification of
classical programs to quantum programs. I thus present my contributions to the
study of quantum properties with the end goal of formalizing them. I study in
particular quantum entanglement and quantum contextuality. These properties
allow to classify quantum states and protocols, and in particular to
differentiate them from classical ones. My study of entanglement is based more
specifically on the evolution of entanglement during two quantum algorithms: the
Grover algorithm and the Quantum Fourier Transform. To quantify entanglement
along those algorithms, I use Mermin’s polynomials, which have the advantage of
being implementable in actual quantum computers. My study of contextuality, on
the other hand, relies on finite geometries representing experiments, which are
said to be contextual when no non-contextual classical theory can predict the
results. These geometries are associated with the binary symplectic polar
spaces. We study their structure, and eventually use this structure to get
insights on quantum protocols using contextuality. The study of these properties
led to interesting conjectures which we started to formalize in proof
environments, such as Coq and Why3, but are left as perspective as these works
have not reach a conclusion yet.

More information related to this paper can be found in the
[Computational_studies_of_entanglement_and_contextuality](Computational_studies_of_entanglement_and_contextuality) 
part of this GitHub repository.

## Three-Qubit-Embedded Split Cayley Hexagon is Contextuality Sensitive

_Paper abstract:_<br/> 
It is known that there are two non-equivalent embeddings of the split Cayley
hexagon of order two into W(5, 2), the binary symplectic polar space of rank
three, called classical and skew. Labelling the 63 points of W(5, 2) by the 63
canonical observables of the three-qubit Pauli group subject to the symplectic
polarity induced by the (commutation relations between the elements of the)
group, the two types of embedding are found to be quantum contextuality
sensitive. In particular, we show that the complement of a classically-
embedded hexagon is not contextual, whereas that of a skewly-embedded one is.

More information related to these papers can be found in the
[Magma-contextuality](Magma-contextuality) part of this GitHub repository.

## Multi-qubit doilies: enumeration for all ranks and classification for ranks four and five

_Paper abstract:_<br/> 
For $N \geq 2$, an $N$-qubit doily is a doily living in the $N$-qubit symplectic
polar space. These doilies are related to operator-based proofs of
quantum contextuality. Following and extending the strategy of Saniga et al.
(Mathematics 9 (2021) 2272) that focused exclusively on three-qubit doilies, we
first bring forth several formulas giving the number of both linear and
quadratic doilies for any $N > 2$. Then we present an effective algorithm for
the generation of all $N$-qubit doilies. Using this algorithm for $N=4$ and
$N=5$, we provide a classification of $N$-qubit doilies in terms of types of
observables they feature and number of negative lines they are endowed with. We
also list several distinguished findings about $N$-qubit doilies that are absent
in the three-qubit case, point out a couple of specific features exhibited by
linear doilies and outline some prospective extensions of our approach.

The numerical results related to this paper can be found in the 
[doilies](doilies) part of this GitHub repository.

## Phase sensitivity of entanglement in the Quantum Phase Estimation Algorithm

_Paper abstract:_<br/>
We study entanglement in the pre-QFT part of the Quantum Phase Estimation and the 
Quantum Counting algorithms. In particular we focus on the sensitivity of 
entanglement to the input value (the phase and the ratio of marked elements M/N) 
in some basic cases. One starts from numerical observations and deduce some general 
results in particular regarding the classes of entanglement.

More information related to this paper can be found in the 
[qpea_and_qca](qpea_and_qca) part of this GitHub repository.

## Implementing 2-qubit pseudo-telepathy games on noisy intermediate scale quantum computers

_Paper abstract:_<br/>
It is known that Mermin-Peres like proofs of quantum contextuality can furnish non-
local games with a guaranteed quantum strategy, when classically no such guarantee 
can exist. This phenomenon, also called quantum pseudo-telepathy, has been studied 
in the case of the so-called Mermin Magic square game. In this paper we review in 
detail two different ways of implementing on a quantum computer such a game and 
propose a new Doily game based on the geometry of 2-qubit Pauli group. We show that
the quantumness of these games are almost revealed when we play them on the IBM 
Quantum Experience, however the inherent noise in the available quantum machines 
prevents a full demonstration of the non-classical aspects.   

More information related to this paper can be found in the 
[quantum_game](quantum_game) part of this GitHub repository.

## Contextuality degree of multi-qubit configurations

All information related to this subject can be found in the
[ContextualityDegree](ContextualityDegree) part of this GitHub repository.
