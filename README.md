# Projects related to the quantum certification

This Git repository lists projects related to our work around quantum
certification.

## Method of constructing a block design by Key and Moori

Key and Moori present a method of building a block designs completed with a part
on graphs that we are not interested in for our work. The linked source files
are present in the [Designs](Designs) part of this GitHub repository.

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

## Automated detection of contextuality proofs with intermediate numbers of observables

_Paper abstract:_<br/> 
Quantum contextuality takes an important place amongst the concepts of quantum
computing that bring an advantage over its classical counterpart. For a large
class of contextuality proofs, aka. observable-based proofs of the
Kochen-Specker Theorem, we first formulate the contextuality property as the
absence of solutions to a linear system. Then we explain why subgeometries of
binary symplectic polar spaces are candidates for contextuality proofs. We
report first results of a software that generates these subgeometries and
decides their contextuality. The proofs we consider involve more contexts and
observables than the smallest known proofs. This intermediate size property of
those proofs is interesting for experimental tests, but could also be
interesting in quantum game theory.

More information related to this paper can be found in the
[Magma-contextuality](Magma-contextuality) part of this GitHub repository.

## Automated detection of contextuality proofs with intermediate numbers of observables _and_ Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank

This projects resulted in two papers, described by the two following abstracts.

_Paper abstract:_<br/> 
Quantum contextuality takes an important place amongst the
concepts of quantum computing that bring an advantage over its classical
counterpart. For a large class of contextuality proofs, aka. observable-based
proofs of the Kochen-Specker Theorem, we first formulate the contextuality
property as the absence of solutions to a linear system. Then we explain why
subgeometries of binary symplectic polar spaces are candidates for contextuality
proofs. We report first results of a software that generates these subgeometries
and decides their contextuality. The proofs we consider involve more contexts
and observables than the smallest known proofs. This intermediate size property
of those proofs is interesting for experimental tests, but could also be
interesting in quantum game theory.

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
