# Contextuality degree of quadrics in multi-qubit symplectic polar spaces <br> _and_ <br> Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank <br> _and_ <br> Three-Qubit-Embedded Split Cayley Hexagon is Contextuality Sensitive

Copyright (C) 2021-2022 Henri de Boutray

Contact: henrideboutray[at]hotmail.fr

## About

The code has been developed to study quantum geometries generated with
symplectic polar spaces in correspondence with the Pauli group. In particular,
contextuality for the article [[dHG+22]](#dHG22) presented as a poster
[[dHG+21]](#dHG21) for QPL'21, subspaces structures for the article
[[SdHG21]](#SdHG21) and the Cayley hexagon for the article [[HdS22]](#HdS22). 
The language chosen was [Magma](http://magma.maths.usyd.edu.au){:target="_blank"} 
since it is a well-established language for mathematics.

These articles are closely related to the content of the program files.

## Installation

We give a heads up to the user: due to the way Magma is achitectured, installation 
is not trivial so some manual work is needed in order to run our code.

Two types of files are used in this project: 
- files containing *intrinsics*, called *packages*
- and files containing "normal" function, called *scripts*

The *intrinsics* are typed and compiled function, their use allow us a greater
flexibility (as overloading), as well as the usual benefits of typing (earlier
error detection, safer usage, ...). But since packages are compiled, they must
be treated differently to scripts. In order to deal with this, the simplest
solution is to put these files in a Magma source code folder. These folders are 
listed in the system environment variable `MAGMA_SYSTEM_SPEC` accessible through 
the Magma prompt with the command `GetEnv("MAGMA_SYSTEM_SPEC");`. But this may 
be impossible for many reasons, in this case please refer to
[this](http://magma.maths.usyd.edu.au/magma/handbook/text/24) documentation
page to understand how to attach a package in Magma. You can obviously contact 
me directly in case you cannot find out
how to attach the packages on your system. All packages are in the 
[src/intrinsics](https://github.com/quantcert/quantcert.github.io/tree/master/Magma-contextuality/src/intrinsics) 
folder.

Once the intrinsics attached, you can call them in scripts or in the Magma
shell. The `Main_***.m` files are such examples of scripts. All scripts are in the 
[src/mains](https://github.com/quantcert/quantcert.github.io/tree/master/Magma-contextuality/src/mains) 
folder.

Assuming Magma is installed on your computer, here is a list of commands you may
use to run our programs. Please keep in mind that they are not the only ones that 
will work and they are only meant to work on Unix-like systems. (You should 
understand commands copied from the web before running them!)

```bash
git clone https://github.com/quantcert/quantcert.github.io.git
cd quantcert.github.io/Magma-contextuality/src
touch ~/.bashrc
echo "export \"MAGMA_USER_SPEC=$(pwd)/intrinsics.spec\"" >> ~/.bashrc
source ~/.bashrc
cd mains
magma Main_AllContextualityChecks.m
```

## Scripts to articles matching

As stated previously, the code follows closely the content of the corresponding 
papers.

The results described in Table 2 of [[dHG+22]](#dHG22) are given by the 
script [Main_AllContextualityChecks.m](src/mains/Main_AllContextualityChecks.m).

[[SdHG21]](#SdHG21) presenting a variety of results, several scripts were 
written to generate them. In particular:
- Sec. 4 is covered by the script [Main_GeometryIntersections_3qubits.m](src/mains/Main_GeometryIntersections_3qubits.m);
- Sec. 5 is covered by the script [Main_Heptads.m](src/mains/Main_Heptads.m);
- Sec. 6 is covered by the script [Main_GeometryIntersections_4qubits.m](src/mains/Main_GeometryIntersections_4qubits.m).

[[HdS22]](#HdS22)'s results (number of embeddings and contextuality of their 
complement) are obtainable by running the script 
[Main_Skewed-Cayley-hexagon.m](src/mains/Main_Skewed-Cayley-hexagon.m).

## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

|                         |                                                    |
|-------------------------|----------------------------------------------------|
|<a id="dHG22"/>[dHG+22]  |Henri de Boutray, Frédéric  Holweck, Alain Giorgetti, Pierre-Alain Masson and Metod Saniga. *Contextuality degree of quadrics in multi-qubit symplectic polar spaces*. [https://doi.org/10.1088/1751-8121/aca36f](https://doi.org/10.1088/1751-8121/aca36f), [arXiv:2105.13798](https://arxiv.org/abs/2105.13798){:target="_blank"}|
|<a id="dHG21"/>[dHG+21]|Henri de Boutray, Frédéric  Holweck, Alain Giorgetti and Pierre-Alain Masson. *Automated detection of contextuality proofs with intermediate numbers of observables*. QPL'21, [Poster](poster-landscape.pdf){:target="_blank"}|
|<a id="SdHG21"/>[SdHG21] |Metod Saniga, Henri de Boutray, Frédéric Holweck and Alain Giorgetti. *Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank*. [doi:10.3390/math9182272](https://doi.org/10.3390/math9182272){:target="_blank"}|
|<a id="HdS22"/>[HdS22]   |Frédéric Holweck, Henri de Boutray and Metod Saniga. *Three-Qubit-Embedded Split Cayley Hexagon is Contextuality Sensitive*. [doi:10.1038/s41598-022-13079-3](https://doi.org/10.1038/s41598-022-13079-3){:target="_blank"}|
