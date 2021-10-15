# Automated detection of contextuality proofs with intermediate numbers of observables _and_ Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank

Copyright (C) 2021 Henri de Boutray

Contact: henri.de_boutray[at]univ-fcomte.fr

## About

The code has been developed to study quantum geometries generated with
symplectic polar spaces in correspondence with the Pauli group. In particular,
contextuality for the article [[dJG+21]](#dJG21) presented as a poster
[[dHG+21']](#dHG21b) for QPL'21 and subspaces structures for the article
[[SdHG21]](#SdHG21). The language chosen was
[Magma](http://magma.maths.usyd.edu.au) since it is a well-established language 
for theoretical mathematics.

Both articles are closely related to the content of the program files.

## Installation

Unfortunately, Magma is a proprietary software, so unlike our previous projects,
we are not able to release a docker image containing all the elements necessary
to run the code out of the box. Some manual work will be necessary.

Two types of files are used in this project: 
- files containing *intrinsics*, called *packages*
- and files containing "normal" function, called *scripts*

The *intrinsics* are typed and compiled function, their use allow us a greater
flexibility (as overloading), as well as the usual benefits of typing (earlier
error detection, safer usage, ...). But since packages are compiled, they must
be treaded differently to scripts. In order to deal with this, the simplest
solution is to put these files in a Magma source code folder. These folders are 
listed in the system environment variable `MAGMA_SYSTEM_SPEC` accessible through 
the Magma prompt with the command `GetEnv("MAGMA_SYSTEM_SPEC");`. But this may 
be impossible for many reasons, in this case please refer to
[this](https://magma.maths.usyd.edu.au/magma/handbook/text/24) documentation
page to understand how to attach a package in Magma. Since this project was my
first use of Magma, I had the occasion of tinkering quite a bit with its various
components, so you can obviously contact me directly in case you cannot find out
how to attach the packages on your system. All packages are in the 
[src/intrinsics](https://github.com/quantcert/quantcert.github.io/tree/master/Magma-contextuality/src/intrinsics) 
folder.

Once the intrinsics attached, you can call them in scripts or in the Magma
shell. The `Main_***.m` are such examples of scripts. All scripts are in the 
[src/mains](https://github.com/quantcert/quantcert.github.io/tree/master/Magma-contextuality/src/mains) 
folder.

## Scripts to articles matching

As stated previously, the code follows closely the content of the corresponding 
papers. 

The main results of [[dJG+21]](#dJG21), descibed in Table 2 are given by the 
script [Main_AllContextualityChecks-V2.m](src/mains/Main_AllContextualityChecks-V2.m).

[[SdHG21]](#SdHG21) presenting a variety of results, several scripts were 
written to generate them. In particular:
- Sec. 4 is covered by the script [Main_GeometryIntersections_3qubits.m](src/mains/Main_GeometryIntersections_3qubits.m);
- Sec. 5 is covered by the script [Main_Heptads.m](src/mains/Main_Heptads.m);
- Sec. 6 is covered by the script [Main_GeometryIntersections_4qubits.m](src/mains/Main_GeometryIntersections_4qubits.m).

## Copyright

This program is distributed under the GNU GPL 3. See the enclosed file 
[LICENSE](LICENSE).

## References

|                         |                                                    |
|-------------------------|----------------------------------------------------|
|<a id="dHG21"/>[dHG+21]  |Henri de Boutray, Frédéric  Holweck, Alain Giorgetti and Pierre-Alain Masson. *Automated synthesis of contextuality proofs from subspaces of symplectic polar spaces*. [arXiv:2105.13798](https://arxiv.org/abs/2105.13798)|
|<a id="dHG21b"/>[dHG+21']|Henri de Boutray, Frédéric  Holweck, Alain Giorgetti and Pierre-Alain Masson. *Automated detection of contextuality proofs with intermediate numbers of observables*. QPL'21, [Poster](poster-landscape.pdf)|
|<a id="SdHG21"/>[SdHG21] |Metod Saniga, Henri de Boutray, Frédéric Holweck and Alain Giorgetti. *Taxonomy of Polar Subspaces of Multi-Qubit Symplectic Polar Spaces of Small Rank*. [doi:10.3390/math9182272](https://doi.org/10.3390/math9182272)|
