**Method of building block designs**

In experimental mathematics, the programs used to achieve various results rarely follow good software engineering practices, making these results difficult to replicate and evaluate. This section presents a formalization of a method for constructing combinatorial structures (called "block designs") and a validation of their properties.

***Context***

This work is part of research into finite geometries called quantum geometries, because they are linked to _quantum_ contextuality. Planat et al. have shown how to construct these geometries, which are special cases of _block designs_, from groups of permutations, but without publishing a program for this construction.

A literature search on the origin of this method led to an earlier reference, which presents a simpler method of constructing block systems from primitive permutation groups. This article is completed by a program, but the latter does not contain a code to validate this method. This method was formalised and then validated by enumeration. More precisely, block systems built according to this method are validated if they have all the characteristics announced in a proposal of the Key and Moori article. We use the Magma environment, which consists of a structured imperative language and an extensive library of mathematical functions, particularly group theory and _designs_. This code is present in the [./Designs/src](https://github.com/quantcert/quantcert.github.io/tree/master/Designs/src) part of this GitHub directory.

A paper has been written presenting good practices to be taken when developing a program to corroborate proposals and apply them to a block designs construction method.
