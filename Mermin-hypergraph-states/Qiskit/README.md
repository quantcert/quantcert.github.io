# Mermin on Qiskit

*With participation of Henri de Boutray*

This project aims to use the Qiskit library to study entanglement in various
cases. It contains files to easily run quantum algorithms (a
light overhead has been added to the default Qiskit method used to submit jobs,
in order to simplify this process), and evaluate states resulting of these
algorithms using optimized Mermin operators.

## Content of this project

The pipeline of the evaluation process is as follows:
1. Compute the optimized Mermin operator. The result will be stored as a list of
  parameters defining the monomials of the Mermin polynomial:
  `((a_1,a_2,...,a_n),(a'_1,a'_2,...,a'_n))` where each `a_i` or `a'_i` is a
  triplet of reals (see the article for more details).
2. Run the algorithm (either using the default Qiskit method, or better, by
  feeding the algorithm to the function that will also perform the evaluation
  and that wrap the Qiskit methods)
3. Evaluate the output state.

Two examples are given here on how to use this module:
1. In the QFT, each state of the execution is studied, but since we are in a
  quantum context, one cannot directly study the states from a single run, so
  the algorithm is split in a series of runs, for which each one contains only
  the first k gates, for k between 0 and the number of gates. In this process,
  the optimization was performed externally, by the `mermin_eval` module
  previously developed for sage. This module is not included in this project,
  but one can fetch it and use it externally, or ask the authors for more
  information.
2. For the graphstates, only the entanglement value of a graphstate was studied,
  so only one evaluation was needed. The optimization process is included in the
  module, in the submodule `graphstates_opti`.

## Documentation

The code documentation can be found
[in pdf format](app/doc/build/latex/Mermin-hypergraphs.pdf) or 
[as a website](app/doc/build/html). It is generated from the code using 
[Sphinx](http://www.sphinx-doc.org).

## Requirements

This module requires:
- Python 3.7+
- Qiskit 0.20.0

In order to build the documentation, the software Sphinx (version 1.8.5) has
been used.

## Run it

### Docker

A *Dockerfile* has been produced in order to easily reproduce the results
associated with this project. In order to use it, one must have Docker
installed, and run the command `docker build --tag sagemath:qiskit .` from the
folder containing the *Dockerfile*.

Once the image built, one can run if using the command:
```bash
docker run -it --mount type=bind,source=$(pwd)/app,target=/home/sage/app \
    sagemath:qiskit "sh -c 'cd /home/sage/app; bash'"
```
(Change the `$(pwd)/app` with a path compatible with you OS, for example on
Windows it would be `%cd%/app`)

with this command, your environment is all set to run everything you want to,
and the files are dynamically mounted so you can edit them on your favorite
editor, and run them from the container.

Since the process of running the various files is made transparent, the
rest of the description on how to run the project will be provided in the
following section.

### On your machine directly

If you installed the requirements, you can also skip Docker and run the various
files directly on your computer.

We encourage you to explore the various files provided. And for an example on
how to use them, you can go to the folder `examples` and run
```bash
python QFT-4-1-1.py
```
or
```bash
python graphstates.py
```

To rebuild the documentation, just type `make doc-html` or `make doc-latexpdf`
in the `app` folder, the result will be in `doc/build`. 