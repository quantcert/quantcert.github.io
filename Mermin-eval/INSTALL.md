Installation instructions
=========================

The package is currently developed and maintained only for Linux. It can either 
be executed on the native OS of the user, if [SageMath](http://www.sagemath.org/)
 is installed, or within a Docker container, as detailed below.

Installation with Docker
------------------------

1. If Docker is not installed, follow the instructions at
   https://docs.docker.com/install/ (choose the server installation instructions
   corresponding to your distribution) and 
   https://docs.docker.com/engine/installation/linux/linux-postinstall/ to 
   install Docker and run Docker commands without using sudo.

2. Build the Docker image (all Docker commands are encapsulated in Makefile 
   entries):

    `make docker-build`

    Warnings: The Docker image is a large file, its construction can be very long.
    A network connection is required. If you already have a Docker image with this
    name, either remove it or change the image name in Makefile.

4. At this point, you have several choices, you can either run the scripts that
   generated the results in the article [[BJH+19]](Mermin_eval/#BJH19), or open 
   an interactive session to try out the tools from the module `mermin_eval`.

    * If you want to run the scripts, you can simply use:

      `make grover` or `make qft`.

    * If you want to open an interactive session, use:

      `make interactive`

      This will start a container and a SageMath interactive session.
      See in the `Execution` section how to run the tools.

5. Finally, to quit the container:

    `exit` a first time will close the SageMath interactive session and `exit`
    a second time will quit the container.

`Successfully tested in October 2018 under Linux Ubuntu 16.04, with Docker
18.03.1-ce.`

`Direct Installation`
-------------------

To install Sagemath without going through docker, please follow the instructions
from [SageMath installation page](https://doc.sagemath.org/html/en/installation/).

Execution instructions
======================

Once the installation over, you can enter a SageMath interactive session using 
the command `sage` in your terminal.

Once in the sage session (either in the Docker container or directly on your
computer), you can import the modules, as in the following example. From this 
point, you need some basic knowledges on Python since SageMath is built on 
Python 2.7.

For specific function uses, see the documentation provided 
[in pdf format](doc/build/latex/Mermin-evaluation.pdf) or as a more complete 
[website](doc/build/html/).

Example of usage:

You can run the Grover's algorithm as well as the Mermin evaluation process with
the following commands once in the sage environment:

```python
>>> from mermin_eval.grover import *
>>> v = target_state_ket_string_to_vector("0000")
>>> grover(v)
[0.0963683533242808,1.21302328780365,0.909118107915599,-0.530436080194821]
```