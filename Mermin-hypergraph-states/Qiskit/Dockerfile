FROM sagemath/sagemath:9.0

RUN sudo apt-get update && \
 sudo apt-get install -y make

RUN sage --pip install qiskit

COPY .bashrc .bashrc
