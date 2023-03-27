# source activate qiskit_env

build:
	docker build --tag sagemath:qiskit .

interactive:
	@docker run -it --mount type=bind,source=$(shell pwd)/app,target=/home/sage/app \
		sagemath:qiskit "sh -c 'cd /home/sage/app; bash'"