TMP-SRC=mermin_eval/*.pyc examples/*.sage.py

# Build a Docker machine with Sage
docker-build:
	docker build --tag sagemath:make .

docker-image-clean: docker-container-clean
	docker image rm sagemath:make
	docker image rm sagemath/sagemath:8.6

docker-container-clean:
	docker ps -aq --filter=ancestor=sagemath:make | xargs -I {} docker rm {}

grover:
	sh examples/run-scripts/experiences.sh grover

qft:
	sh examples/run-scripts/experiences.sh qft

# Launch the container for an interactive session
interactive:
	sh examples/run-scripts/interactive.sh

clean-src:
	rm -f $(TMP-SRC)

doc-latex:
	sage -sh -c "$(MAKE) -C doc latex"
	make -C doc/build/latex

doc-latex-clean:
	make -C doc/build/latex clean

doc-html:
	sage -sh -c "make -C doc html"

doc-clean:
	sage -sh -c "make -C doc clean"
