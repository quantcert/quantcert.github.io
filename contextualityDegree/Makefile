compile:
	make -C lib/BCpackage-0.40 bc2cnf > /dev/null 
	cd ./lib/kissat_gb && \
	if [ ! -d ./build ]; then \
		./configure; \
	fi
	make -C ./lib/kissat_gb > /dev/null

clean: 
	make -C ./lib/kissat_gb/ clean
	make -C ./lib/BCpackage-0.40/ clean 
