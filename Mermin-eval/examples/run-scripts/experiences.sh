_uid="$(id -u)"
_group="$(id -g)"

cd ..;
docker run \
  -v "$(pwd)"/Mermin-eval:/home/sage/Mermin-eval \
  sagemath:make \
    "sudo chown -R sage:sage /home/sage/Mermin-eval; \
    cd Mermin-eval/examples; \
    sage $1_experiences.sage; \
    cd ../..; \
    sudo chown -R ${_uid}:${_group} /home/sage/Mermin-eval"
