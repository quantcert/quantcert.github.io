import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_eval.run_circuit import *

circuit = [[('I',2),('H',1)],[('H',1),('H',1),('I',1)],[('I',1),('S',2)]]
print "Circuit easly readable in a terminal"
print print_circuit(circuit, to_latex=False)
print "Circuit to copy paste in a LaTeX document using the qcircuit package"
"""
The content of the LaTeX document using the qcircuit package should be something 
like:

"\\documentclass[preview]{standalone}
\\usepackage{qcircuit}
\\usepackage{amsmath}
\\usepackage{physics}
\\begin{document}
\\begin{align*}
  \\Qcircuit @C=1em @R=.7em {
    & \\qw & \\gate{H}& \\qw & \\qw \\\\
    & \\qw & \\gate{H}& \\multigate{1}{S} & \\qw \\\\
    & \\gate{H}& \\qw & \\ghost{S} & \\qw \\\\
  }
\\end{align*}
.
\\end{document}"
"""
print print_circuit(circuit, to_latex=True)