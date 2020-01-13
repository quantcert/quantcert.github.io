import csv

import sys, os
sys.path.append(os.getcwd() + "/../")
from mermin_eval.qft import *

# l can go up to 14, r goes from 1 up to 15-l
for l in range(1):
  for r in range(1,2):
    print "[l,r,n]:"
    print [l,r,4]
    v = periodic_state(l,r,4)
    file_name="qft_experience_result/period_"+str(l)+"-"+str(r)+".csv"
    qft_main(v, file_name=file_name, verbose=True)
