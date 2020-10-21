import datetime
import sys

#Store the arguments sent from GAP into variables to use later.
gap_var = sys.argv #Store arguments sent from GAP into the list 'gap_var'
time_file = gap_var[1]

now = datetime.datetime.now()

with open(time_file, "w") as f:
    f.write(str(now))
