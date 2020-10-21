import sys
import os
import datetime
from gurobipy import *
from contextlib import contextmanager

##### USER DEFINED FUNCTIONS #####
@contextmanager #This surpresses output to the console
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

#Store the arguments sent from GAP into variables to use later.
gap_var = sys.argv #Store arguments sent from GAP into the list 'gap_var'
lp_file = gap_var[1]
group_name = gap_var[2]
solution_count = int(gap_var[3])
variable_sharing_directory = gap_var[4]
gurobi_log_directory = gap_var[5]
c = int(gap_var[6])

#Solution file
solution_file =  variable_sharing_directory + "Solution.txt"

#Analyze model
print("Gurobi is looking for a solution...") #optional - it's nice to see the current process (Gurobi or forbidden restrictions) for long computations.
with suppress_stdout():
    model = read(lp_file)
    model.Params.LogFile = "{0}{1}_{2}".format(gurobi_log_directory, group_name, "gurobi.log")
    model.Params.NodeLimit = 1000000 #Assume 0 solution if optimization has explored 1 million nodes with no solution found
    #model.Params.TimeLimit = 1
    
   # model.Params.TimeLimit = 18000 #Assume 0 solution if optimization takes 5 hours
    # if c in [2,3]:
        # model.Params.TimeLimit = 1
    # elif c == 4 and solution_count-1 == 2:
        # model.Params.TimeLimit = 1
    # elif c == 5 and solution_count-1 == 3:
        # model.Params.TimeLimit = 1
    
    #model.Params.NodefileStart = 0.5 #When the amount of memory used to store nodes (measured in GBytes) exceeds the specified parameter value, nodes are compressed and written to disk. Performance penalty < 10% (usually).
    #model.Params.Threads = 4 #Reducing the number of threads will reduce the amount of memory used (significantly decreases performance)
    #model.Params.NumericFocus = 3 #0-3; Try this if you get "warning: very big kappa". Higher is more accurate but slower.
    
    model.optimize()

optimal_solution = model.SolCount

if optimal_solution == 0:
    print("No solution.\n")
    with open(solution_file, "w") as f:
        f.write("solution := 0;") #This puts the solution into a format that GAP will be able to read and store easily.
else:
    #Find optimal solution from models and track at what date/time the solution was found
    gurobi_solution = [j+1 for j, x in enumerate(model.getAttr("x")) if x > 0.5]
    now = datetime.datetime.now()

    #Print the most recently found solution
    print("Solution {0} (length {3}) found at {2}:\n{1}\n".format(solution_count, gurobi_solution, now.strftime("%Y-%m-%d %H:%M"), len(gurobi_solution)))

    #Write solution to text file for GAP to read
    with open(solution_file, "w") as f:
        f.write("solution := {};".format(gurobi_solution)) #This puts the solution into a format that GAP will be able to read and store easily.


