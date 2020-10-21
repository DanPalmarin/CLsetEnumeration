import sys
import os
import fileinput
import shutil
from contextlib import contextmanager
from collections import Counter

#Store the arguments sent from GAP into variables to use later.
gap_var = sys.argv #Store arguments sent from GAP into the list 'gap_var'
n = int(gap_var[1])
group_name = gap_var[2]
order = int(gap_var[3])
solution_file = gap_var[4]
solutions_summary_directory = gap_var[5]
temp_sum_file = gap_var[6]

final_file = solutions_summary_directory + group_name + "_Solutions_Summary.txt"

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

def solutions_summary(degree, group_order, grp_name, sol_file, temp_file):
    length_summary = []
    orbit_summary = []
    CL_sets = 0
    
    #This reads in the length (size of solution) and size of the family from the orbit_size.txt.
    #It stores these values in two lists (length_summary and orbit_summary).
    with open(sol_file) as f:
        for line in f:
            if line.find("Solution size") != -1:
                solution_length = [int(i) for i in line.split() if i.isdigit()]
                length_summary.append(solution_length[0])
            elif line.find("Size of solution family") != -1:
                orbit_size = [int(i) for i in line.split() if i.isdigit()]
                orbit_summary.append(orbit_size[0])
            elif line.find("Cameron-Liebler set: true") != -1:
                CL_sets += 1

    #Creates a list of size two tuples. For each tuple, the first digit is the length and the second digit is the orbit size.
    combined_with_repeats = list(zip(length_summary, orbit_summary))
    combined_no_repeats = sorted(set(combined_with_repeats)) #Removes all repeated tuples from the above list of tuples.

    lines_output = []
    total = 0
    
    for i in combined_no_repeats:
        amount = combined_with_repeats.count(i) #Counts how many repeated tuples there are. This gives us how many solutions there are of a certain length and orbit size.
        output_text = "Size {1}: {0} solution(s) with family size {2}.\n".format(amount, i[0], i[1])
        lines_output.append(output_text)
        total += amount

    length_summary.sort()
    length_dict = Counter(length_summary)
    length_total_output = []
    
    for i in length_dict:
        length_output_text = "Total of {0} size {1} solution(s).\n".format(length_dict[i], i)
        length_total_output.append(length_output_text)

    totals_text = "Total number of solutions: {0}\n".format(total) #Tracks the total number of solutions
    CL_totals = "Total number of solutions that are Cameron-Liebler sets: {0}".format(CL_sets)
    
    with open(temp_file, 'w') as f:
        f.writelines(totals_text)
        f.writelines(CL_totals)
        f.write("\n\n")
        f.writelines(length_total_output)
        f.write("\n")
        f.writelines(lines_output)
        f.write("\n\n")

def merge_files(file1, file2, new_file):
    with open(new_file, 'wb') as wfd:
        for f in [file2, file1]:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd) #Add 'wfd.write(b"\n")' after this line if a space is required (shouldn't be necessary)

######################
    
#Prepare the solution and orbit summary text file
solutions_summary(n, order, group_name, solution_file, temp_sum_file)

#Merge the two files (list of solutions and tallied solutions)
merge_files(solution_file, temp_sum_file, final_file)
