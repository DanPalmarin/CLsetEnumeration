import sys


def delete_end(lp_file, lp_length, end_length): #This removes the "ending" part of the GAP output file
    new_file_length = lp_length - end_length
    with open(lp_file, "r") as f:
        lines = f.readlines()

    with open(lp_file, 'w') as f:
        f.writelines(lines[:new_file_length])


def file_len(lp_file): #This determines how many lines are in the .lp text file (GAP output file/Gurobi input file)
    with open(lp_file) as f:
        for w, l in enumerate(f):
            pass
    return w + 1


#Store the arguments sent from GAP into variables to use later.
gap_var = sys.argv #Store arguments sent from GAP into the list 'gap_var'
starting_LP = int(gap_var[1])
ending_LP = int(gap_var[2])
gurobi_output_directory = gap_var[3]
group_name = gap_var[4]
lp_file_end = int(gap_var[5])
c = int(gap_var[6])

#Store lp files
lp_files = []
for i in range(starting_LP, (ending_LP+1)):
    file_ending = "_c_{0}.lp".format(str(i))
    temp_file = "{0}{1}{2}".format(gurobi_output_directory, group_name, file_ending)
    lp_files.append(temp_file)


#We delete the ending of each lp file that is and will be used
for d in range(c, ending_LP+1):
    lp_file_length = file_len(lp_files[d-1])
    delete_end(lp_files[d-1], lp_file_length, lp_file_end)
