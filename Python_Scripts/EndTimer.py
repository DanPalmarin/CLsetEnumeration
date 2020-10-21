import datetime
import sys
import os

def time_formatting(end_time):
    if end_time >= 3600:
        end_time = round(end_time/3600, 2)
        print("Total time of program execution: {} hours".format(end_time))
    elif end_time >= 60:
        end_time = round(end_time/60, 2)
        print("Total time of program execution: {} minutes".format(end_time))
    else:
        end_time = round(end_time, 2)
        print("Total time of program execution: {} seconds".format(end_time))


#Store the arguments sent from GAP into variables to use later.
gap_var = sys.argv #Store arguments sent from GAP into the list 'gap_var'
time_file = gap_var[1]
solution_file = gap_var[2]
solutions_file = gap_var[3]
solutions_tallied_file = gap_var[4]

now = datetime.datetime.now()

with open(time_file) as f:
    start_time_string = f.readlines()


start_time = datetime.datetime.strptime(start_time_string[0], "%Y-%m-%d %H:%M:%S.%f") #Read in date string and convert to datetime object
end_time = datetime.datetime.now()

total_time = (end_time - start_time).total_seconds() #Convert the total time to seconds

now = datetime.datetime.now()

print("\nProcess Completed: {}".format(now.strftime("%Y-%m-%d %H:%M"))) #Signal the end of computation
time_formatting(total_time) #Print the total time in a nice format

if os.path.isfile(solutions_tallied_file):
    os.remove(solutions_tallied_file)

if os.path.isfile(solutions_file):
    os.remove(solutions_file)

if os.path.isfile(solution_file):
    os.remove(solution_file)

if os.path.isfile(time_file):
    os.remove(time_file)
