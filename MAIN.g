#Reference the 'Cameron-Lieber Set Enumeration Manual' before continuing.
#Note: Open GAP from powershell/terminal using the '-o' command line option with your total RAM; i.e. 'gap -o 64g'
#GAP terminal input to run this file (Linux OS): Read("/home/daniel/Documents/Graduate_Research/Github/CLsetEnumeration/MAIN.g");
#GAP powershell/command prompt to run this file (Windows OS): Read("C:/Users/Daniel/Documents/Graduate_Research/Github/CLsetEnumeration/MAIN.g");


######## TABLE OF CONTENTS #######

# SECTION 0 - PREAMBLE
# SECTION 1 - USER DEFINED FUNCTIONS
# SECTION 2 - PRELIMINARY CHECKS
# SECTION 3 - DETERMINING THE MAXIMUM CLIQUES
# SECTION 4 - FINDING THE MAXIMUM COCLIQUES VIA LINEAR PROGRAMMING
# SECTION 5 - TEX OUTPUT: ANALYSES OF GROUP AND MAXIMUM COCLIQUES
# SECTION 6 - FINDING ALL CAMERON-LIEBLER SETS
# SECTION 7 - TEX OUTPUT: REPORT CAMERON-LIEBLER SETS

##################################


##### SECTION 0 - PREAMBLE #####

#Home directory (location of MAIN.g and required folders) and location of python (used to call Gurobi).
home_dir := "C:/Users/Daniel/Documents/Graduate_Research/Github/CLsetEnumeration/";
python3 := "C:/Users/Daniel/Python310/python.exe";

### USER INPUT (REQUIRED) ###

#Option 1 (suggested): Reference the group using its index in the PrimitiveGroups library or the TransitiveGroups library.
grp := PrimitiveGroup(4,1);
group_name := Name(grp);
RemoveCharacters(group_name, " "); #Removes whitespace
#group_name := "Index(25,19)"; #Use this if the group name for the primitive group has illegal characters for text files.

#Option 2: You can simply specify the group by generators in 'grp' and its group name in 'group_name'.
#grp := Group([ (3,7)(4,8)(5,6), (1,5,3,7,8,9)(2,6,4) ]);
#group_name := "AGammaL(1,9)Subg";

#NOTE: IF "group_name" ISN"T CHANGED BETWEEN DIFFERENT PROGRAM EXECUTIONS, IT MAY OVERWRITE PREVIOUS LP AND SUMMARY FILES


### TWEAKING THE OUTPUT (OPTIONAL) ###
#NOTES:
#   If you wish to find all Cameron-Liebler sets then 'ending_LP' should equal 'Int(n/2)'.
#   If you only wish to only find all maximum cocliques, 'ending_LP' should equal '1'.
#   If you wish to start from a different LP file, specify it using 'starting_LP'.

n := LargestMovedPoint(grp); #degree of grp
order := Size(grp); #order of grp
starting_LP := 1; #LP to start at
ending_LP := Int(n/2);

#CONTINUE FROM PREVIOUS COMPUTATION
#0 = Default. This will create all new files (LP, TeX, txt, log, etc.)
#1 = Use this if you had to stop the computation for some reason, but now you wish to continue without having to remake the LP files.
    # NOTE: This assumes that you've run the program with cont_program = 0 at least once (so all files required are generated already).
cont_program := 0;

#OPTIONS FOR THE NUMBER OF MAXIMUM CLIQUES THAT WILL BE GENERATED:
#0 = Default.
    #We automatically determine if the maximum cliques that are subgroups span a space of dimension |G| - (n-1)^2. 
    #If yes, it uses option #1; if no, it uses option #2.
#1 = Only maximum cliques that are subgroups (and their cosets) are generated - uses the 'sonata' package, which is preloaded (see README). THIS MAY CAUSE ERRORS WHEN IN ANALYSES.
#2 = All maximum cliques are generated - uses the 'grape' package, which is preloaded (see README). THIS MAY CAUSE ERRORS WHEN IN ANALYSES.
clique_type := 0;

#OPTIONS FOR SOLUTION RESTRICTIONS:
#0 = Default. The stabilizers and their cosets (i.e. all mappings that send i -> j for all i,j) are NOT restricted (Gurobi will find them as a solution).
#1 = The ij-mappings ARE restricted, so Gurobi won't find them. NOTE: If '1' is selected, Gurobi may output 0 solutions.
ij_mappings := 0;

#EXTRA TESTS (ANALYSES):
#Checks for the following:
    #1 group degree, order, transitivity, Frobenius or non-Frobenius, and the number of components in Gamma_G;
    #2 the strict EKR property and other relevant maximum coclique properties;
    #3 the dimensions of C', C, and W;
    #4 the spectrum of Gamma_G and whether or not the ratio bound holds with equality.
    
#Specify which questions you want to perform by inserting its number in the following list.
extra_tests := [1,2,3,4]; 

#TEX OUTPUT
#0 = default - a .tex file is generated with all information that is displayed to console (and a bit more).
#1 = no .tex file is generated.
tex_output := 0;


#Python scripts
python_gurobi_script := Concatenation(home_dir, "Python_Scripts/GurobiReader.py");
python_file_edit_script := Concatenation(home_dir, "Python_Scripts/LPEndingDel.py");
python_summary_script := Concatenation(home_dir, "Python_Scripts/SolutionSummary.py");
python_start_timer_script := Concatenation(home_dir, "Python_Scripts/StartTimer.py");
python_end_timer_script := Concatenation(home_dir, "Python_Scripts/EndTimer.py");

#Directories
gurobi_output_directory := Concatenation(home_dir, "Gurobi_LP_Models/");
gurobi_log_directory := Concatenation(home_dir, "Gurobi_Logs/");
solutions_summary_directory := Concatenation(home_dir, "Solutions_Summary_Files/");
temp_textfile_directory := Concatenation(home_dir, "Temp_Textfiles/");
latex_output_directory := Concatenation(home_dir, "LaTeX_Output/");

#Important files
solution_file := Concatenation(temp_textfile_directory, "Solution.txt");
start_time_file := Concatenation(temp_textfile_directory, "StartTime.txt");
solutions_summary_file := Concatenation(temp_textfile_directory, group_name, "_Solutions.txt");
solutions_tallied_file := Concatenation(temp_textfile_directory, group_name, "_Tallied.txt");
latex_file := Concatenation(latex_output_directory, group_name, "_Summary.tex");
analyses_file := Concatenation(home_dir, "Analyses.g");

#Use Python to time the program.
Exec(python3, python_start_timer_script, start_time_file);

Print("\n"); #Nicer console output

lp_files := [];
#Store lp files
for i in [1..ending_LP] do
    file_ending := Concatenation("_c_", String(i), ".lp");
    temp_file := Concatenation(gurobi_output_directory, group_name, file_ending);
    Add(lp_files, temp_file);
od;

#Determine the length of the "end" of each lp file generated
lp_file_end := 3; #This never changes for any group. It's used to remove the last three lines: "Binary, (variables), End"

#This stores each individual group element in a list called 'elements'
elements := Elements(grp);
trans := Transitivity(grp); #Determines the maximal nonnegative integer k such that the action of G is k-transitive.
aut := AutomorphismGroup(grp);
inn := InnerAutomorphismsAutomorphismGroup(aut); #Note that "aut" can potentially generate sets isomorphic to CL sets that are in the wrong module. Use "inn" instead.
###############################################################################################################


##### SECTION 1 - USER DEFINED FUNCTIONS #####
initializeLPFiles := function(G, degree, num, cliques, els, file, ij_map, stab_forb)
	local output, y, x, i, j, L, cur_pos, constraint;
	output := OutputTextFile(file, false);
	
	
	### Objective Section ###
  	WriteAll(output, "Maximize\n");
  	WriteAll(output, Concatenation(Concatenation(List([1..Size(els)], x -> Concatenation("x", String(x), " + "))), "NULL\n"));
	
	### Constraint Section ###
  	WriteAll(output, "Subject To\n");
	WriteAll(output, "x1 = 1\n");
	WriteAll(output, "NULL = 0\n");
 
	#Checks for integer intersection with cliques
  	for y in [1..Size(cliques)] do
        #Converts each maximum clique to a list of indices of 'elements'.
        #Because the cliques are sorted, we can tell GAP to search through 'elements' from the index of the previously indexed permutation.
        L := [];
        cur_pos := Position(els, cliques[y][1]);
        Add(L, cur_pos);
        for i in [2..Size(cliques[y])] do
            cur_pos := Position(els, cliques[y][i], cur_pos);
            Add(L, cur_pos);
        od;
        
        WriteAll(output, Concatenation("C", String(y), ": ", Concatenation(List(L, x -> Concatenation("x", String(x), " + "))), "NULL = ", String(num),"\n"));
  	od;
	
    if ij_map = 1 then
        for y in [1..Size(stab_forb)] do
            #Converts each permutation in stab_forb into a list of indices of 'elements'.
            L := [];
            cur_pos := Position(els, stab_forb[y][1]);
            Add(L, cur_pos);
            for i in [2..Size(stab_forb[y])] do
                cur_pos := Position(els, stab_forb[y][i], cur_pos);
                Add(L, cur_pos);
            od;
            
            WriteAll(output, Concatenation("M", String(y), ": ", Concatenation(List(L, x -> Concatenation("x", String(x), " + "))), "NULL <= ", String((Size(G)/degree)-1),"\n"));
        od;
    fi;

	CloseStream(output);
end;

endLPFiles := function(els, file)
	local output;
	output := OutputTextFile(file, true);
	
	### Integrality Section ###
    WriteAll(output, "Binary\n");
	WriteAll(output, Concatenation(Concatenation(List([1..Size(els)], x -> Concatenation("x", String(x), " "))), "NULL\n"));
	
	### Close Gurobi LP ###
	WriteAll(output, "End\n");

	CloseStream(output);
end;

solution_family := function(sol, els, inn_els, G, sol_sum_file, tallied_file, degree, ord, python_file, sol_dir, CL)
    #We find the orbit of the solution set on the inner automorphism group.
    #i.e. g in inn acts on the solution set via conjugation of each of the e in elements of the set: geg^(-1).
    local solution_orbit, family, x, coclique, output;
    output := OutputTextFile(sol_sum_file, true);
    
    solution_orbit := Orbit(inn_els, els{sol}, OnSets);
    family := Orbits(grp, solution_orbit, {S, g} -> AsSet(S*g));
    family := Concatenation(family);
    family := AsSet(List(family, AsSet)); #Sorting the solutions seems to speed of Gurobi's optimization
    
    #Print solution info to a text file as we find the solutions.
    #This is useful in the event you need to stop the program early. No solution found is lost.
    PrintToFormatted(output, "Solution: {1}\nSolution size: {2}\nSize of solution family: {3}\nCameron-Liebler set: {4}\n\n", sol, Size(sol), Size(family), CL);
    CloseStream(output);
    
    Exec(python3, python_file, String(degree), group_name, String(ord), sol_sum_file, sol_dir, tallied_file);
    
    return family;
end;

forbidLPFiles := function(els, forb, S_num, file, cont)
    local output, y, L, cur_pos, i;
    output := OutputTextFile(file, true);
    
    #This converts forbidden_final into the indices of els and writes it to the LP file in the proper format.
    for y in [1..Size(forb)] do
        #Converts each maximum coclique to a list of indices of 'elements'.
        #Because the cocliques are sorted, we can tell GAP to search through 'elements' from the index of the previously indexed permutation.
        L := [];
        cur_pos := Position(els, forb[y][1]);
        Add(L, cur_pos);
        for i in [2..Size(forb[y])] do
            cur_pos := Position(els, forb[y][i], cur_pos);
            Add(L, cur_pos);
        od;
        #The solutions are title "S". If we are continuing from a previous computation, we call them "F".
        if cont = 0 then
            WriteAll(output, Concatenation("S", String(S_num), ": ", Concatenation(List(L, x -> Concatenation("x", String(x), " + "))), "NULL <= ", String(Size(forb[y])-1),"\n"));
        else
            WriteAll(output, Concatenation("F", String(S_num), ": ", Concatenation(List(L, x -> Concatenation("x", String(x), " + "))), "NULL <= ", String(Size(forb[y])-1),"\n"));
        fi;
        
        S_num := S_num + 1;
    od;

    CloseStream(output);
end;

cliqueCheck := function(set, degree)
	local clique, x;

	clique := true;
    
    #We iterate over all pairs (a,b) of elements of set. Elements a and b are adjacent if ab^(-1) is a derangement.
    for x in IteratorOfCombinations(set, 2) do
        if NrMovedPoints(x[1]*x[2]^(-1)) < degree then
            clique := false;
            break;
        fi;
    od;
    
    return clique;
end;

cocliqueCheck := function(set, degree)
	local coclique, x;
    
    coclique := true;
    
    #We iterate over all pairs (a,b) of elements of set. Elements a and b are adjacent if ab^(-1) is a derangement.
    for x in IteratorOfCombinations(set, 2) do
        if NrMovedPoints(x[1]*x[2]^(-1)) = degree then
            coclique := false;
            break;
        fi;
    od;
    
    return coclique;
end;

moduleCheck := function(set, G)
    #We check if each maximum coclique lies in the same module.
    #For each irreducible representation there is a module. 
    #For a set (coclique) if the sum of the values of the character on the elements in the set is non-zero, then the “set is in the module”, if the sum is 0, then the set is not in the module.
    #For every coclique we should get that it is in two modules, the trivial module and the dimension (n-1)^2 module. (For any set S, the sum of S over the trivial module is |S|.)
    local M, i, temp_module, c, x, sum;
    M := [];
    for i in Irr(G) do
        temp_module := [];
        for c in set do
            sum := 0;
            for x in c do
                sum := sum + x^i;
            od;
            
            if sum <> 0 then
                Add(temp_module, AsSet(c));
            fi;
        od;
        
        Add(M, AsSet(temp_module));
    od;
    
    return M;
end;

cliqueSpan := function(set, G)
    local M, dim, i, temp_module, c, x, sum, M_dim;
    M := [];
    dim := [];
    M_dim := [];
    
    for i in Irr(G) do
        temp_module := [];
        for c in set do
            sum := 0;
            for x in c do
                sum := sum + x^i;
            od;
            
            if sum <> 0 then
                Add(temp_module, AsSet(c));
            fi;
        od;
        
        Add(M, AsSet(temp_module));
        
        #We keep track of the dimension of the modules that the cliques don't span.
        if Size(AsSet(temp_module)) = 0 then
            Add(dim, ()^i);
        fi;
    od;
    
    Add(M_dim, [M, dim]);
    return(M_dim);
end;

eValues := function(G, degree)
    ### For reference: See Theorem 11.12.3 in EKR Theorems: Algebraic Approaches (Godsil and Meagher). ###
	
    local chi, sum, der, c, eigenv_chi_chi_squared, evals, mult, evals_mult, eigen_mult, x, y, temp_mult, eigen, i;
    
    #We determine Der(G).
	der := Filtered(ConjugacyClasses(G), x -> NrMovedPoints(Representative(x)) = degree);
	
    #We find the eigenvalues of \Gamma_G
    eigenv_chi_chi_squared := [];
	for chi in Irr(G) do
		sum := 0;
		for c in der do
			sum := sum + Size(c)*(Representative(c))^chi;
		od;
		Add(eigenv_chi_chi_squared, [sum/()^chi, chi, (()^chi)^2]);
	od;
	Sort(eigenv_chi_chi_squared);
    
    #We compile the eigenvalues and their multiplicities in a nice way.
    evals := List(eigenv_chi_chi_squared, x -> x[1]);
    mult := List(eigenv_chi_chi_squared, x -> x[3]);
    evals_mult := [];
    for i in [1..Size(evals)] do
        Add(evals_mult, [evals[i], mult[i]]);
    od;

    eigen_mult := [];
    for x in evals_mult do
        eigen := x[1];
        temp_mult := 0;
        for y in evals_mult do
            if y[1] = eigen then
                temp_mult := temp_mult + y[2];
            fi;
        od;
        Add(eigen_mult, [eigen, temp_mult]);
    od;

    eigen_mult := AsSet(eigen_mult);

	return(eigen_mult);
end;

frobeniusCheck := function(spectrum)
    local frobenius;
    
    if Size(spectrum) = 2 then
        frobenius := true;
    else
        frobenius := false;
    fi;
    
    return(frobenius);
end;

split := function(size, num_cpu)
	local L, zp, pp, L_final, running_total, i;
	
    L := [];
    # If we cannot split the number into exactly 'N' parts
    if(size < num_cpu) then
        Print(-1);
    # If x % n == 0 then the minimum difference is 0 and all numbers are x / n
    elif (RemInt(size, num_cpu) = 0) then
        for i in [1..num_cpu] do
			Add(L, Int(size/num_cpu));
		od;
    else
        # up to n-(x % n) the values will be x / n and after that the values will be x / n + 1
        zp := num_cpu - RemInt(size, num_cpu);
        pp := Int(size/num_cpu);
        for i in [1..num_cpu] do
            if(i>= zp) then
				Add(L, pp+1);
            else
				Add(L, pp);
			fi;
		od;
	fi;

    L_final := [0];
    running_total := 0;
    for i in L do
        running_total := running_total + i;
	    Add(L_final, running_total);
	od;
    
    #Forces the last entry to be the order of the group. 
    #Without this, sometimes the last entry can be one too big.
    L_final[num_cpu+1] := size;
    return L_final;
end;

#####################################

if cont_program = 0 then

##### SECTION 2 -  PRELIMINARY CHECKS #####

Print("--- PRELIMINARY CHECKS (CONSOLE ONLY) ---\n\n");

### CHECK IF THE GROUP IS FROBENIUS ###

#We report the group name, degree, and size.
PrintFormatted("1. {1}; degree {2}; order {3}\n", group_name, n, order);

#We report the transitivity of the group.
PrintFormatted("2. The group is {1}-transitive.\n", trans);

#We check if G has two distinct eigenvalues. If it does, then it is a Frobenius group.
eMult := eValues(grp, n);
frobenius_group := frobeniusCheck(eMult);

if frobenius_group = true then
    Print("3. The group is Frobenius.\n");
else
    Print("3. The group is non-Frobenius.\n");
fi;

### CHECK IF \Gamma_G HAS MULTIPLE COMPONENTS ###

#The multiplicity of the largest eigenvalue is the number of connected components of \Gamma_G.
num_components := Maximum(eMult)[2];
if num_components > 1 then
    if frobenius_group = true then
        PrintFormatted("4. Gamma_G has {1} components, where each component is K_{2}.\n", num_components, n);
    else
        adj_rels := Filtered(ConjugacyClasses(grp), x -> NrMovedPoints(Representative(x))=n);
        H := GroupByGenerators(Union(adj_rels));
        H_gen := SmallGeneratingSet(H);
        H_trans := Transitivity(H);
        PrintFormatted("4. Gamma_G has {1} components. The subgroup H of derangements of G is {2}-transitive and its generators are: {3}\n", num_components, H_trans, H_gen);
        PrintFormatted("NOTE: We recommend terminating the current computation and restarting with: grp := H; group_name := {1}_Subgroup\n", group_name);
    fi;
else
    Print("4. Gamma_G has 1 component.\n");
fi;

### CHECK IF THE CLIQUES THAT ARE SUBGROUPS SPAN A SPACE OF DIMENSION |G| - (n-1)^2 ###

#The nature of this check depends on the user parameter "clique_type".

#If clique_type = 0, we determine whether to use just the cliques that are subgroups (and their cosets) or all maximum cliques.
#If clique_type = 1, we simply find the cliques that are subgroups (and their cosets) whether they span |G| - (n-1)^2 or not.
if clique_type <> 2 then
    #We determine which subgroups of grp form a clique for Gamma_G, while we don't actually compute Gamma_G.
    clique_subgroups := Filtered(Subgroups(grp), x -> Size(x) = n and cliqueCheck(AsList(x), n));; #The function 'Subgroups' uses the 'sonata' package.
    
    #We check if the cliques that are subgroups span a module of dimension |G| - (n-1)^2.
    #This is done by calling the function 'cliqueSpan'. We then generate the list 'zeroes_check_span'.
    modules_dimension := cliqueSpan(clique_subgroups, grp);
    module_cliques_span := modules_dimension[1][1];
    clique_vector_span := List(module_cliques_span, x -> Size(x));
    zeroes_check_span := Number(clique_vector_span, x -> x = 0);
    
    #The dimensions of the modules that the cliques don't span.
    dim_clique_vector := modules_dimension[1][2];
    
    #Add up squares of values of dim_clique and subtract from |G|.
    dim_clique := 0;
    for x in dim_clique_vector do
        dim_clique := dim_clique + x^2;
    od;
    
    dim_clique := order - dim_clique;
    dim_required := order - (n-1)^2;
     
    PrintFormatted("5. The cliques that are subgroups and their right cosets span a space of size {1}.\n", dim_clique);
    PrintFormatted("6. |G| - (n-1)^2 = {1}\n", dim_required);
    
    #If zeroes_check_span = 1, the subgroups span; otherwise, they don't span.
    if clique_type <> 1 then
        if zeroes_check_span = 1 then
            clique_type := 1; #This implies that the cliques that are subgroups span the space.
        else
            clique_type := 2; #This implies that the cliques that are subgroups don't span.
        fi;
    fi;
fi;

###############################################################################################################

##### SECTION 3 - DETERMINING THE MAXIMUM CLIQUES #####

#The nature of this section depends on: "clique_type".

#To find all maximum cliques, we form Gamma_G and use the "grape" package.
if clique_type = 2 then
    # the derangement graph of grp; this uses the GRAPE package
    Cay := CayleyGraph(grp, Filtered(grp, x -> NrMovedPoints(x) = n));

    # compute a set of maximal cliques of size n in Cay which is guaranteed to contain
    # at least one representative from each orbit of maximal cliques;
    # returned as lists of indices into Cay
    max_clique_indices := CompleteSubgraphs(Cay,n,1);

    # convert the vertices of Cay into permutations of grp.
    max_clique_perms := List(max_clique_indices, i -> AsSet(Cay.names{i}));

    # we want all maximum cliques, so compute the orbits of the orbit representatives;
    # we act on sets by right multiplication
    # the second (right) argument is an element of the acting group, and we find the image of the first argument under the action of the second
    maximum_clique_orbs := Orbits(grp, max_clique_perms, {cl,g} -> AsSet(cl*g));

    # finally merge all the orbits into one
    maximum_cliques := Concatenation(maximum_clique_orbs);
    maximum_cliques := AsSet(List(maximum_cliques, AsSet)); #Sorting the cliques seems to speed up Gurobi's optimization
else
    #Compute the right cosets of each subgroup in 'clique_subgroups'
    #Make use of 'RightCosets' GAP function, since all cliques are subgroups.
    maximum_cliques := [];
    for cl in clique_subgroups do
       r_cosets := RightCosets(grp, AsGroup(cl));
       for x in r_cosets do
           Add(maximum_cliques, AsList(x));
       od;
    od;
    maximum_cliques := AsSet(List(maximum_cliques, AsSet)); #This may not be necessary
fi;

###############################################################################################################



#The following list will hold all of the maximum cocliques at some point.
cocliques := [];

#If ij_mappings = 1, that means Gurobi won't find it. We add it as a solution.
    #Restricts all of the stabilizers and their cosets (mappings from i->j). These are all known cocliques.
    #Sometimes, the orbit + cosets of the a stabilizer contains many more than the standard n^2 solutions.
    #Sometimes the n^2 ij_mappings are contained in the orbit + cosets of more solutions.
    #To be safe, we manually compute the stabilizer of 1, find its orbit on inn, and compute its cosets.
stab1_orbit_cosets := [];
if ij_mappings = 1 then
    stab1_orbit := Orbit(inn, Elements(Stabilizer(grp, 1)), OnSets);
    
    #Compute the right cosets of each stabilizer in stab1_orbit.
    #Make use of 'RightCosets' GAP function, since all stabilizers are subgroups.
    for stab in stab1_orbit do
       r_cosets := RightCosets(grp, AsGroup(stab));
       for x in r_cosets do
           Add(stab1_orbit_cosets, AsList(x));
       od;
    od;
    
    stab1_orbit_cosets := AsSet(List(stab1_orbit_cosets, AsSet));
    Add(cocliques, stab1_orbit_cosets);
    
    #Convert the stabilizer of 1 into indices of elements (for text file).
    stab_sol := [];
    for x in stab1_orbit[1] do
        Add(stab_sol, Position(elements, x));
    od;
    
    #Write to the solutions text file.
    solutions_summary_output := OutputTextFile(solutions_summary_file, true); 
    PrintToFormatted(solutions_summary_output, "Solution: {1}\nSolution size: {2}\nSize of solution family: {3}\nCameron-Liebler set: true\n\n", stab_sol, Size(stab_sol), Size(stab1_orbit_cosets));
    CloseStream(solutions_summary_output);
    
    Exec(python3, python_summary_script, String(n), group_name, String(order), solutions_summary_file, solutions_summary_directory, solutions_tallied_file);
fi;

#Initialize LP File (and optionally restrict ij-mappings)
for c in [starting_LP..ending_LP] do
    initializeLPFiles(grp, n, c, maximum_cliques, elements, lp_files[c], ij_mappings, stab1_orbit_cosets);
od;

#End LP File
for c in [starting_LP..ending_LP] do
    endLPFiles(elements, lp_files[c]);
od;


fi; #Ends the cont_program condition.

sol_data := []; #This will hold all solutions and information pertaining to each.

if ij_mappings = 1 and cont_program = 0 then
    temp_data := [Size(stab_sol), Size(stab1_orbit_cosets), stab_sol, stab1_orbit_cosets];
    Add(sol_data, temp_data);
fi;

solutions := [];
zeroes_check := 0; #This initializes the count of the number of modules that don't contain a clique or coclique.
solution_count := 1;
starting_S := 1; #When writing to the LP files, we start with an 'S' number of 1 (this is incremented as we go).


### CALL GUROBI VIA PYTHON SCRYPT ITERATIVELY ###
for c in [starting_LP..ending_LP] do
    if c = 1 then
        ##### SECTION 4 - FINDING THE MAXIMUM COCLIQUES VIA LINEAR PROGRAMMING #####
        Print("\n--- FINDING THE MAXIMUM COCLIQUES ---\n\n");
    elif starting_LP = 1 and c = 2 then
        ##### SECTION 5 - TEX OUTPUT: ANALYSES OF GROUP AND MAXIMUM COCLIQUES #####
        Read(analyses_file);

        #group_name := ReplacedString(group_name, " ", "_"); #Replaces spaces in group name for writing to text files
        
        #In case group_name contains backslashes or whitespaces, we remove them.
        RemoveCharacters(group_name, " "); #Removes whitespace
        RemoveCharacters(group_name, "\\"); #Removes backslashes

        ##### SECTION 6 - FINDING ALL CAMERON-LIEBLER SETS #####
        Print("\n--- FINDING CAMERON-LIEBLER SETS ---\n\n");
    fi;
    
    Exec(python3, python_gurobi_script, lp_files[c], group_name, String(solution_count), temp_textfile_directory, gurobi_log_directory, String(c));
    Read(solution_file); #This is stored in the variable: 'solution'
    
    while solution <> 0 do
        solution_count := solution_count + 1;
        
        #For basic analysis (extra_tests 4 isn't selected and starting_LP > 1) we assume the solution found is a Cameron-Lieber set.
        CL_set := true;
        
        #If starting_LP = 1, analysis question 4 has been enabled, and there exists module(s) that don't contain a clique or coclique, we can check if the solution found is a CL set.
        if starting_LP = 1 and (4 in extra_tests) and zeroes_check > 0 then
            solution_perm := elements{solution}; #Convert the solution to a list of permutations in grp.
            nonspan_irr := Positions(cocl_cl_vector, 0); #Indices of Irr(G) that don't contain a clique or coclique.
            
            for i in nonspan_irr do
                s := 0;
                for x in solution_perm do
                    s := s + x^Irr(grp)[i];
                od;
                
                if s <> 0 then
                    break;
                fi;
            od;
            
            CL_set := false;
        fi;
        
        #Find the family of the solution.
        forbidden_final := solution_family(solution, elements, inn, grp, solutions_summary_file, solutions_tallied_file, n, order, python_summary_script, solutions_summary_directory, CL_set);

        #We use Python to delete the ending of each lp file that is and will be used (this is faster and easier to do in Python).
        Exec(python3, python_file_edit_script, String(starting_LP), String(ending_LP), gurobi_output_directory, group_name, String(lp_file_end), String(c));
            
        #Write the 'S' lines to the current c and larger c LP files
        for d in [c..ending_LP] do
            forbidLPFiles(elements, forbidden_final, starting_S, lp_files[d], cont_program);
        od;

        #End LP File
        for d in [c..ending_LP] do
            endLPFiles(elements, lp_files[d]);
        od;
        
        #Increment starting_S
        starting_S := starting_S + Size(forbidden_final);
        
        #Add forbidden_final to 'cocliques' to be analyzed later; c = 1 corresponds to maximum cocliques.
        if c = 1 then
            Add(cocliques, forbidden_final);
        fi;
        
        #Otherwise, we add it to 'sol_data'; c > 1 corresponds to Cameron-Lieber sets that aren't maximum cocliques.
        if CL_set = true then
            #We store all solutions, their lengths, and the size of their orbits + cosets in the following list of lists; i.e. sol_data := [[solution length, orbit+cosets, solution], etc.].
            temp_data := [Size(solution), Size(forbidden_final), solution, forbidden_final];
            Add(sol_data, temp_data);
            Add(solutions, elements{solution}); #This may be useful for debugging later on. It is unused.
        else
            Print("NOTE: Solution above is NOT a Cameron-Lieber set.\n\n");
        fi;
        
        #Call Gurobi to analyze the same LP file again.
        Exec(python3, python_gurobi_script, lp_files[c], group_name, String(solution_count), temp_textfile_directory, gurobi_log_directory, String(c));
        Read(solution_file); #This is stored in the variable: 'solution'
    od;
od;

#Perform analyses if starting_LP = ending_LP = 1
if starting_LP = 1 and ending_LP = 1 then
    ##### SECTION 5 - TEX OUTPUT: ANALYSES OF GROUP AND MAXIMUM COCLIQUES #####
    Read(analyses_file);
fi;
###############################################################################################################

##### SECTION 7 - TEX OUTPUT: REPORT CAMERON-LIEBLER SETS #####

if tex_output = 0 then
    latex_output := OutputTextFile(latex_file, true); #Continues the analyses .tex file
    #group_name := ReplacedString(group_name, "_", " "); #Removes underscore(s) in group name for nicer tex output
    #In case group_name contains backslashes or whitespaces, we remove them.
    RemoveCharacters(group_name, " "); #Removes whitespace
    RemoveCharacters(group_name, "\\"); #Removes backslashes
    
    #In case group_name contains greek letters, we replace them with nice tex output.
    group_name := ReplacedString(group_name, "Gamma", "\\Gamma ");
    group_name := ReplacedString(group_name, "Sigma", "\\Sigma ");
    
    WriteAll(latex_output, "\\noindent \\textbf{Cameron-Liebler sets} \n\n");

    #Title and value(s) c used
    if starting_LP = ending_LP then
        if Int(n/2) = 1 and ending_LP = 1 and ij_mappings = 0 then
            PrintToFormatted(latex_output, "\\noindent All non-trivial Cameron-Liebler sets (including the canonical cocliques) were searched for and found (${1} \\leq c \\leq {2}$). ", starting_LP, ending_LP);
        elif Int(n/2) = 1 and ending_LP = 1 then
            PrintToFormatted(latex_output, "\\noindent All non-trivial Cameron-Liebler sets were searched for and found (${1} \\leq c \\leq {2}$). ", starting_LP, ending_LP);
        elif ending_LP = 1 then
            PrintToFormatted(latex_output, "\\noindent Only maximum cocliques were searched for and found ($c = {1}$). ", starting_LP);
        else
            PrintToFormatted(latex_output, "\\noindent Only Cameron-Liebler sets with $c = {1}$ were searched for and found. ", starting_LP);
        fi;
    elif ij_mappings = 0 then
        PrintToFormatted(latex_output, "\\noindent All non-trivial Cameron-Liebler sets (including the canonical cocliques) were searched for and found (${1} \\leq c \\leq {2}$). ", starting_LP, ending_LP);
    else
        PrintToFormatted(latex_output, "\\noindent All non-trivial Cameron-Liebler sets were searched for and found (${1} \\leq c \\leq {2}$). ", starting_LP, ending_LP);
    fi;

    PrintToFormatted(latex_output, "The total number of families of solutions $\\mathcal{{F}}$ found: {1}. \n\n", Size(sol_data));

    ### Insert table that summarizes the families of Cameron-Liebler sets found ###
    #Form start of table
    WriteAll(latex_output, "\\FloatBarrier\n");
    WriteAll(latex_output, "\\begin{table}[h] \n\\centering \n\\begin{threeparttable}\n");
    PrintToFormatted(latex_output, "\\caption{{Cameron-Liebler sets of ${1}$}}\n", group_name);
    
    #In case group_name contains greek letters, we replace them back with non-tex input.
    group_name := ReplacedString(group_name, "\\Gamma ", "Gamma");
    group_name := ReplacedString(group_name, "\\Sigma ", "Sigma");
        
    PrintToFormatted(latex_output, "\\label{{Ch4: Table: Cameron-Liebler sets of {1}}}\n", group_name);
    WriteAll(latex_output, "\\begin{tabular}{cccc} \n\\toprule\n");
    WriteAll(latex_output, "$c$ & Weight & Number of families & Size of each family \\\\ \n\\midrule\n");

    #Table data
    data_orb := AsSet(Collected(List(sol_data, x -> [x[1],x[2]])));
    data_len := Collected(List(sol_data, x -> x[1]));
    for c in [starting_LP..ending_LP] do
        empty_line := true;
        for L in data_len do
            if c * (order/n) = L[1] then
                #We isolate the orbit/coset data with length L[1]
                orb_coset := Filtered(data_orb, x -> x[1][1] = L[1]);
                PrintToFormatted(latex_output, "{1} & {2} & {3} & \\makecell[ct]{{", c, L[1], L[2]);
                for i in [1..Size(orb_coset)] do
                    if i <> Size(orb_coset) then
                        PrintToFormatted(latex_output, "{1} of size {2} \\\\ \n", orb_coset[i][2], orb_coset[i][1][2]);
                    else
                        if c = Int(n/2) then
                            PrintToFormatted(latex_output, "{1} of size {2}}} \\\\ \n", orb_coset[i][2], orb_coset[i][1][2]);
                        else
                            PrintToFormatted(latex_output, "{1} of size {2}}} \\\\ \\\\ \n", orb_coset[i][2], orb_coset[i][1][2]);
                        fi;
                    fi;
                od;
                empty_line := false;
            fi;
        od;
        if empty_line = true then
            if c = Int(n/2) then
                PrintToFormatted(latex_output, "{1} & {2} & {3} & {4} \\\\ \n", c, c * (order/n), 0, 0);
            else
                PrintToFormatted(latex_output, "{1} & {2} & {3} & {4} \\\\ \\\\ \n", c, c * (order/n), 0, 0);
            fi;
        fi;
    od;

    #End table
    WriteAll(latex_output, "\\bottomrule \n\\end{tabular} \n\\end{threeparttable} \n\\end{table}\n\n");
    WriteAll(latex_output, "\\FloatBarrier\n");
    WriteAll(latex_output, "\\end{document}\n");
    CloseStream(latex_output);
fi;
###############################################################################################################

#We use Python to time the entire program and signal the end of computation.
Exec(python3, python_end_timer_script, start_time_file, solution_file, solutions_summary_file, solutions_tallied_file);





# ### OPTIONAL CHECK ###
# #Results of this check can be found in the Temp_Textfiles folder.
# #We check if each solution is the union of the cosets of a subgroup contained in the solution


# solutions_summary_output := OutputTextFile(solutions_tallied_file, true); 
# G_subgroups := List(Filtered(Subgroups(grp), x -> (Size(x) < order and Size(x) > 2)), x -> Elements(x));
# number_generated_subg := 0;

# for solution in solutions do
    # subg_bool := false;
    
    # #Filter out subgroups that aren't subsets of the solution
    # contained := [];
    # for i in [1..Size(G_subgroups)] do
        # if IsSubset(solution, G_subgroups[i]) then
            # Add(contained, G_subgroups[i]);
        # fi;
    # od;

    # num_success := 0;
    # subg_size := [];
    # num_fam_success := 0;
    # for H in contained do
        # # BUILD SOLUTION FROM A SUBGROUP THAT INTERSECTS THE SOLUTION AND ITS RIGHT COSETS
        # H := AsGroup(H);
        # r_cosets := RightCosets(grp, H);
        # L := Filtered(r_cosets, x -> Size(Intersection(x, solution)) = Size(H));
        # if solution = Union(L) then
            # num_success := num_success + 1;
            # Add(subg_size, Size(H));
            # subg_bool := true;
        # fi;
    # od;
    
    # if subg_bool then
        # number_generated_subg := number_generated_subg + 1;
    # fi;
# od;

# Print("\n\n--- GENERATING CL-SETS RESULTS --- \n");
# PrintFormatted("Total number of solutions {1}.\n", Size(solutions));
# PrintFormatted("Number of solutions generated successfully using cosets of a subgroup: {1}\n", number_generated_subg);

# WriteAll(solutions_summary_output, "Summary:\n\n");
# PrintToFormatted(solutions_summary_output,"Total number of solutions {1}.\n", Size(solutions));
# PrintToFormatted(solutions_summary_output,"Number of solutions generated successfully using cosets of a subgroup: {1}\n", number_generated_subg);
# CloseStream(solutions_summary_output);
















