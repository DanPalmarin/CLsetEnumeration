
##### SECTION 5 - TEX OUTPUT: ANALYSIS OF GROUP AND MAXIMUM COCLIQUES #####
#These can be toggled off with the 'extra_tests' variable at the top of MAIN.g.


### ANALYSIS OF GROUP AND MAXIMUM COCLIQUES ###

if Size(extra_tests) > 0 then
    Print("--- PERFORMING ANALYSIS OF MAXIMUM COCLIQUES ---\n\n");
    
    #We start by confirming that each set in 'cocliques' is indeed a maximum coclique. If it isn't, we filter it out.
    maximum_cocliques_unchecked := AsSet(List(Concatenation(cocliques), AsSet)); #Note: 'Concatenation' flattens 'cocliques'; i.e. Concatentation([[1],[2]]) = [1,2].
    maximum_cocliques := Filtered(maximum_cocliques_unchecked, x -> cocliqueCheck(x, n));
    
    if Size(maximum_cocliques_unchecked) <> Size(maximum_cocliques) then
        num_other_solutions := Size(maximum_cocliques_unchecked) - Size(maximum_cocliques);
    fi;
    
    #Questions 1 and 2 both require the maximum cocliques that are subgroups, cosets, or neither.
    #Hence, we find them here.
    if 1 in extra_tests or 2 in extra_tests then
        #Manually compute the stabilizers and their cosets to compare with 'maximum_cocliques'
        stabs := [];
        for i in [1..n] do
            Add(stabs, Elements(Stabilizer(grp, i)));
        od;

        stabs_cosets := [];
        for cocl in stabs do
           r_cosets := RightCosets(grp, AsGroup(cocl));
           for x in r_cosets do
               Add(stabs_cosets, AsList(x));
           od;
        od;

        stabs_cosets := AsSet(List(stabs_cosets, AsSet));
        
        #Use Sonata to find all subgroups of size |grp|/n and then check for intersection with 'cocliques'.
        cocliques_subg := AsSet(List(Filtered(Subgroups(grp), x -> Size(x) = (order / n) and AsSet(x) in maximum_cocliques), x -> Elements(x)));
        cocliques_unknown := AsSet(Difference(maximum_cocliques, cocliques_subg));

        #For those maximum cocliques that aren't subgroups, we check if they are cosets of a subgroup.
        #We do this by pulling out some element (we choose the 1st one) from the set, find its inverse, and multiply every element of the set by it.
        #We check if this new set is in 'cocliques_subg'.
        cocliques_cosets := [];
        for i in [1..Size(cocliques_unknown)] do
            inv := Inverse(cocliques_unknown[i][1]);
            temp := AsSet(inv*cocliques_unknown[i]);
            if temp in cocliques_subg then
                Add(cocliques_cosets, cocliques_unknown[i]);
            fi;
        od;

        #Sort and remove duplicates (maintain mutability).
        cocliques_cosets := Set(List(cocliques_cosets, AsSet));
        cocliques_unknown := Set(List(cocliques_unknown, AsSet));

        #We store the maximum_cocliques that are neither groups nor cosets (if any) as 'other'.
        cocliques_other := Difference(cocliques_unknown, cocliques_cosets);
    
        #Partitions out the solutions that aren't in stabs_cosets. Further determines which ones are subgroups.
        unknown_solutions := AsSet(List(Difference(maximum_cocliques, stabs_cosets), AsSet));
        unknown_solutions_subg := AsSet(List(Intersection(unknown_solutions, cocliques_subg), AsSet));
        unknown_solutions_cosets := AsSet(List(Difference(unknown_solutions, unknown_solutions_subg), AsSet));
    fi;
fi;

### Question 1 ###
    #1 group degree, order, transitivity, Frobenius or non-Frobenius, and the number of components in Gamma_G;

if 1 in extra_tests then
    Print("Question 1.\n");
    #1a
    PrintFormatted("The group has degree {1} and order {2}.\n", n, order);

    #1b
    PrintFormatted("The group is {1}-transitive.\n", trans);

    #1c
    if frobenius_group = true then
        PrintFormatted("The group is Frobenius; Gamma_G has {1} components, where each component is K_{2}.\n", num_components, n);
    else
        PrintFormatted("The group is non-Frobenius; Gamma_G has {1} component(s).\n", num_components);
    fi;
fi;
    
### QUESTION 2 ###
    #2 the strict EKR property and other relevant maximum coclique properties;

if 2 in extra_tests then
    Print("\nWorking on Question 2...\n");
    #1a is done above and stored in "maximum_cocliques".
    #1b is also mostly done above in "stabs_cosets".
    #Does the group have the strict EKR property (stabilizers and their cosets are the only maximum cocliques)?
    strict_EKR := maximum_cocliques = stabs_cosets; #true or false
    
    #The following is a fringe case that is interesting to study.
    if strict_EKR = false and frobenius_group = false and num_components = 1 then
        #If G is not Frobenius (and has more than 1 solution) we make no conclusions regarding the unknown cocliques.
        #We organize according to solution count. If it is one, then that implies that there is a single orbit + cosets that generates all maximum cocliques.
        #In this case, we attempt to see if there are 'isomorphic' families that make up the cocliques. Each family would probably be isomorphic to the stabilizers and their cosets.
        #Ex: ASL(2,4) of order 960 has 1024 = 4(16^2) maximum cocliques; i.e. |stabs_cosets| = 16^2 and three other families 'isomorphic' to it. 
        
        #Checks if the size of stabs_cosets divides the size of maximum_cocliques evenly.
        num_pieces := Size(unknown_solutions) / Size(stabs_cosets);
        chk := false;
        
        if IsInt(num_pieces) = true then
            unk_sol_subg_cosets := [];
            for cocl in unknown_solutions_subg do
                temp_cosets := [];
                r_cosets := RightCosets(grp, AsGroup(cocl));
                for x in r_cosets do
                    Add(temp_cosets, AsList(x));
                od;
                Add(unk_sol_subg_cosets, temp_cosets);
            od;
            
            unk_sol_subg_cosets := AsSet(List(unk_sol_subg_cosets, AsSet));
            
            #Attempt to evenly partition maximum_cocliques with the subgroups and cosets
            maximum_cocliques_partitions := [];
            unk_sol_subg_partition := split(Size(unknown_solutions_subg), num_pieces);
            for i in [1..num_pieces] do
                Add(maximum_cocliques_partitions, AsSet(Concatenation(unk_sol_subg_cosets{[(unk_sol_subg_partition[i]+1)..unk_sol_subg_partition[i+1]]})));
            od;
            
            Add(maximum_cocliques_partitions, stabs_cosets);
            maximum_cocliques_partitions := AsSet(List(maximum_cocliques_partitions, AsSet));
            size_cl_partitions := Size(maximum_cocliques_partitions); #For convenience when writing to files later
            
            #Perform various checks on each partition in maximum_cocliques_partitions.
            for x in IteratorOfCombinations(maximum_cocliques_partitions, 2) do
                #Check that each partition is the same size.
                if Size(x[1]) = Size(x[2]) then
                    chk := true;
                else
                    chk := false;
                fi;
                
                #Check that each partition has the same number of subgroups. Assume true and change size_chk to false otherwise.
                if Size(Intersection(cocliques_subg, AsSet(List(x[1], AsSet)))) = Size(Intersection(cocliques_subg, AsSet(List(x[2], AsSet)))) then
                    chk := true;
                else
                    chk := false;
                fi;
            od;
        fi;
    fi;
    
    ### Reporting Question 2 results ###
    Print("\nQuestion 2.\n");
    if maximum_cocliques_unchecked <> maximum_cocliques then
        PrintFormatted("a. There are {1} maximum cocliques and {2} other solutions of the same length that are not maximum cocliques.\n", Size(maximum_cocliques), num_other_solutions);
    else
        PrintFormatted("a. There are {1} maximum cocliques.\n", Size(maximum_cocliques));
    fi;
    
    if strict_EKR = true then
        Print("b. The strict EKR property holds; i.e. the canonical cocliques are the only maximum cocliques of Gamma_G.\n\n");
    else
        Print("b. The strict EKR property does not hold; i.e. there are other maximum cocliques besides the canonical cocliques. ");
        if frobenius_group = true then
            PrintFormatted("The group is Frobenius; hence, Gamma_G is copies of K_{1} and so the maximum cocliques are known.\n\n", n);
        elif num_components > 1 then
            PrintFormatted("As Gamma_G has {1} components, the composition of the non-canonical cocliques is known.\n\n", num_components);
        else
            if chk = true then
                PrintFormatted("Maximum cocliques are split into {1} parts, each of size {2} containing {3} subgroups and {4} cosets.\n\n", size_cl_partitions, Size(stabs_cosets), Size(cocliques_subg)/size_cl_partitions, Size(cocliques_cosets)/size_cl_partitions);
            elif chk = false then
                PrintFormatted("Besides the canonical cocliques, it is unclear what the {1} other maximum cocliques are.\n\n", Size(unknown_solutions));
            fi;
        fi;
    fi;
fi;

# Question 2 continued...
if 2 in extra_tests then
    #2a. Answered in the "fringe" case of 1b.
    #2b.
    Print("Working on Question 2 (part 2)...\n");
    #We check if all of the cocliques that are subgroups are isomorphic.
    subg_iso := [];
    subg_not_iso := [];
    for x in IteratorOfCombinations(cocliques_subg, 2) do
        if IsIsomorphicGroup(AsGroup(x[1]), AsGroup(x[2])) = true then
            Add(subg_iso, x);
        else
            Add(subg_not_iso, x);
        fi;
    od;
    
    iso := Size(subg_iso) = NrCombinations(cocliques_subg, 2); #true or false
    structure := StructureDescription(AsGroup(cocliques_subg[1]));
    
    #2c.
    #We first determine how many conjucacy classes there are among the maximum cocliques.
    all_conj_classes := ConjugacyClassesSubgroups(grp);
    conj_classes := Filtered(all_conj_classes, x -> Size(Representative(x)) = order / n and cocliqueCheck(AsList(Representative(x)), n));
    
    ### Reporting Question 2 results ###
    
    #2c
    PrintFormatted("c. The number of maximum cocliques that are subgroups: {1}; cosets: {2}; neither: {3}.\n", Size(cocliques_subg), Size(cocliques_cosets), Size(cocliques_other));

    #2d
    if iso = true then
        PrintFormatted("d. All maximum cocliques that are subgroups are isomorphic. Structure description of these subgroups: {1}\n", structure);
    else
        Print("d. Not all maximum cocliques that are subgroups are isomorphic. ");
        PrintFormatted("Number of isomorphic pairs of maximum cocliques that are subgroups: {1}. ", Size(subg_iso));
        PrintFormatted("Number of non-isomorphic pairs of maximum cocliques that are subgroups: {1}.\n", Size(subg_not_iso));
    fi;

    #2e
    PrintFormatted("e. The number of conjucacy classes among the maximum cocliques that are subgroups: {1}\n", Size(conj_classes));
fi;  

### QUESTION 3 ###
    #3 the dimensions of C', C, and W;

if 3 in extra_tests then
    Print("\nWorking on Question 3...\n");
    if trans >= 2 and clique_type <> 2 and dim_clique = dim_required then
        zeroes_check := 0;
        module_cocliques_bool := true;
        module_cliques_bool := true;
    elif trans >= 2 then
        module_cocliques := cliqueSpan(maximum_cocliques, grp);
        module_cliques := cliqueSpan(maximum_cliques, grp);
        
        #The dimensions of the modules that the cliques and cocliques don't span.
        dim_clique_vector := module_cliques[1][2];
        dim_coclique_vector := module_cocliques[1][2];
        
        #Add up squares of values of dim_all_cliques and subtract from |G|.
        dim_all_cliques := 0;
        for x in dim_clique_vector do
            dim_all_cliques := dim_all_cliques + x^2;
        od;
        
        dim_all_cliques := order - dim_all_cliques;
        
        dim_coclique := 0;
        for x in dim_coclique_vector do
            dim_coclique  := dim_coclique  + x^2;
        od;
        
        dim_coclique := order - dim_coclique;
        
        
        dim_required := order - (n-1)^2;
        
        if dim_coclique = dim_required + 1 then
            module_cocliques_bool := true;
        else
            module_cocliques_bool := false;
        fi;
        
        if dim_all_cliques = dim_required then
            module_cliques_bool := true;
            zeroes_check := 0;
        else
            module_cliques_bool := false;
            zeroes_check := Size(module_cliques[1][2])-1;
        fi;
    else
        #Here, it is not guaranteed that the maximum cocliques are in the same module.
        
        module_cocliques := moduleCheck(maximum_cocliques, grp);
        module_cliques := moduleCheck(maximum_cliques, grp);
        
        #We store the size of each member of module_cocliques and module_cliques in two separate lists (for reference).
        mod_cocliques := List(module_cocliques, x -> Size(x));
        mod_cliques := List(module_cliques, x -> Size(x));
        
        #We check which members of Irr(grp) yielded maximum cocliques with a sum <> 0.
        module_count := 0;
        not_module_count := 0;
        for x in module_cocliques do
            if x = maximum_cocliques then
                module_count := module_count + 1;
            elif Size(x) = 0 then
                not_module_count := not_module_count + 1;
            fi;
        od;
        
        #We check which members of Irr(grp) yielded maximum cliques with a sum <> 0.
        clique_zero := 0;
        for x in module_cliques do
            if Size(x) = 0 then
                clique_zero := clique_zero + 1;
            fi;
        od;
        
        coclique_vector := List(module_cocliques, x -> Size(x));
        clique_vector := List(module_cliques, x -> Size(x));
        cocl_cl_vector := coclique_vector + clique_vector; #Adds corresponding elements in both lists
        zeroes_check := Number(cocl_cl_vector, x -> x = 0);
        #We check if the modules that contain cocliques are the modules that correspond to the irreducibles representations in the permutation module
        #(so the rep that gives the number of fixed points).
        irr_perm := [];
        for x in Irr(grp) do
            sum := 0;
            for g in ConjugacyClasses(grp) do
                sum := sum + Size(g) * ((Representative(g)^x) * (n-NrMovedPoints(Representative(g))));
            od;
            Add(irr_perm, sum/order);
        od;
        
        if Number(irr_perm + mod_cocliques, x -> x = 0) = Number(irr_perm, x -> x = 0) then
            irr_perm_bool := true;
        else
            irr_perm_bool := false;
        fi;
    fi;
    
    ### Reporting Question 3 results ###

    Print("\nQuestion 3.\n");
    if trans >= 2 then
        #3a
        PrintFormatted("a. There are {1} modules that contain neither a maximum clique nor a maximum coclique.\n", zeroes_check); 
        
        #3b
        if module_cocliques_bool = true then
            Print("b. All maximum cocliques are contained in the same module.\n");
        else
            Print("b. Not all maximum cocliques lie in the same module.\n");
        fi;
        
        #3c
        if module_cliques_bool = true then
            Print("c. All maximum cliques span every module except the non-trivial module that contains the maximum cocliques.\n");
        else
            if clique_type = 1 then
                Print("c. Excluding the non-trivial module that contains the maximum cocliques, the maximum cliques do not span all of the other modules. Note that only the maximum cliques that are subgroups (along with their cosets) were used.\n");
            else
                Print("c. Excluding the non-trivial module that contains the maximum cocliques, the maximum cliques do not span all of the other modules.\n");
            fi;
        fi;
        
        #3d
        if dim_clique = dim_required then
            PrintFormatted("d. The maximum cliques that are subgroups span a space of dimension {1} and so they span a module of dimension |G| - (n-1)^2 = {2}; hence, only these maximum cliques are required to find all Cameron-Liebler sets.\n", dim_clique, dim_required);
        else
            PrintFormatted("d. The maximum cliques that are subgroups do not span a module of dimension |G| - (n-1)^2 = {1}; hence, more (perhaps all) maximum cliques are required to find all Cameron-Liebler sets.\n", dim_required);
        fi;
    else
        if irr_perm_bool = true then
            PrintFormatted("There are {1} modules that contain neither a maximum clique nor a maximum coclique. ", zeroes_check);
            Print("The modules that contain maximum cocliques are the same modules that correspond to the irreducible representations in the permutation module.\n");
        else
            PrintFormatted("There are {1} modules that contain neither a maximum clique nor a maximum coclique. ", zeroes_check);
            Print("There is at least one module that contains maximum cocliques that does not correspond to any of the irreducible representations in the permutation module.\n");
        fi;
    fi;
fi;

### QUESTION 4 ###
    #4 the spectrum of Gamma_G and whether or not the ratio bound holds with equality.

if 4 in extra_tests then
    Print("\nWorking on Question 4...\n");
    #4a (contained in the variable "eMult")
    
    #4b: Ratio bound: \alpha(G) <= n(-Min(evals)) / d - Min(evals), where d is the regularity = Max(evals).
    alpha_X := order/n; #Size of a maximum coclique in grp
    ratio_bound := ((order * (-1*Minimum(eMult)[1])) / (Maximum(eMult)[1] - Minimum(eMult)[1]));
    
    if alpha_X = ratio_bound then
        ratio := 1;
    elif alpha_X < ratio_bound then
        ratio := 0;
    else
        ratio := -1;
    fi;
    
    #4c: Inertia bound test: \alpha(G) <= n_0 + min{n_+, n_-}, where n_0, n_+, and n_- represent the number of positive, negative, and zero eigenvalues of \Gamma_G.
    alpha_X := order/n; 
    pos_evals := 0;
    neg_evals := 0;
    zero_evals := 0;
    for x in eMult do
        if x[1] > 0 then
            pos_evals := pos_evals + x[2];
        elif x[1] < 0 then
            neg_evals := neg_evals + x[2];
        elif x[1] = 0 then
            zero_evals := zero_evals + x[2];
        fi;
    od;
    
    inertia_bound := zero_evals + Minimum(pos_evals, neg_evals);
    if alpha_X = inertia_bound then
        inertia := 1;
    elif alpha_X < inertia_bound then
        inertia := 0;
    else
        inertia := -1;
    fi;
    
    ### Reporting Question 4 results ###

    Print("\nQuestion 4.\n");
    #4a
    PrintFormatted("a. The eigenvalues (and their multiplicities) of Gamma_G: {1}\n", eMult);
    
    #4b
    if ratio = 1 then
        Print("b. The ratio bound holds with equality.\n");
    elif ratio = 0 then
        PrintFormatted("b. The ratio bound holds but not with equality: {1} < {2}\n", alpha_X, ratio_bound);
    else
        Print("b. The ratio bound does not hold.\n");
    fi;
    
    #4c
    if inertia = 1 then
        Print("b. The inertia bound holds with equality.\n");
    elif inertia = 0 then
        PrintFormatted("b. The inertia bound holds but not with equality: {1} < {2}\n", order/n, inertia_bound);
    else
        Print("b. The inertia bound does not hold.\n");
    fi;
fi;


### TEX OUTPUT ###
if tex_output = 0 then
    latex_output := OutputTextFile(latex_file, false); #Starts a new .tex file
    #group_name := ReplacedString(group_name, "_", " "); #Removes underscore(s) in group name for nicer tex output
    
    #In case group_name contains greek letters, we replace them with nice tex output.
    group_name := ReplacedString(group_name, "Gamma", "\\Gamma ");
    group_name := ReplacedString(group_name, "Sigma", "\\Sigma ");
    
    #TeX preamble
    WriteAll(latex_output, "\\documentclass[12pt]{report}\n\n");
    WriteAll(latex_output, "\\usepackage[margin=1.0in]{geometry}\n");
    WriteAll(latex_output, "\\usepackage{booktabs}\n");
    WriteAll(latex_output, "\\usepackage{makecell}\n");
    WriteAll(latex_output, "\\usepackage{threeparttable}\n");
    WriteAll(latex_output, "\\usepackage{placeins}\n");
    WriteAll(latex_output, "\\usepackage{amsmath}\n");
    WriteAll(latex_output, "\\usepackage{enumitem}\n");
    WriteAll(latex_output, "\\usepackage{enumitem}\n");
    WriteAll(latex_output, "\\newcommand{\\abs}[1]{\\lvert#1\\rvert}\n\n");
    
    WriteAll(latex_output, "\\begin{document}\n");
    PrintToFormatted(latex_output, "\\begin{{large}} \\noindent  $\\boldsymbol{{{1}}}$ \\end{{large}}\n\n", group_name);

    if Size(extra_tests) > 0 then
        WriteAll(latex_output, "\\noindent \\textbf{Properties} \\vspace{-8pt}\n");
        WriteAll(latex_output, "\\begin{enumerate}[leftmargin=*, label=(\\arabic*)]\n");
    fi;
    
    #Question 1
    if 1 in extra_tests then
        if frobenius_group then
            PrintToFormatted(latex_output, "\\item $\\abs{{G}}={1}$; $\\deg(G)={2}$; {3}-transitive; Frobenius; $\\Gamma_G$ has {4} components, where each component is $K_{2}$.\n", order, n, trans, num_components);
        else
            if num_components > 1 then
                PrintToFormatted(latex_output, "\\item $\\abs{{G}}={1}$; $\\deg(G)={2}$; {3}-transitive; non-Frobenius; $\\Gamma_G$ has {4} components.\n", order, n, trans, num_components);
            else
                PrintToFormatted(latex_output, "\\item $\\abs{{G}}={1}$; $\\deg(G)={2}$; {3}-transitive; non-Frobenius; $\\Gamma_G$ has {4} component.\n", order, n, trans, num_components);
            fi;
        fi;
    fi;
    
    #Question 2
    if 2 in extra_tests then
        if strict_EKR then
            PrintToFormatted(latex_output, "\\item The strict EKR property holds ({1} maximum cocliques).", Size(maximum_cocliques));
            if maximum_cocliques_unchecked <> maximum_cocliques then
                PrintToFormatted(latex_output, " There are {1} other solutions of the same length that are not maximum cocliques.", num_other_solutions);
            fi;
            WriteAll(latex_output, "\n");
        else
            PrintToFormatted(latex_output, "\\item The strict EKR property does not hold ({1} maximum cocliques). ", Size(maximum_cocliques));
            if frobenius_group = true then
                PrintToFormatted(latex_output, "The group is Frobenius; hence, $\\Gamma_G$ is copies of $K_{1}$ and so the maximum cocliques are known.\n", n);
            elif num_components > 1 then
                PrintToFormatted(latex_output, "As $\\Gamma_G$ has {1} components, the composition of the non-canonical cocliques is known.\n", num_components);
            else
                if chk = true then
                    PrintToFormatted(latex_output, "The maximum cocliques are split into {1} parts, each of size {2} containing {3} subgroups and {4} cosets. ", size_cl_partitions, Size(stabs_cosets), Size(cocliques_subg)/size_cl_partitions, Size(cocliques_cosets)/size_cl_partitions);
                elif chk = false then
                    PrintToFormatted(latex_output, "Besides the canonical cocliques, it is unclear what the {1} other maximum cocliques are. ", Size(unknown_solutions));
                fi;
            fi;
            
            PrintToFormatted(latex_output, "The number of maximum cocliques that are subgroups: {1}; cosets: {2}; neither: {3}. ", Size(cocliques_subg), Size(cocliques_cosets), Size(cocliques_other));
        
            if iso = true then
                PrintToFormatted(latex_output, "All maximum cocliques that are subgroups are isomorphic. The structure description of these subgroups: {1}. ", structure);
            else
                WriteAll(latex_output, "Not all maximum cocliques that are subgroups are isomorphic. ");
                PrintToFormatted(latex_output, "The number of isomorphic pairs of maximum cocliques that are subgroups: {1}. ", Size(subg_iso));
                PrintToFormatted(latex_output, "The number of non-isomorphic pairs of maximum cocliques that are subgroups: {1}. ", Size(subg_not_iso));
            fi;
            
            PrintToFormatted(latex_output, "The number of conjucacy classes among the maximum cocliques that are subgroups: {1}.\n", Size(conj_classes));
        fi;
    fi;
    
    #Question 3
    if 3 in extra_tests then
        if clique_type = 1 then
            dim_all_cliques := dim_clique;
        fi;
        
        if zeroes_check = 0 then
            PrintToFormatted(latex_output, "\\item $\\dim(\\mathcal{{C}}^\\prime) = {1}$; $\\dim(\\mathcal{{C}}) = {2}$; $\\dim(W) = {3}$.\n", dim_clique, dim_all_cliques, dim_required);
        else
            if zeroes_check = 1 then
                PrintToFormatted(latex_output, "\\item There is {1} module that contain neither a maximum clique nor a maximum coclique; $\\dim(\\mathcal{{C}}^\\prime) = {2}$; $\\dim(\\mathcal{{C}}) = {3}$; $\\dim(W) = {4}$.\n", zeroes_check, dim_clique, dim_all_cliques, dim_required);
            else
                PrintToFormatted(latex_output, "\\item There are {1} modules that contain neither a maximum clique nor a maximum coclique; $\\dim(\\mathcal{{C}}^\\prime) = {2}$; $\\dim(\\mathcal{{C}}) = {3}$; $\\dim(W) = {4}$.\n", zeroes_check, dim_clique, dim_all_cliques, dim_required);
            fi;
        fi;
    fi;
    
    #Question 4
    if 4 in extra_tests then
        if ratio = 1 then
            WriteAll(latex_output, "\\item The ratio bound holds with equality; ");
        elif ratio = 0 then
            WriteAll(latex_output, "\\item The ratio bound holds but not with equality; ");
        else
            WriteAll(latex_output, "\\item The ratio bound does not hold; ");
        fi;
        
        WriteAll(latex_output, "$\\sigma(\\Gamma_{G})=\\{$");
        for i in [1..Size(eMult)] do
            if i <> Size(eMult) then
                PrintToFormatted(latex_output, "${1}^{{({2})}}$, ", eMult[i][1], eMult[i][2]);
            else
                PrintToFormatted(latex_output, "${1}^{{({2})}}\\}}$.\n", eMult[i][1], eMult[i][2]);
            fi;
        od;
    fi;
    
    if Size(extra_tests) > 0 then
        WriteAll(latex_output, "\\end{enumerate}\n");
        CloseStream(latex_output);
    fi;
fi;