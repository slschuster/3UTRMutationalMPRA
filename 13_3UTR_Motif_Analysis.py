#! usr/bin/env/python3

from helper import make_tri_nt_context_dic, gen_random_muts, score
import random
import statistics
from datetime import datetime

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

# number of simulations to run
iters = 100

print("Doing " + str(iters) + " iterations...")

# open and read mutations file
mutations_file = open("all_3utr_seqs_01062021.txt", 'r')
mutations_lines = (mutations_file.read().splitlines())[1:]
mutations_file.close()
###################################################################################
# ___________________Collect mutation locations, wts, and muts_____________________
###################################################################################
mut_locs = []
wts = []
muts = []

# where count is 
index1 = 20
# where wt is 
index2 = 23
# where mut is 
index3 = 24

for i in range(len(mutations_lines)):
    line = mutations_lines[i].split("\t")
    
    # find mutation location
    counter = 0
    while line[index2][counter] == line[index3][counter]:
        counter += 1
    mut_locs.append(counter)

    wts.append(line[index2])
    muts.append(line[index3])
    
    # check that locations are right     
    if line[index2][counter] == line[index3][counter]:
        print("problem: didn't find mutation")
    if line[index2][counter-1] != line[index3][counter-1]:
        print("problem: other nt don't match")
print(str(len(mut_locs)) + " total mutations")

###############################################################################
# ___________________Make Trinucleotide Context Dictionary_____________________
###############################################################################
dic = make_tri_nt_context_dic()
OTHERS_LIST = []

for i in range(len(mutations_lines)): 
    mutant = muts[i]
    x = mut_locs[i]
    wildtype = wts[i]

    # count is number of patients to weight 3n context
    count = mutations_lines[i].split("\t")[index1]
    
    # dic of where mutations lie
    if ((x+2) <= len(wildtype) and ((x-1) >= 0)):
        mut_trip = wildtype[(x-1):(x+2)]
        new_base = mutant[x]
        if mut_trip in dic.keys(): 
            dic[mut_trip]["count"] = dic[mut_trip]["count"] + int(count)
            if new_base == "A":
                (dic[mut_trip]["to_a"]) = (dic[mut_trip]["to_a"]) + int(count)
            elif new_base == "C":
                (dic[mut_trip]["to_c"]) = (dic[mut_trip]["to_c"]) + int(count)
            elif new_base == "T":
                (dic[mut_trip]["to_t"]) = (dic[mut_trip]["to_t"]) + int(count)
            elif new_base == "G":
                (dic[mut_trip]["to_g"]) = (dic[mut_trip]["to_g"]) + int(count)
            else: 
                print("ERROR in making trinucleotide context dictionary")
        else: 
            OTHERS_LIST.append(mut_trip)
    elif ((x+2) > len(wildtype)): 
        OTHERS_LIST.append(wildtype[(x-1):(x+1)])
    else: 
        OTHERS_LIST.append(wildtype[x:(x+2)])

# calculate background 3nt freq

dic2 = make_tri_nt_context_dic()
for j in range(len(wts)): 
    count = int(mutations_lines[j].split("\t")[index1])
    for i in (range(len(wts[j]) - 2)):
        dic2[(wts[j])[i:i+3]]["count"] += count

for i in dic.keys(): 
    for j in dic[i]: 
        dic[i][j] = dic[i][j]/dic2[i]["count"]

mut_sets_list = []

# mutate randomly
for _ in range(iters):
    muts_list  = gen_random_muts(wts, dic)
    mut_sets_list.append(muts_list)

print(str((sum([len(f) for f in mut_sets_list]))/iters) + " unique mutations introduced on average in permutations")

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

##################################################################################
# ______________________How many mutations affect miRNA sites_____________________
##################################################################################
muts_that_impact_mirnas = 0
muts_that_add_mirnas = 0
muts_that_remove_mirnas = 0

mirna_changed_file = open("mirna_mutations.txt", "w+")

# load miRNA into a dictionary 
mirna_lines_f = open("motif_lists/miRNA_list.txt", "r")
mirna_lines = (mirna_lines_f.read().splitlines())
mirna_lines_f.close()
mirna_seqs = {}
for i in range(len(mirna_lines)//2):
    # note: many mirna sequences have duplicates, but we removed these in process_initial_motif_file
    # create dictionary of mrna sequences
    if mirna_lines[i*2+1] not in mirna_seqs: 
        mirna_seqs[mirna_lines[i*2+1]] = (mirna_lines[2*i])
    

# how many add, remove, both mirna things
for i in range(len(mutations_lines)-1):
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  
    count = int(mutations_lines[i].split("\t")[index1])
    # find start and end indexes
    if mut_loc < 7: 
        start = 0
    else: 
        start = mut_loc - 6
    if mut_loc > len(wt_seq) - 7: 
        end = len(wt_seq)
    else: 
        end = mut_loc + 8

    # iterate through 13 nt and check if any in mirna
    wt_mirna_list = []
    mut_mirna_list = []
    #print(" ")
    #print(wt_seq[start:end])
    for j in range (end-start-7):
        #print(wt_seq[start+j:end-6+j])
        #print(mut_seq[start+j:end-6+j])
        if wt_seq[start+j:end-7+j] in mirna_seqs: 
            wt_mirna_list.append(mirna_seqs[wt_seq[start+j:end-7+j]])
        if mut_seq[start+j:end-7+j] in mirna_seqs:
            mut_mirna_list.append(mirna_seqs[mut_seq[start+j:end-7+j]])
    if wt_mirna_list != mut_mirna_list:
        muts_that_impact_mirnas += count
    if len(wt_mirna_list) > len(mut_mirna_list):
        muts_that_add_mirnas += count
    elif len(mut_mirna_list) > len(wt_mirna_list):
        muts_that_remove_mirnas += count
    mirna_changed_file.write(str(mutations_lines[i].split('\t')[6]) + "\t" + str(mutations_lines[i].split('\t')[19]) + "\t" + str(wt_mirna_list) + "\t" + str(mut_mirna_list) + "\n")
print(str(muts_that_impact_mirnas) + " mutations affect mirnas")
print(str(muts_that_add_mirnas) + " mutations add mirnas")
print(str(muts_that_remove_mirnas) + " mutations remove mirnas")
####################################################################################
# ______________________Simulate random mutation of miRNA sites_____________________
####################################################################################
mutate = []
added = []
removed = []

for i in mut_sets_list:
    muts_that_remove_mirnas = 0
    muts_that_add_mirnas = 0
    muts_that_impact_mirnas = 0
    count = int(mutations_lines[j].split("\t")[index1])
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]
        if mut_loc < 7: 
            start = 0
        else: 
            start = mut_loc - 6
        if mut_loc > len(wt_seq) - 7: 
            end = len(wt_seq)
        else: 
            end = mut_loc + 8

        # iterate through 13 nt and check if any in mirna
        wt_mirna_list = []
        mut_mirna_list = []
        #print(" ")
        #print(wt_seq[start:end])
        for j in range (end-start-7):
            #print(wt_seq[start+j:end-6+j])
            #print(mut_seq[start+j:end-6+j])
            if wt_seq[start+j:end-7+j] in mirna_seqs: 
                wt_mirna_list.append(mirna_seqs[wt_seq[start+j:end-7+j]])
            if mut_seq[start+j:end-7+j] in mirna_seqs:
                mut_mirna_list.append(mirna_seqs[mut_seq[start+j:end-7+j]])
        if wt_mirna_list != mut_mirna_list:
            muts_that_impact_mirnas += count
        if len(wt_mirna_list) > len(mut_mirna_list):
            muts_that_add_mirnas += count
        elif len(mut_mirna_list) > len(wt_mirna_list):
            muts_that_remove_mirnas += count
    added.append(muts_that_add_mirnas)
    removed.append(muts_that_remove_mirnas)
    mutate.append(muts_that_impact_mirnas)
print(str(sum(mutate)/iters) + " mutations affect mirnas on average")
print(str(statistics.stdev(mutate))  + " is stdev of mutations that affect mirnas")
print(str(sum(added)/iters) + " mutations add mirnas on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add mirnas")
print(str(sum(removed)/iters) + " mutations remove mirnas on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove mirnas")

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

################################################################################
# ______________________How many mutations affect M6A sites_____________________
################################################################################

M6A_mutated = 0
M6A_added = 0
M6A_removed = 0
for i in range(len(mutations_lines)):
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    # find start and end indexes
    if mut_loc < 5: 
        start = 0
    else: 
        start = mut_loc - 4
    if mut_loc > len(wt_seq) - 5: 
        end = len(wt_seq)
    else: 
        end = mut_loc + 6

    inwt = False
    inmut = False

    if "GAACT" in wt_seq[start:end]: 
        inwt = True
    if "GAACT" in mut_seq[start:end]:
        inmut = True
    
    if inwt or inmut: 
        M6A_mutated += 1
    if inwt and (not inmut):
        M6A_removed += 1
    elif inmut and (not inwt):
        M6A_added += 1

print(str(M6A_mutated) + " mutations affect M6A")
print(str(M6A_added) + " mutations add M6A")
print(str(M6A_removed) + " mutations remove M6A")

######################################################################################################
# ______________________Permutation to determine statistics of mutating M6A sites_____________________
######################################################################################################
mutated = []
added = []
removed = []
for i in mut_sets_list:
    M6A_mutated = 0
    M6A_added = 0
    M6A_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]
        # find start and end indexes
        if mut_loc < 5: 
            start = 0
        else: 
            start = mut_loc - 4
        if mut_loc > len(wt_seq) - 5: 
            end = len(wt_seq)
        else: 
            end = mut_loc + 6

        inwt = False
        inmut = False

        if "GAACT" in wt_seq[start:end]: 
            inwt = True
        if "GAACT" in mut_seq[start:end]:
            inmut = True
        
        if inwt or inmut: 
            M6A_mutated += 1
        if inwt and (not inmut):
            M6A_removed += 1
        elif inmut and (not inwt):
            M6A_added += 1
    mutated.append(M6A_mutated)
    added.append(M6A_added)
    removed.append(M6A_removed)
print(str(sum(mutated)/iters) + " mutations affect M6A on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect M6A")
print(str(sum(added)/iters) + " mutations add M6A on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add M6A")
print(str(sum(removed)/iters) + " mutations remove M6A on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove M6A")

########################################################################################
# ______________________How many mutations affect 2' O Methylation?_____________________
########################################################################################
'''
mutated = 0
added = 0
removed = 0

for i in range(len(mutations_lines)):
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    # find start and end indexes
    if mut_loc < 10: 
        start = 0
    else: 
        start = mut_loc - 9
    if mut_loc > len(wt_seq) - 10: 
        end = len(wt_seq)
    else: 
        end = mut_loc + 6

    inwt = False
    inmut = False

    if "AGATC" in wt_seq[start:end]: 
        counter  = 0
        while wt_seq[start+counter:start+counter+5] != "AGATC":
            counter += 1
        #print(wt_seq[start+counter:start+counter+11])  -- there seem to be many that are close, should consider adjustment 
        # if there's a second agatc, use that as the starting point
        if "AGATC" in mut_seq[start+counter+5:end]: 
            counter += 1
            while mut_seq[start+counter:start+counter+5] != "AGATC":
                counter += 1
        five_purines = True
        for j in range(5):
            if start+counter+6+j < len(wt_seq):
                if wt_seq[start+counter+6+j] != "G" and wt_seq[start+counter+6+j] != "A":
                    five_purines = False
            else:
                five_purines = False
        if five_purines: 
            inwt = True 
    
    if "AGATC" in mut_seq[start:end]: 
        counter  = 0
        while mut_seq[start+counter:start+counter+5] != "AGATC":
            counter += 1
        # if there's a second agatc, use that as the starting point
        if "AGATC" in mut_seq[start+counter+5:end]: 
            counter += 1
            while mut_seq[start+counter:start+counter+5] != "AGATC":
                counter += 1
        five_purines = True
        for j in range(5):
            if start+counter+6+j < len(mut_seq):
                if mut_seq[start+counter+6+j] != "G" and mut_seq[start+counter+6+j] != "A":
                    five_purines = False
            else:
                five_purines = False
        if five_purines: 
            inwt = True 

    if inwt or inmut: 
        mutated += 1
    if inwt and (not inmut):
        removed += 1
    elif inmut and (not inwt):
        added += 1

print(str(mutated) + " mutations affect 2'O Methylation")
print(str(added) + " mutations add 2'O Methylation")
print(str(removed) + " mutations remove 2'O Methylation")

######################################################################################################
# ________________Permutation to determine statistics of mutating 2'O Methylation sites_______________
######################################################################################################

mutated = []
added = []
removed = []
for i in mut_sets_list:
    O2_mutated = 0
    O2_added = 0
    O2_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]

        # find start and end indexes
        if mut_loc < 10: 
            start = 0
        else: 
            start = mut_loc - 9
        if mut_loc > len(wt_seq) - 10: 
            end = len(wt_seq)
        else: 
            end = mut_loc + 6

        inwt = False
        inmut = False

        if "AGATC" in wt_seq[start:end]: 
            counter  = 0
            while wt_seq[start+counter:start+counter+5] != "AGATC":
                counter += 1
            #print(wt_seq[start+counter:start+counter+11])  -- there seem to be many that are close, should consider adjustment 
            # if there's a second agatc, use that as the starting point
            if "AGATC" in mut_seq[start+counter+5:end]: 
                counter += 1
                while mut_seq[start+counter:start+counter+5] != "AGATC":
                    counter += 1
            five_purines = True
            for j in range(5):
                if start+counter+6+j < len(wt_seq):
                    if wt_seq[start+counter+6+j] != "G" and wt_seq[start+counter+6+j] != "A":
                        five_purines = False
                else:
                    five_purines = False
            if five_purines: 
                inwt = True 
        
        if "AGATC" in mut_seq[start:end]: 
            counter  = 0
            while mut_seq[start+counter:start+counter+5] != "AGATC":
                counter += 1
            # if there's a second agatc, use that as the starting point
            if "AGATC" in mut_seq[start+counter+5:end]: 
                counter += 1
                while mut_seq[start+counter:start+counter+5] != "AGATC":
                    counter += 1
            five_purines = True
            for j in range(5):
                if start+counter+6+j < len(mut_seq):
                    if mut_seq[start+counter+6+j] != "G" and mut_seq[start+counter+6+j] != "A":
                        five_purines = False
                else:
                    five_purines = False
            if five_purines: 
                inwt = True 

        if inwt or inmut: 
            O2_mutated += 1
        if inwt and (not inmut):
            O2_removed += 1
        elif inmut and (not inwt):
            O2_added += 1
    mutated.append(O2_mutated)
    added.append(O2_added)
    removed.append(O2_removed)
print(str(sum(mutated)/iters) + " mutations affect 2'O Methylation on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect 2'O Methylation")
print(str(sum(added)/iters) + " mutations add 2'O Methylation on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add 2'O Methylation")
print(str(sum(removed)/iters) + " mutations remove 2'O Methylation on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove 2'O Methylation")

'''
####################################################################################################
# ________________use CISBP motifs to determine how many of our mutations affect them_______________
####################################################################################################

mutated = 0
added = 0
removed = 0

cisbp_changed_file = open("cisbp_mutations.txt", "w+")

motif_file = open("motif_lists/cisbp_motif_list_wmax.txt", "r")
motif_lines = motif_file.read().splitlines()
motif_file.close()
percent = 0
for i in range(len(mutations_lines)):
    if ((i+8)%600 == 0):
        #print(str(percent) + " percent done with CISBP native analysis")
        percent = percent + 10
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    inwt = []
    inmut = []
    for j in range(len(motif_lines)//6):
        if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inwt.append(motif_lines[j*6])
        if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inmut.append(motif_lines[j*6])
    
    something_added = 0
    something_removed = 0
    something = 0
    for j in inwt:
        something = 1 
        if j not in inmut: 
            something_added = 1
    for j in inmut: 
        something = 1
        if j not in inwt: 
            something_removed = 1
    added += something_added
    removed += something_removed
    mutated += something
    cisbp_changed_file.write(str(mutations_lines[i].split('\t')[6]) + "\t" + str(mutations_lines[i].split('\t')[19]) + "\t" + str(inwt) + "\t" + str(inmut) + "\n")

print(str(mutated) + " mutations affect Direct CISBP Motifs")
print(str(added) + " mutations add Direct CISBP Motifs")
print(str(removed) + " mutations remove Direct CISBP Motifs")
cisbp_changed_file.close()

############################################################################################
# ________________Permutation to determine statistics of mutating CISBP sites_______________
############################################################################################   


mutated = []
added = []
removed = []
iteration = 0

for i in mut_sets_list:
    print("on the " + str(iteration) + "th iteration of simulating cisbp mutations")
    iteration += 1
    cisbp_mutated = 0
    cisbp_added = 0
    cisbp_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]

        inwt = []
        inmut = []
        for j in range(len(motif_lines)//6):
            if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inwt.append(motif_lines[j*6])
            if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inmut.append(motif_lines[j*6])
        something_added = 0
        something_removed = 0
        something = 0
        for j in inwt:
            something = 1 
            if j not in inmut: 
                something_added = 1
        for j in inmut: 
            something = 1
            if j not in inwt: 
                something_removed = 1
        cisbp_added += something_added
        cisbp_removed += something_removed
        cisbp_mutated += something

    mutated.append(cisbp_mutated)
    added.append(cisbp_added)
    removed.append(cisbp_removed)

print(str(sum(mutated)/iters) + " mutations affect CISBP on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect CISBP")
print(str(sum(added)/iters) + " mutations add CISBP on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add CISBP")
print(str(sum(removed)/iters) + " mutations remove CISBP on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove CISBP")





now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)




####################################################################################################
# ________________use eCLIP motifs to determine how many of our mutations affect them_______________
####################################################################################################

mutated = 0
added = 0
removed = 0

eclip_changed_file = open("eclip_mutations.txt", "w+")

motif_file = open("motif_lists/eclip_motif_list_wmax.txt", "r")
motif_lines = motif_file.read().splitlines()
motif_file.close()
percent = 0
for i in range(len(mutations_lines)):
    if ((i+8)%600 == 0):
        #print(str(percent) + " percent done with CISBP native analysis")
        percent = percent + 10
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    inwt = []
    inmut = []
    for j in range(len(motif_lines)//6):
        if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inwt.append(motif_lines[j*6])
        if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inmut.append(motif_lines[j*6])
    
    something_added = 0
    something_removed = 0
    something = 0
    for j in inwt:
        something = 1 
        if j not in inmut: 
            something_added = 1
    for j in inmut: 
        something = 1
        if j not in inwt: 
            something_removed = 1
    added += something_added
    removed += something_removed
    mutated += something
    eclip_changed_file.write(str(mutations_lines[i].split('\t')[6]) + "\t" + str(mutations_lines[i].split('\t')[19]) + "\t" + str(inwt) + "\t" + str(inmut) + "\n")


print(str(mutated) + " mutations affect Direct eCLIP Motifs")
print(str(added) + " mutations add Direct eCLIP Motifs")
print(str(removed) + " mutations remove Direct eCLIP Motifs")
eclip_changed_file.close()

############################################################################################
# ________________Permutation to determine statistics of mutating eCLIP sites_______________
############################################################################################   


mutated = []
added = []
removed = []
iteration = 0

for i in mut_sets_list:
    print("on the " + str(iteration) + "th iteration of simulating eCLIP mutations")
    iteration += 1
    eCLIP_mutated = 0
    eCLIP_added = 0
    eCLIP_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]

        inwt = []
        inmut = []
        for j in range(len(motif_lines)//6):
            if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inwt.append(motif_lines[j*6])
            if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inmut.append(motif_lines[j*6])
        something_added = 0
        something_removed = 0
        something = 0
        for j in inwt:
            something = 1 
            if j not in inmut: 
                something_added = 1
        for j in inmut: 
            something = 1
            if j not in inwt: 
                something_removed = 1
        eCLIP_added += something_added
        eCLIP_removed += something_removed
        eCLIP_mutated += something

    mutated.append(eCLIP_mutated)
    added.append(eCLIP_added)
    removed.append(eCLIP_removed)

print(str(sum(mutated)/iters) + " mutations affect eCLIP on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect eCLIP")
print(str(sum(added)/iters) + " mutations add eCLIP on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add eCLIP")
print(str(sum(removed)/iters) + " mutations remove eCLIP on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove eCLIP")






now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)






####################################################################################################
# ________________use RBN Berger motifs to determine how many of our mutations affect them_______________
####################################################################################################

mutated = 0
added = 0
removed = 0

rbn_changed_file = open("RBN_mutations.txt", "w+")

motif_file = open("motif_lists/RBN_Berger_list_wmax.txt", "r")
motif_lines = motif_file.read().splitlines()
motif_file.close()
percent = 0
for i in range(len(mutations_lines)-1):
    if ((i+8)%600 == 0):
        #print(str(percent) + " percent done with CISBP native analysis")
        percent = percent + 10
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    inwt = []
    inmut = []
    for j in range(len(motif_lines)//6):
        if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inwt.append(motif_lines[j*6])
        if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inmut.append(motif_lines[j*6])
    
    something_added = 0
    something_removed = 0
    something = 0
    for j in inwt:
        something = 1 
        if j not in inmut: 
            something_added = 1
    for j in inmut: 
        something = 1
        if j not in inwt: 
            something_removed = 1
    added += something_added
    removed += something_removed
    mutated += something
    rbn_changed_file.write(str(mutations_lines[i].split('\t')[6]) + "\t" + str(mutations_lines[i].split('\t')[19]) + "\t" + str(inwt) + "\t" + str(inmut) + "\n")


print(str(mutated) + " mutations affect Direct RBN Motifs")
print(str(added) + " mutations add Direct RBN Motifs")
print(str(removed) + " mutations remove Direct RBN Motifs")

rbn_changed_file.close()

############################################################################################
# ________________Permutation to determine statistics of mutating RBN sites_______________
############################################################################################   


mutated = []
added = []
removed = []
iteration = 0

for i in mut_sets_list:
    print("on the " + str(iteration) + "th iteration of simulating RBN mutations")
    iteration += 1
    RBN_mutated = 0
    RBN_added = 0
    RBN_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]

        inwt = []
        inmut = []
        for j in range(len(motif_lines)//6):
            if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inwt.append(motif_lines[j*6])
            if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inmut.append(motif_lines[j*6])
        something_added = 0
        something_removed = 0
        something = 0
        for j in inwt:
            something = 1 
            if j not in inmut: 
                something_added = 1
        for j in inmut: 
            something = 1
            if j not in inwt: 
                something_removed = 1
        RBN_added += something_added
        RBN_removed += something_removed
        RBN_mutated += something

    mutated.append(RBN_mutated)
    added.append(RBN_added)
    removed.append(RBN_removed)

print(str(sum(mutated)/iters) + " mutations affect RBN on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect RBN")
print(str(sum(added)/iters) + " mutations add RBN on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add RBN")
print(str(sum(removed)/iters) + " mutations remove RBN on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove RBN")


####################################################################################################
# ________________use PAS PWM to determine how many of our mutations affect them____________________
####################################################################################################

mutated = 0
added = 0
removed = 0

pas_changed_file = open("pas_mutations.txt", "w+")

motif_file = open("motif_lists/PAS_PWM_wmax.txt", "r")
motif_lines = motif_file.read().splitlines()
motif_file.close()
percent = 0
for i in range(len(mutations_lines)):
    if ((i+8)%600 == 0):
        #print(str(percent) + " percent done with CISBP native analysis")
        percent = percent + 10
    wt_seq = wts[i]
    mut_seq = muts[i]
    mut_loc = mut_locs[i]  

    inwt = []
    inmut = []
    for j in range(len(motif_lines)//6):
        if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inwt.append(motif_lines[j*6])
        if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
            inmut.append(motif_lines[j*6])
    
    something_added = 0
    something_removed = 0
    something = 0
    for j in inwt:
        something = 1 
        if j not in inmut: 
            something_added = 1
    for j in inmut: 
        something = 1
        if j not in inwt: 
            something_removed = 1
    added += something_added
    removed += something_removed
    mutated += something
    pas_changed_file.write(str(mutations_lines[i].split('\t')[6]) + "\t" + str(mutations_lines[i].split('\t')[19]) + "\t" + str(inwt) + "\t" + str(inmut) + "\n")

print(str(mutated) + " mutations affect PAS")
print(str(added) + " mutations add PAS")
print(str(removed) + " mutations remove PAS")
pas_changed_file.close()

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

############################################################################################
# ________________Permutation to determine statistics of mutating PAS sites_________________
############################################################################################   


mutated = []
added = []
removed = []
iteration = 0

for i in mut_sets_list:
    iteration += 1
    cisbp_mutated = 0
    cisbp_added = 0
    cisbp_removed = 0
    for j in i:
        wt_seq = j[0]
        mut_seq = j[1]
        mut_loc = j[2]

        inwt = []
        inmut = []
        for j in range(len(motif_lines)//6):
            if score(wt_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inwt.append(motif_lines[j*6])
            if score(mut_seq, mut_loc, motif_lines[1+j*6:6+j*6]):
                inmut.append(motif_lines[j*6])
        something_added = 0
        something_removed = 0
        something = 0
        for j in inwt:
            something = 1 
            if j not in inmut: 
                something_added = 1
        for j in inmut: 
            something = 1
            if j not in inwt: 
                something_removed = 1
        cisbp_added += something_added
        cisbp_removed += something_removed
        cisbp_mutated += something

    mutated.append(cisbp_mutated)
    added.append(cisbp_added)
    removed.append(cisbp_removed)

print(str(sum(mutated)/iters) + " mutations affect PAS on average")
print(str(statistics.stdev(mutated))  + " is stdev of mutations that affect PAS")
print(str(sum(added)/iters) + " mutations add PAS on average")
print(str(statistics.stdev(added))  + " is stdev of mutations that add PAS")
print(str(sum(removed)/iters) + " mutations remove PAS on average")
print(str(statistics.stdev(removed))  + " is stdev of mutations that remove PAS")




now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
print("job finished")