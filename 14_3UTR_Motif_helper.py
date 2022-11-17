#! usr/bin/env/python3 
import random
def make_tri_nt_context_dic():

    dic = {}
    #STARTING WITH A
    dic_AAA = {}
    dic_AAA['count'] = 0
    dic_AAA["to_t"] = 0
    dic_AAA["to_c"] = 0
    dic_AAA["to_g"] = 0
    dic['AAA'] = dic_AAA
    dic_AAC = {}
    dic_AAC['count'] = 0
    dic_AAC["to_t"] = 0
    dic_AAC["to_c"] = 0
    dic_AAC["to_g"] = 0
    dic['AAC'] = dic_AAC
    dic_AAG = {}
    dic_AAG['count'] = 0
    dic_AAG["to_t"] = 0
    dic_AAG["to_c"] = 0
    dic_AAG["to_g"] = 0
    dic['AAG'] = dic_AAG
    dic_AAT = {}
    dic_AAT['count'] = 0
    dic_AAT["to_t"] = 0
    dic_AAT["to_c"] = 0
    dic_AAT["to_g"] = 0
    dic['AAT'] = dic_AAT

    dic_ACA = {}
    dic_ACA['count'] = 0
    dic_ACA["to_t"] = 0
    dic_ACA["to_a"] = 0
    dic_ACA["to_g"] = 0
    dic['ACA'] = dic_ACA
    dic_ACC = {}
    dic_ACC['count'] = 0
    dic_ACC["to_t"] = 0
    dic_ACC["to_a"] = 0
    dic_ACC["to_g"] = 0
    dic['ACC'] = dic_ACC
    dic_ACG = {}
    dic_ACG['count'] = 0
    dic_ACG["to_t"] = 0
    dic_ACG["to_a"] = 0
    dic_ACG["to_g"] = 0
    dic['ACG'] = dic_ACG
    dic_ACT = {}
    dic_ACT['count'] = 0
    dic_ACT["to_t"] = 0
    dic_ACT["to_a"] = 0
    dic_ACT["to_g"] = 0
    dic['ACT'] = dic_ACT

    dic_AGA = {}
    dic_AGA['count'] = 0
    dic_AGA["to_t"] = 0
    dic_AGA["to_a"] = 0
    dic_AGA["to_c"] = 0
    dic['AGA'] = dic_AGA
    dic_AGC = {}
    dic_AGC['count'] = 0
    dic_AGC["to_t"] = 0
    dic_AGC["to_a"] = 0
    dic_AGC["to_c"] = 0
    dic['AGC'] = dic_AGC
    dic_AGG = {}
    dic_AGG['count'] = 0
    dic_AGG["to_t"] = 0
    dic_AGG["to_a"] = 0
    dic_AGG["to_c"] = 0
    dic['AGG'] = dic_AGG
    dic_AGT = {}
    dic_AGT['count'] = 0
    dic_AGT["to_t"] = 0
    dic_AGT["to_a"] = 0
    dic_AGT["to_c"] = 0
    dic['AGT'] = dic_AGT

    dic_ATA = {}
    dic_ATA['count'] = 0
    dic_ATA["to_g"] = 0
    dic_ATA["to_a"] = 0
    dic_ATA["to_c"] = 0
    dic['ATA'] = dic_ATA
    dic_ATC = {}
    dic_ATC['count'] = 0
    dic_ATC["to_g"] = 0
    dic_ATC["to_a"] = 0
    dic_ATC["to_c"] = 0
    dic['ATC'] = dic_ATC
    dic_ATG = {}
    dic_ATG['count'] = 0
    dic_ATG["to_g"] = 0
    dic_ATG["to_a"] = 0
    dic_ATG["to_c"] = 0
    dic['ATG'] = dic_ATG
    dic_ATT = {}
    dic_ATT['count'] = 0
    dic_ATT["to_g"] = 0
    dic_ATT["to_a"] = 0
    dic_ATT["to_c"] = 0
    dic['ATT'] = dic_ATT

    ##STARGING WITH C

    dic_CAA = {}
    dic_CAA['count'] = 0
    dic_CAA["to_t"] = 0
    dic_CAA["to_g"] = 0
    dic_CAA["to_c"] = 0
    dic['CAA'] = dic_CAA
    dic_CAC = {}
    dic_CAC['count'] = 0
    dic_CAC["to_t"] = 0
    dic_CAC["to_g"] = 0
    dic_CAC["to_c"] = 0
    dic['CAC'] = dic_CAC
    dic_CAG = {}
    dic_CAG['count'] = 0
    dic_CAG["to_t"] = 0
    dic_CAG["to_g"] = 0
    dic_CAG["to_c"] = 0
    dic['CAG'] = dic_CAG
    dic_CAT = {}
    dic_CAT['count'] = 0
    dic_CAT["to_t"] = 0
    dic_CAT["to_g"] = 0
    dic_CAT["to_c"] = 0
    dic['CAT'] = dic_CAT

    dic_CCA = {}
    dic_CCA['count'] = 0
    dic_CCA["to_t"] = 0
    dic_CCA["to_g"] = 0
    dic_CCA["to_a"] = 0
    dic['CCA'] = dic_CCA
    dic_CCC = {}
    dic_CCC['count'] = 0
    dic_CCC["to_t"] = 0
    dic_CCC["to_g"] = 0
    dic_CCC["to_a"] = 0
    dic['CCC'] = dic_CCC
    dic_CCG = {}
    dic_CCG['count'] = 0
    dic_CCG["to_t"] = 0
    dic_CCG["to_g"] = 0
    dic_CCG["to_a"] = 0
    dic['CCG'] = dic_CCG
    dic_CCT = {}
    dic_CCT['count'] = 0
    dic_CCT["to_t"] = 0
    dic_CCT["to_g"] = 0
    dic_CCT["to_a"] = 0
    dic['CCT'] = dic_CCT

    dic_CGA = {}
    dic_CGA['count'] = 0
    dic_CGA["to_t"] = 0
    dic_CGA["to_c"] = 0
    dic_CGA["to_a"] = 0
    dic['CGA'] = dic_CGA
    dic_CGC = {}
    dic_CGC['count'] = 0
    dic_CGC["to_t"] = 0
    dic_CGC["to_c"] = 0
    dic_CGC["to_a"] = 0
    dic['CGC'] = dic_CGC
    dic_CGG = {}
    dic_CGG['count'] = 0
    dic_CGG["to_t"] = 0
    dic_CGG["to_c"] = 0
    dic_CGG["to_a"] = 0
    dic['CGG'] = dic_CGG
    dic_CGT = {}
    dic_CGT['count'] = 0
    dic_CGT["to_t"] = 0
    dic_CGT["to_c"] = 0
    dic_CGT["to_a"] = 0
    dic['CGT'] = dic_CGT

    dic_CTA = {}
    dic_CTA['count'] = 0
    dic_CTA["to_g"] = 0
    dic_CTA["to_c"] = 0
    dic_CTA["to_a"] = 0
    dic['CTA'] = dic_CTA
    dic_CTC = {}
    dic_CTC['count'] = 0
    dic_CTC["to_g"] = 0
    dic_CTC["to_c"] = 0
    dic_CTC["to_a"] = 0
    dic['CTC'] = dic_CTC
    dic_CTG = {}
    dic_CTG['count'] = 0
    dic_CTG["to_g"] = 0
    dic_CTG["to_c"] = 0
    dic_CTG["to_a"] = 0
    dic['CTG'] = dic_CTG
    dic_CTT = {}
    dic_CTT['count'] = 0
    dic_CTT["to_g"] = 0
    dic_CTT["to_c"] = 0
    dic_CTT["to_a"] = 0
    dic['CTT'] = dic_CTT

    #STARTING WITH G
    dic_GAA = {}
    dic_GAA['count'] = 0
    dic_GAA["to_g"] = 0
    dic_GAA["to_c"] = 0
    dic_GAA["to_t"] = 0
    dic['GAA'] = dic_GAA
    dic_GAC = {}
    dic_GAC['count'] = 0
    dic_GAC["to_g"] = 0
    dic_GAC["to_c"] = 0
    dic_GAC["to_t"] = 0
    dic['GAC'] = dic_GAC
    dic_GAG = {}
    dic_GAG['count'] = 0
    dic_GAG["to_g"] = 0
    dic_GAG["to_c"] = 0
    dic_GAG["to_t"] = 0
    dic['GAG'] = dic_GAG
    dic_GAT = {}
    dic_GAT['count'] = 0
    dic_GAT["to_g"] = 0
    dic_GAT["to_c"] = 0
    dic_GAT["to_t"] = 0
    dic['GAT'] = dic_GAT

    dic_GCA = {}
    dic_GCA['count'] = 0
    dic_GCA["to_g"] = 0
    dic_GCA["to_a"] = 0
    dic_GCA["to_t"] = 0
    dic['GCA'] = dic_GCA
    dic_GCC = {}
    dic_GCC['count'] = 0
    dic_GCC["to_g"] = 0
    dic_GCC["to_a"] = 0
    dic_GCC["to_t"] = 0
    dic['GCC'] = dic_GCC
    dic_GCG = {}
    dic_GCG['count'] = 0
    dic_GCG["to_g"] = 0
    dic_GCG["to_a"] = 0
    dic_GCG["to_t"] = 0
    dic['GCG'] = dic_GCG
    dic_GCT = {}
    dic_GCT['count'] = 0
    dic_GCT["to_g"] = 0
    dic_GCT["to_a"] = 0
    dic_GCT["to_t"] = 0
    dic['GCT'] = dic_GCT

    dic_GGA = {}
    dic_GGA['count'] = 0
    dic_GGA["to_c"] = 0
    dic_GGA["to_a"] = 0
    dic_GGA["to_t"] = 0
    dic['GGA'] = dic_GGA
    dic_GGC = {}
    dic_GGC['count'] = 0
    dic_GGC["to_c"] = 0
    dic_GGC["to_a"] = 0
    dic_GGC["to_t"] = 0
    dic['GGC'] = dic_GGC
    dic_GGG = {}
    dic_GGG['count'] = 0
    dic_GGG["to_c"] = 0
    dic_GGG["to_a"] = 0
    dic_GGG["to_t"] = 0
    dic['GGG'] = dic_GGG
    dic_GGT = {}
    dic_GGT['count'] = 0
    dic_GGT["to_c"] = 0
    dic_GGT["to_a"] = 0
    dic_GGT["to_t"] = 0
    dic['GGT'] = dic_GGT

    dic_GTA = {}
    dic_GTA['count'] = 0
    dic_GTA["to_c"] = 0
    dic_GTA["to_a"] = 0
    dic_GTA["to_g"] = 0
    dic['GTA'] = dic_GTA
    dic_GTC = {}
    dic_GTC['count'] = 0
    dic_GTC["to_c"] = 0
    dic_GTC["to_a"] = 0
    dic_GTC["to_g"] = 0
    dic['GTC'] = dic_GTC
    dic_GTG = {}
    dic_GTG['count'] = 0
    dic_GTG["to_c"] = 0
    dic_GTG["to_a"] = 0
    dic_GTG["to_g"] = 0
    dic['GTG'] = dic_GTG
    dic_GTT = {}
    dic_GTT['count'] = 0
    dic_GTT["to_c"] = 0
    dic_GTT["to_a"] = 0
    dic_GTT["to_g"] = 0
    dic['GTT'] = dic_GTT

    ##STARGING WITH T

    dic_TAA = {}
    dic_TAA['count'] = 0
    dic_TAA["to_c"] = 0
    dic_TAA["to_t"] = 0
    dic_TAA["to_g"] = 0
    dic['TAA'] = dic_TAA
    dic_TAC = {}
    dic_TAC['count'] = 0
    dic_TAC["to_c"] = 0
    dic_TAC["to_t"] = 0
    dic_TAC["to_g"] = 0
    dic['TAC'] = dic_TAC
    dic_TAG = {}
    dic_TAG['count'] = 0
    dic_TAG["to_c"] = 0
    dic_TAG["to_t"] = 0
    dic_TAG["to_g"] = 0
    dic['TAG'] = dic_TAG
    dic_TAT = {}
    dic_TAT['count'] = 0
    dic_TAT["to_c"] = 0
    dic_TAT["to_t"] = 0
    dic_TAT["to_g"] = 0
    dic['TAT'] = dic_TAT

    dic_TCA = {}
    dic_TCA['count'] = 0
    dic_TCA["to_a"] = 0
    dic_TCA["to_t"] = 0
    dic_TCA["to_g"] = 0
    dic['TCA'] = dic_TCA
    dic_TCC = {}
    dic_TCC['count'] = 0
    dic_TCC["to_a"] = 0
    dic_TCC["to_t"] = 0
    dic_TCC["to_g"] = 0
    dic['TCC'] = dic_TCC
    dic_TCG = {}
    dic_TCG['count'] = 0
    dic_TCG["to_a"] = 0
    dic_TCG["to_t"] = 0
    dic_TCG["to_g"] = 0
    dic['TCG'] = dic_TCG
    dic_TCT = {}
    dic_TCT['count'] = 0
    dic_TCT["to_a"] = 0
    dic_TCT["to_t"] = 0
    dic_TCT["to_g"] = 0
    dic['TCT'] = dic_TCT

    dic_TGA = {}
    dic_TGA['count'] = 0
    dic_TGA["to_a"] = 0
    dic_TGA["to_t"] = 0
    dic_TGA["to_c"] = 0
    dic['TGA'] = dic_TGA
    dic_TGC = {}
    dic_TGC['count'] = 0
    dic_TGC["to_a"] = 0
    dic_TGC["to_t"] = 0
    dic_TGC["to_c"] = 0
    dic['TGC'] = dic_TGC
    dic_TGG = {}
    dic_TGG['count'] = 0
    dic_TGG["to_a"] = 0
    dic_TGG["to_t"] = 0
    dic_TGG["to_c"] = 0
    dic['TGG'] = dic_TGG
    dic_TGT = {}
    dic_TGT['count'] = 0
    dic_TGT["to_a"] = 0
    dic_TGT["to_t"] = 0
    dic_TGT["to_c"] = 0
    dic['TGT'] = dic_TGT

    dic_TTA = {}
    dic_TTA['count'] = 0
    dic_TTA["to_a"] = 0
    dic_TTA["to_g"] = 0
    dic_TTA["to_c"] = 0
    dic['TTA'] = dic_TTA
    dic_TTC = {}
    dic_TTC['count'] = 0
    dic_TTC["to_a"] = 0
    dic_TTC["to_g"] = 0
    dic_TTC["to_c"] = 0
    dic['TTC'] = dic_TTC
    dic_TTG = {}
    dic_TTG['count'] = 0
    dic_TTG["to_a"] = 0
    dic_TTG["to_g"] = 0
    dic_TTG["to_c"] = 0
    dic['TTG'] = dic_TTG
    dic_TTT = {}
    dic_TTT['count'] = 0
    dic_TTT["to_a"] = 0
    dic_TTT["to_g"] = 0
    dic_TTT["to_c"] = 0
    dic['TTT'] = dic_TTT
    return dic 


def gen_random_muts(wts, dic):
    muts_list = []
    for i in wts: 
        for j in (range(len(i)-2)):
            if random.random()<dic[i[j:j+3]]["count"]:
                # to select what nucleotide to mutate to
                total_norm = 0
                for k in [f for f in (dic[i[j:j+3]].keys()) if f != "count" ]:
                    total_norm += dic[i[j:j+3]][k]
                randvar = random.random()*total_norm
                sumvar = 0
                for k in [f for f in (dic[i[j:j+3]].keys()) if f != "count" ]:
                    sumvar += dic[i[j:j+3]][k]
                    if sumvar>=randvar: 
                        new_nt = k[-1].upper()
                muts_list.append((i, (i[0:j] + new_nt + i[j+1:len(i)]),j))
    return muts_list

# greater than 0 counts 
def score(seq, mut_loc, matrix):

    dic = {"A": 1, "C": 2, "G": 3, "T": 4}
    length = len(matrix[1].split())
    seq_length = len(seq)
    barrier = float(matrix[0]) * .9

    # for the length of the motif -- every possible starting position
    for i in range(length):
        # for the length of the motif -- score
        temp_score = 0
        for j in range(length):
            if mut_loc-i+j >= seq_length:
                temp_score -= 10000
            else:
                temp_score += float((matrix[dic[seq[mut_loc-i+j]]].split())[j])
        if temp_score >= barrier: 
            return True 
    return False

def max_score(matrix):
    sc = 0
    for i in range(len(matrix[0].split())):
        sc += max(float((matrix[0].split())[i]),float((matrix[1].split())[i]),float((matrix[2].split())[i]),float((matrix[3].split())[i]))
    return sc

def revcomp(mirna): 
    rc = ""
    for i in range(len(mirna)):
        if mirna[-i] == "U":
            rc+=("A")
        elif mirna[-i] == "C":
            rc+=("G")
        elif mirna[-i] == "G":
            rc+=("C")
        elif mirna[-i] == "A":
            rc+=("T")
        else: 
            print("error, bad base in revcomp")
    return rc
        