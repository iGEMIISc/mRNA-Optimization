import numpy as np
from arnie.bpps import bpps
from arnie.mfe import mfe
from arnie.utils import get_expected_accuracy, get_mean_base_pair_propensity
#from RiboGraphViz import RGV
from scipy.stats import gmean
import time
from DegScore import DegScore


# C2W was created using the values from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
# Maps codon -> weight = frequency of codon / max frequency of synonymous codons
# Only for calculating CAI w.r.t. human codon frequencies
C2W = {'UUU': 0.85185, 'UUC': 1.0, 'UUA': 0.2, 'UUG': 0.325, 'UCU': 0.79167, 
       'UCC': 0.91667, 'UCA': 0.625, 'UCG': 0.20833, 'UAU': 0.78571, 'UAC': 1.0, 
       'UAA': 0.6383, 'UAG': 0.51064, 'UGU': 0.85185, 'UGC': 1.0, 'UGA': 1.0, 
       'UGG': 1.0, 'CUU': 0.325, 'CUC': 0.5, 'CUA': 0.175, 'CUG': 1.0, 
       'CCU': 0.90625, 'CCC': 1.0, 'CCA': 0.875, 'CCG': 0.34375, 'CAU': 0.72414, 
       'CAC': 1.0, 'CAA': 0.36986, 'CAG': 1.0, 'CGU': 0.38095, 'CGC': 0.85714, 
       'CGA': 0.52381, 'CGG': 0.95238, 'AUU': 0.76596, 'AUC': 1.0, 'AUA': 0.3617, 
       'AUG': 1.0, 'ACU': 0.69444, 'ACC': 1.0, 'ACA': 0.77778, 'ACG': 0.30556, 
       'AAU': 0.88679, 'AAC': 1.0, 'AAA': 0.75439, 'AAG': 1.0, 'AGU': 0.625, 
       'AGC': 1.0, 'AGA': 1.0, 'AGG': 1.0, 'GUU': 0.3913, 'GUC': 0.52174, 
       'GUA': 0.26087, 'GUG': 1.0, 'GCU': 0.675, 'GCC': 1.0, 'GCA': 0.575, 
       'GCG': 0.275, 'GAU': 0.85185, 'GAC': 1.0, 'GAA': 0.72414, 'GAG': 1.0, 
       'GGU': 0.47059, 'GGC': 1.0, 'GGA': 0.73529, 'GGG': 0.73529}

UTR_list = [["GGGACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACC",
"UUCUAGAGCGGCCGCUUCGAGCCGGUUGAAUCGCUGAUCUCACGCCGUGGUGAGCUCGCUUUCUUGCUGUCCAAUUUCUAUUAAAGGUUCCUUUGUUCCCUAAGUCCAACUACUAAACUGGGGGAUAUUAUGAAGGGCCUUGAGCAUCUGGAUUCUGCCUAAUAAAAAACAUUUAUUUUCAUUGCAAAGUUCCGCGUACGUACGGCGUC"],
["GGGUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCUUGCUUAGUGCACUCACGCAGUAUAAUUAAUAACUAAUUACUGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGCACAUCUAGGUUUCGUCCGGGUGUGACCGAAAGGUAAGUUGGAGAGCCUUGUCCCUGGUUUCAACGAGAAAAC",
"UUCUAGAGCGGCCGCUUCGAGCCGAGACAAUCGCUGAUCUCACGCCGUGGUGAAAAGCAAAACUAACAUGAAACAAGGCUAGAAGUCAGGUCGGAUUAAGCCAUAGUACGGAAAAAACUAUGCUACCUGUGAGCCCCGUCCAAGGACGUUAAAAGAAGUCAGGCCAUCAUAAAUGCCAUAGCUUGAGUAAACUAUGCAGCCUGUAGCUCCACCUGAGAAGGUGUAAAAAAUCCGGGAGGCCACAAACCAUGGAAGCUGUACGCAUGGCGUAGUGGACUAGCGGUUAGAGGAGACCCCUCCCUUACAAAUCGCAGCAACAAUGGGGGCCCAAGGCGAGAUGAAGCUGUAGUCUCGCUGGAAGGACUAGAGGUUAGAGGAGACCCCCCCGAAACAAAAAACAGCAUAUUGACGCUGGGAAAGACCAGAGAUCCUGCUGUCUCCUCAGCAUCAUUCCAGGCACAGAACGCCAGAAAAUGGAAUGGUGCUGUUGAAUCAACAGGUUCUAGUUCCGCGUACGUACGGCGUC"],
["GGGACUCCUCCCCAUCCUCUCCCUCUGUCCCUCUGUCCCUCUGACCCUGCACUGUCCCAGCACC",
"UUCUAGAGCGGCCGCUUCGAGCUGCGACAAUCGCUGAUCUCACGCCGUGGUGACGUCUGCAUAACUUUUAUUAUUUCUUUUAUUAAUCAACAAAAUUUUGUUUUUAACAUUUCAGUUCCGCGUACGUACGGCGUC"]]

def calc_cai(seq):
    
    # Fixes sequence format. Convert to RNA
    seq = seq.upper().replace('T','U')
    if len(seq) % 3 != 0:
        raise ValueError("Not a valid coding sequence. Length is not a multiple of 3.")
    
    # Gets all the weights per codon
    w_list = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        # Do not count W or M codon since there is only one that encodes them
        if codon not in ['UGG', 'AUG']:
            w_list.append(C2W[codon])
            
    # return the geometric mean
    return gmean(w_list)


def get_MFE_open(seq, feature, package, param_file, linearfold, N=50):
    for i in range(0,N):
        try:
            ss = mfe(seq, package=package, viterbi=True, param_file=param_file, linear=linearfold)

            if "open" in feature:
                v = (14 - ss[:14].count('.') + ss[14:].count('.')) / len(ss)
            else:
                v = ss.count('.') / len(ss)
            return v
        except:
            time.sleep(0.25)
    raise ValueError(f"Could not compute MFEopen score for {seq}")

def get_degscore(seq, by_position=False, feature="", linearfold=False, package="eternafold", mask_U=False, N=50):
    for i in range(0,N):
        try:
            if "all" in feature:
                mdl_list = []
                for UTR in UTR_list:
                    mdl = DegScore(UTR[0] + seq + UTR[1], package=package, linear=linearfold, mask_U=mask_U)
                    mdl_list.append(mdl)
            else:
                mdl = DegScore(seq, package=package, linear=linearfold, mask_U=mask_U)
            if "all" in feature:
                output = np.mean([mdl.degscore_by_position[len(UTR[0]):-len(UTR[1])] for mdl, UTR in zip(mdl_list,UTR_list)], axis=0)
            else:
                output =  mdl.degscore_by_position

            # open 14
            if "open" in feature:
                output[:14] = 1 - output[:14]

            if by_position:
                return output
            else:
                return np.sum(output)
        except:
            time.sleep(0.25)
    raise ValueError(f"Could not compute degscore for {seq}")

def get_rg_feature(ss_constraint, feature, get_RNA_seq):
    # Fix secondary structures by removing hairpins with loop length 0
    if RNA_seq is None:
        ss = mfe(seq, package=package, constraint=ss_constraint, viterbi=True, param_file=param_file)
    else:
        ss = RNA_seq.get_mfe(ss_constraint)
    ss = ss.replace('()','..')
    rg = RGV(ss)

    if feature == 'mld':
        return rg.MLD
    elif feature == 'hp':
        return rg.n_hairpins
    elif feature == '3wj':
        return rg.n_3WJs
    elif feature == '4wj':
        return rg.n_4WJs
    elif feature == '5wj':
        return rg.n_5WJs_up
    elif feature == 'wj':
        return rg.n_3WJs + rg.n_4WJs + rg.n_5WJs_up
    elif feature == 'hp/3wj':
        # If no 3WJs then treat as if there is 0.1
        return rg.n_hairpins/max(rg.n_3WJs, 0.1)
    else:
        raise ValueError(f"Unkown rg feature {feature}.")

def bpps_try_wrapper(seq, package, param_file, linearfold, N=50):
    for i in range(0,N):
        try:
            bpm = bpps(seq, package=package, param_file=param_file, linear=linearfold)
            return bpm
        except:
            time.sleep(0.25)
    raise ValueError(f"Unable to compute bpps for {seq}")

def get_bpm_feature(seq, ss_constraint, feature, package, param_file, linearfold, RNA_seq):
    if RNA_seq is None:
        bpm = bpps_try_wrapper(seq, package, param_file, linearfold)
        
    else:
        bpm = RNA_seq.get_bpm()

    if feature == 'bpsum':
        return np.sum(bpm.flatten())/bpm.shape[0]

    elif feature == 'bpunpaired':
        bp_prob = np.sum(bpm, axis=0)
        return 1-np.mean(bp_prob)

    elif feature == 'bppaired':
        bp_prob = np.sum(bpm, axis=0)
        # For open 14 nts
        bp_prob[:14] = 1 - bp_prob[:14]
        return np.mean(bp_prob)

    elif feature == 'mcc':
        ss = mfe(seq, package=package, constraint=ss_constraint, viterbi=True, param_file=param_file)
        return get_expected_accuracy(ss, bpm, mode='mcc')
    else:
        raise ValueError(f"Unknown bpm feature {feature}.")

def get_util_feature(seq, ss_constraint, feature, RNA_seq):
    if feature == 'mbp':
        if RNA_seq is None:
            ss = mfe(seq, package=package, constraint=ss_constraint, viterbi=True, param_file=param_file)
        else:
            ss = RNA_seq.get_mfe(ss_constraint)
        return get_mean_base_pair_propensity(ss)
    if feature == 'cai':
        return calc_cai(seq)
    else:
        raise ValueError(f"Unkown util feature {feature}.")

def get_feature_value(seq, ss_constraint, feature, package, param_file, linearfold, RNA_seq, mod_U=False):
    # mld = maximal ladder distance; hp = hairpin
    # Nwj = N way junction; 5wj includes >= 5wj
    rg_feature_list = ['mld', 'hp', '3wj', '4wj', '5wj', 'wj', 'hp/3wj']
    
    # bpsum = sum of all bp probabilities
    # bpunpaired = # of unpaired based with probabilities > 50%
    # mcc = expected MCC
    bp_feature_list = ['bpsum', 'bpunpaired', 'bppaired', 'mcc']

    # mbp = mean base pair propensity
    # cai = codon adaption index (HUMAN)
    util_feature_list = ['mbp', 'cai'] 
    if "degscore" in feature:
        return get_degscore(seq, feature=feature, linearfold=linearfold, package=package, mask_U=mod_U)
    if "mfeunpaired" in feature:
        return get_MFE_open(seq, feature, package, param_file, linearfold)
    if feature in rg_feature_list:
        return get_rg_feature(ss_constraint, feature, RNA_seq)
    elif feature in bp_feature_list:
        return get_bpm_feature(seq, ss_constraint, feature, package, param_file, linearfold, RNA_seq)
    elif feature in util_feature_list:
        return get_util_feature(seq, ss_constraint, feature, RNA_seq)
    else:
        raise ValueError(f"Unknown feature {feature}. Please look at function get_feature_value for details.")
