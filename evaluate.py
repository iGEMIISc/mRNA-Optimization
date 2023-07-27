from itertools import chain, combinations
import numpy as np
import re
from arnie.pfunc import pfunc
from arnie.free_energy import free_energy
import utils_features as uf
import pickle

# Set limit for printing floats
np.set_printoptions(precision=4)

# GLOBAL CONSTANTS:
DELTA_PROB = 1e-100 # offset to preservevent log(0) when probability is 0

def mrna_score_func(seq, RNA_seq, condition, package, dna=False, Temp=37, linearfold=False, mod_U=False, **kwargs):
    """ Returns a float that represents the score of the sequence """

    # Calculate dG
    if dna:
        param_file = 'DNA'
    else:
        param_file = None

    constraint_binary = []
    constraint_float = [] # list of dG/dG_goal
    prob_list = [] # list of dG
    condition_list = [condition[f"condition_{x+1}"] for x in range(len(condition))]

    total_dG = 0
    for condition in condition_list:
        for condition_key, condition_value in condition:
            condition_key = condition_key.lower()
            # condition_value = float(condition_value)

            if 'dgopen' in condition_key:
                dG = free_energy(seq, package=package, constraint=RNA_seq.get_ss_open(), T=Temp,  param_file=param_file)
            elif 'dg' in condition_key:
                if linearfold:
                    Z = pfunc(seq, package=package, linear=True)
                    rt = 1.987e-3*(T + 273.15)
                    dG = -rt*np.log(Z)
                else:
                    dG = free_energy(seq, package=package, T=Temp,  param_file=param_file)
            else:
                ss_constraint = RNA_seq.get_ss_open()
                feature = condition_key.split('.')[0]
                feature_value = uf.get_feature_value(seq, ss_constraint, feature, package, param_file, linearfold, RNA_seq, mod_U=mod_U)
                dG = feature_value

            if condition_value is None:
                condition[0][1] = dG
                condition_value = dG
            else:
                condition_value = float(condition_value)

            # Min or Max the value
            if 'min' in condition_key:
                norm_dG = 1 - dG / condition_value
            else:
                norm_dG = dG / condition_value
            total_dG += norm_dG
            prob_list.append(dG)
    constraint_binary.append(0)
    constraint_float.append(total_dG)
    return [constraint_binary], [constraint_float], [prob_list]

def ensemble_score_func(args_d, RNA_seq = None, seq = ''):
    """ returns an ensemble score of a sequence using NUPACK to test different conditions """

    if RNA_seq:
        seq = RNA_seq.get_seq()
    return mrna_score_func(seq, RNA_seq, **args_d)

#-----------------------------------------------------------------------------------

def ensemble_score_comparison(score1, score2, T, **kwargs):
    """ compares two ensemble score to calculate probability of accepting """

    binary_score2 = np.array([x for y in score2[0] for x in y])
    if np.prod(binary_score2):
        print('Activation ratio goal is met')
        return 1, True

    dG_1 = -score1[1][0][0]
    dG_2 = -score2[1][0][0]
    delta_G = dG_2 - dG_1

    p = np.exp(-delta_G/(1.987e-3 * T))
    return p, False

def ensemble_score_dG(node):
    """ calculates the delta G for an ensemble score """
    score = node.ensemble_score

    # If ligand
    if node.RNA_seq.check_ratio_mode():
        P = score[1] # idx 1 is the float scores
        P = -np.mean(P)
    elif node.RNA_seq.check_mrna_mode():
        P = score[1]
        P = -np.mean(P)
    else: # If not ligands
        P = score[1] #idx 1 is the float scores
        P = np.array([min(x) for x in P])
        if min(P) < 0.5:
            P = np.clip(P,0,0.5) # Threshold probability at 0.5 because that's all that matters

        # Used for openTB... not sure if relevant here
        #P = np.array([P[0], P[1], min(P[2:-1]), P[-1]]) # Take the min prob for 3rd condition
        P = -sum(P)
    return P*node.RNA_seq.get_scale()