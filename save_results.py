import csv
import numpy as np
import pandas as pd
import time
from plot_tree import plot_pickle_save
import os
from utils_features import calc_cai
from DegScore import DegScore

def append_full_sequence(df, constant_5_prime, constant_3_prime):
    idx = df.columns.get_loc("Accept") + 1 # insert after "Accept"
    num_col = sum(['frag_' in x for x in df.columns])
    sequence = constant_5_prime + np.sum(df.iloc[:,idx:idx+num_col].values,axis=1) + constant_3_prime
    df.insert(idx, 'Sequence', sequence)
    return df

def save_results(args_d, root, sol_node, MC_results, save_full=True):

    if save_full:
        full_outfile = _generate_outfile_name(args_d, 'FULL', timestamp=False)
        save_full_MC_results_to_csv(root, sol_node, args_d, full_outfile, MC_results)

    best_outfile = _generate_outfile_name(args_d, 'RUNNING_BEST', timestamp=False)
    save_best_MC_results_to_csv(sol_node, args_d, best_outfile)

def save_best_MC_results_to_csv(sol_node, args_d, outfile):
    """ saves the Monte Carlo results as a tab-delimited txt """
    MC_results = save_node(sol_node, MC_results = [])

    #print('saving best result\n')
    column_headers = _compile_column_headers(args_d, sol_node, IDs = True)

    # Appending 5' and 3' prime sequences
    c5prime = sol_node.RNA_seq.c5prime
    c3prime = sol_node.RNA_seq.c3prime
    #sequence = MC_results[0][2]
    sequence = sol_node.RNA_seq.get_seq()
    #if sequence != MC_results[0][2]:
    #    print('sequence was,',sequence,'MC_results[0][2] was,',MC_results[0][2])

    bp_matrix = sol_node.RNA_seq.get_bpm()
    mfe_struct = sol_node.RNA_seq.get_mfe()
    RiboGraphViz = sol_node.RNA_seq.get_RGV()

    AUP = np.mean(1 - np.sum(bp_matrix, axis=0))
    AUP_init14 = np.mean(1 - np.sum(bp_matrix, axis=0)[:14])

    #full_sequence = c5prime+sequence+c3prime
    degscore_mdl = DegScore(sequence, structure=mfe_struct, mask_U=args_d['mod_U'])
    degscore = degscore_mdl.degscore

    rg_mdl = sol_node.RNA_seq.get_RGV()
    MLD = rg_mdl.MLD

    dG_MFE = sol_node.RNA_seq.get_dG_MFE()

    if len(c3prime)==0:
        CDS_sequence = sequence[len(c5prime):]
    else:
        CDS_sequence = sequence[len(c5prime):-1*len(c3prime)]

    output_dict = {'CDS_sequence': [CDS_sequence],'full_sequence':[sequence], 'CAI': [calc_cai(CDS_sequence)],\
    'AUP': [AUP], 'AUP_init14': [AUP_init14], 'MFE Structure': [mfe_struct], 'dG(MFE)': [dG_MFE],\
    'DegScore': [degscore],'MLD': [MLD], 'c5prime': [c5prime], 'c3prime': [c3prime]}

    df = pd.DataFrame(output_dict) 
    # except:
    #     print("Warning: Something went wrong with the column labels")
    #     df = pd.DataFrame(MC_results)

    #append_full_sequence(df, args_d['constant_5_prime'], args_d['constant_3_prime'])

    if os.path.exists(args_d['output'] + outfile):
        df_old = pd.read_csv(args_d['output'] + outfile, delimiter='\t')
        df = pd.concat([df_old, df])

    df.to_csv(args_d['output'] + outfile, sep='\t', quoting=csv.QUOTE_NONE, index=False)


def save_full_MC_results_to_csv(root, sol_node, args_d, outfile, MC_results):
    """ saves the Monte Carlo results as a tab-delimited txt """

    print('saving full result\n')
    column_headers = _compile_column_headers(args_d, sol_node, IDs = True)

    try:
        df = pd.DataFrame(MC_results, columns=column_headers) 
        df = df[['ID','Accept','frag_1','raw_value_1']]
        df = df.rename(columns={'frag_1':'sequence', 'raw_value_1':'loss'})
    except:
        print("Warning: Something went wrong with the column labels")
        df = pd.DataFrame(MC_results)

    if os.path.exists(args_d['output'] + outfile):
        df_old = pd.read_csv(args_d['output'] + outfile, delimiter='\t')
        df = df_old.append(df, ignore_index=True)
        
    df.to_csv(args_d['output'] + outfile, sep='\t', quoting=csv.QUOTE_NONE, index=False)


def _compile_column_headers(args_d, sol_node, IDs):
    RNA_seq = sol_node.get_RNA_seq()

    column = []
    if IDs is True:
        column.append("ID")
    
    column.append('Accept')

    num_frag = len(RNA_seq.seg_order)
    [column.append('frag_{}'.format(i+1)) for i in range(num_frag)]
    [column.append('segment_{}'.format(i+1)) for i in range(num_frag)]

    num_condition = len(args_d['condition'])
    score = RNA_seq.get_score()
    for i in range(len(score[0])): # Iterate through binary values
        for j in range(len(score[0][i])):
            column.append(f"binary_constraint_{i+1}_{j+1}")
    for i in range(len(score[1])): # Iterate through float values
        for j in range(len(score[1][i])):
            column.append(f"float_constraint_{i+1}_{j+1}")

    num_float = len(sol_node.get_MC_result()[0]) - len(column) - 1
    if IDs is True:
        num_float += 1
    for i in range(num_float):
        column.append('raw_value_{}'.format(i+1))

    column.append('move')
    return column

def _generate_outfile_name(args_d, append, timestamp=False):
    if timestamp:
        ts = time.strftime("%Y%m%d-%H%M%S")
        input_file = args_d['input_file'].split('/')[-1].replace('.txt','')
        return '{}.{}-{}.{}.txt'.format(input_file, ts, os.getpid(), append)
    else:
        input_file = args_d['input_file'].split('/')[-1].replace('.txt','')
        return '{}.{}.{}.txt'.format(input_file, args_d['unique_id'], append)

def save_node(node, MC_results):
    ID = node.get_id()
    for result in node.get_MC_result():
        MC_results.append([ID] + result)
    return MC_results
