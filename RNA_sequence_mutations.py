import numpy as np
import random
import pickle
from arnie.bpps import bpps
from arnie.mfe import mfe
from init_starting_RNA_Sequence import generate_simple_seq
from utils_sequence import match_seq_len
from utils_features import get_degscore
from copy import copy
from wrapper_CDSFold import mutate_fold

RNA_BASES = set('ACGU')
DNA_BASES = set('ACGT')

#AA2C = pickle.load(open( "reference/aa2c.pickle", "rb" ))
#C2AA = pickle.load(open( "reference/c2aa.pickle", "rb" ))
AA2C = {'F': [('UUU', 0.46), ('UUC', 0.54)], 'L': [('UUA', 0.08),  ('UUG', 0.13),  ('CUU', 0.13),  ('CUC', 0.2),  ('CUA', 0.07),  ('CUG', 0.4)], 
    'S': [('UCU', 0.19),  ('UCC', 0.22),  ('UCA', 0.15),  ('UCG', 0.05),  ('AGU', 0.15),  ('AGC', 0.24)], 'Y': [('UAU', 0.44), ('UAC', 0.56)], 
    '*': [('UAA', 0.3), ('UAG', 0.24), ('UGA', 0.47)], 'C': [('UGU', 0.46), ('UGC', 0.54)], 'W': [('UGG', 1.0)], 
    'P': [('CCU', 0.29), ('CCC', 0.32), ('CCA', 0.28), ('CCG', 0.11)], 'H': [('CAU', 0.42), ('CAC', 0.58)], 'Q': [('CAA', 0.27), ('CAG', 0.73)], 
    'R': [('CGU', 0.08),  ('CGC', 0.18),  ('CGA', 0.11),  ('CGG', 0.2),  ('AGA', 0.21),  ('AGG', 0.21)], 'I': [('AUU', 0.36), ('AUC', 0.47), ('AUA', 0.17)], 
    'M': [('AUG', 1.0)], 'T': [('ACU', 0.25), ('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.11)], 'N': [('AAU', 0.47), ('AAC', 0.53)], 
    'K': [('AAA', 0.43), ('AAG', 0.57)], 'V': [('GUU', 0.18), ('GUC', 0.24), ('GUA', 0.12), ('GUG', 0.46)], 
    'A': [('GCU', 0.27), ('GCC', 0.4), ('GCA', 0.23), ('GCG', 0.11)], 'D': [('GAU', 0.46), ('GAC', 0.54)], 'E': [('GAA', 0.42), ('GAG', 0.58)], 
    'G': [('GGU', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)]}

C2AA = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 
    'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGU': 'C', 'UGC': 'C', 
    'UGA': '*', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 
    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 
    'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 
    'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 
    'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 
    'GGG': 'G'}

def codon2value(codon):
    """ please optimize. not worth calling every single time """
    value = 0
    for base in codon:
        if base in ['C','G']:
            value += 2
        else:
            value += 1
    return value

def mutate_codon(RNA_seq, seg_id, category, options, mod_U=False):
    """ Only works for RNA """
    seq = RNA_seq.RNA_frag[seg_id]
    counter = 0
    while seq == RNA_seq.RNA_frag[seg_id]:

        seq = list(RNA_seq.RNA_frag[seg_id])
        start_idx, end_idx = RNA_seq.get_segment_bounds(seg_id)

        # Choose codon to mutate
        if options=='p': # choose based on probability unpaired 
            bpm = RNA_seq.get_bpm()

            
            p_idx = 1- np.sum(bpm, axis=0) # P(idx) = P(unpaired) = 1 - P(paired)
            p_idx = np.clip(p_idx,0,1) # sometimes p_idx can be a small negative number
            
            p_idx = p_idx[start_idx:end_idx] # Take slice of only the domain

            idx = np.random.choice(np.arange(len(seq)), p=p_idx/np.sum(p_idx))
            idx = idx - idx%3 # to get the index to align with codons

        elif options=='d': # choose based on degscore  

            structure = RNA_seq.get_mfe()      

            mdl_by_position = get_degscore(''.join(seq), by_position=True, mod_U=mod_U)
            p_idx = np.exp(mdl_by_position)

            p_idx = p_idx[start_idx:end_idx] # Take slice of only the domain

            idx = np.random.choice(np.arange(len(seq)), p=p_idx/np.sum(p_idx))
            idx = idx - idx%3 # to get the index to align with codons


        elif options=='s': # choose unpaired regions in secondary structure
            ss = RNA_seq.get_mfe()

            N_d = ss.count('.')
            N_p = len(ss) - N_d
            p_idx = np.array([0.9/N_d if c=='.' else 0.1/N_p for c in ss])
            
            p_idx = p_idx[start_idx:end_idx] # Take slice of only the domain

            idx = np.random.choice(np.arange(len(seq)), p=p_idx/np.sum(p_idx))
            idx = idx - idx%3 # to get the index to align with codons

        elif isinstance(options,int):
            # There's a restriction site or repeat and we've been specifically tasked
            #to mutate a codon here to make it go away.
            idx = options - options % 3

            # jiggle a little around
            perturbation = np.random.choice([-3,0,3]) 

            idx += perturbation

            idx = max(idx, 3) # not allowed to mutate start codon
            idx = min(idx, (end_idx - end_idx % 3)) # don't go past the final one

        else: # choose randomly
            idx = np.random.randint(len(seq)/3)*3
        
        codon = ''.join(seq[idx:idx+3])

        # Pick out a new codon for the amino acid
        aa = C2AA[codon]
        codon_list = [x[0] for x in AA2C[aa]]

        if 'f' in category: # Choose codon based on natural frequency
            p_list = np.array([x[1] for x in AA2C[aa]])
            codon = np.random.choice(codon_list, 1, p=p_list/sum(p_list))[0]

        elif 'g' in category: # Choose codon based on GC content
            p_list = np.array([codon2value(x) for x in codon_list])
            codon = np.random.choice(codon_list, 1, p=p_list/sum(p_list))[0]

        else: # Choose randomly
            codon = np.random.choice(codon_list, 1)[0]

        seq[idx:idx+3] = list(codon)
        seq = ''.join(seq)

        counter += 1
        if counter > 10000:
            raise ValueError(f"Failed to mutate codon at location %d" % idx)

    RNA_seq.RNA_frag[seg_id] = seq
    return RNA_seq, 'Mutate:{}'.format(seg_id)

def shuffle(num_shuffle, RNA_seq, idx_swap = [], swap = False):
    """ shuffles the segment order array 

    Args:
        num_shuffle: number of segments to shuffle
        RNA_seq: An instance of RNA_Sequence class (contains seg_order)
        idx_swap: Specified segments to swap
        swap: Whether to randomly insert the segments or to only switch them

    Return:
        seg_order_new: A new segment ordering
    """
    seg_order = RNA_seq.seg_order
    seg_order_new = np.array(seg_order)
    if num_shuffle > len(seg_order_new):
        num_shuffle = len(seg_order_new)
    if not idx_swap:
        idx_swap = np.random.choice(seg_order_new, num_shuffle, False)
    else:
        num_shuffle = max(len(idx_swap), num_shuffle)
        temp_arr = list(set(seg_order_new) - set(idx_swap))
        idx_swap += list(np.random.choice(temp_arr, num_shuffle - len(idx_swap), False))
        num_shuffle = len(idx_swap)
    if swap and num_shuffle > 1:
        while np.array_equal(seg_order_new, seg_order):
            np.random.shuffle(idx_swap)
            idx_temp = list(idx_swap)
            seg_order_new = np.array([x if x not in idx_swap else idx_temp.pop(0) for x in seg_order_new])
    else:
        while np.array_equal(seg_order_new, seg_order):
            for x in idx_swap:
                seg_order_new = np.delete(seg_order_new, np.where(seg_order_new == x)[0])
            for x in idx_swap:
                seg_order_new = np.insert(seg_order_new, np.random.randint(0,len(seg_order_new)+1), x)
    return seg_order_new, 'Reorder:{}'.format(','.join(idx_swap))

def mutate_fragment(RNA_seq, seg_id, constant_length, preserve_seq, num_mutate=1, \
        length_limit = 0.8, DNA=False, category=None, options=None, \
        CDSFold_path=None, CDSFold_prob=0, mod_U=False):
    """ mutates the segment in RNA_seq specified by seg_id
    
    This will either (0) mutate or (-1) remove a base from the segment (seg_id, e.g. 'A').
    If a base it removed it will be added to the adjacent segment. Segments cannot be smaller
    than 4 bases.

    Currently not set up to run on DNA sequences (CDSFold is only mRNA)
    """
    RNA_frag = RNA_seq.RNA_frag
    seg_order = RNA_seq.seg_order


    if category in ['aa', 'mrna']: # if mRNA only mutate codon or CDSFold
        if np.random.rand() < (1 - CDSFold_prob):
            for i in range(np.random.randint(1,num_mutate+1)):
                RNA_seq, move = mutate_codon(RNA_seq, seg_id, category, options, mod_U=mod_U)
            return RNA_seq, f"({i+1}){move}"
        else:
            RNA_seq.RNA_frag[seg_id] = mutate_fold(RNA_seq.RNA_frag[seg_id], CDSFold_path)
            move = "CDSFold"
            return RNA_seq, f"{move}"

    option = [0]
    if not constant_length:
        cur_len = len(RNA_frag[seg_id])
        og_len = RNA_seq.seg_len_og[seg_id]

        seg_id_idx = np.where(seg_order == seg_id)[0][0]
        # Determine if base can be removed
        if cur_len > og_len*length_limit and cur_len >= 5 and len(set(seg_order) - set(preserve_seq)) > 1: # I don't want segments to be < 4 bases @@@
            option.append(-1)
    option = np.random.choice(option)

    for i in range(np.random.randint(1,num_mutate+1)):
        if option == 0: # mutate a base
            RNA_seq, move = mutate_seq(RNA_seq, seg_id, DNA=DNA)
        elif option == -1: # remove a base
            RNA_seq, move = mutate_edge_len(RNA_seq, seg_id, preserve_seq, DNA=DNA)
    return RNA_seq, f"({i+1}){move}"

# def mutate_seq(RNA_seq, seg_id, DNA=False):
#     """ mutates the sequence of a RNA segment """
#     if DNA:
#         base_set = DNA_BASES
#     else:
#         base_set = RNA_BASES

#     seq = list(RNA_seq.RNA_frag[seg_id])
#     idx = np.random.randint(0, len(seq))
#     b_new = random.sample(base_set - set(seq[idx]) , 1)[0]
#     seq[idx] = b_new
#     RNA_seq.RNA_frag[seg_id] = ''.join(seq)
#     return RNA_seq, 'Mutate:{}'.format(seg_id)

# def mutate_edge_len(RNA_seq, seg_id, preserve_seq, DNA=False):
#     """ mutates the length of a RNA segment 

#     Removes a base from the segment (seg_id) and a base is added to a random
#     segment. Either the first or last base is moved.
#     """
#     # Remove a base from the segment
#     if DNA:
#         base_set = DNA_BASES
#     else:
#         base_set = RNA_BASES

#     RNA_frag = RNA_seq.RNA_frag
#     seq = list(RNA_seq.RNA_frag[seg_id])
#     seq.pop(np.random.choice([0,-1]))
#     RNA_frag[seg_id] = ''.join(seq)

#     # Add a random base to a random fragment, can't be the original segment
#     seg_id_added = random.sample(set(RNA_seq.seg_order) - set([seg_id]) - set(preserve_seq), 1)[0]
#     assert(seg_id_added!=seg_id)

#     b_new = random.sample(base_set,1)[0]
#     seq = list(RNA_frag[seg_id_added])
#     seq.insert(np.random.choice([0, len(seq)]),b_new) # insert randomly @@@
#     RNA_frag[seg_id_added] = ''.join(seq)
#     return RNA_seq, 'Length:{}->{}'.format(seg_id, seg_id_added)

# def swap_fragment(RNA_seq, seg_id, args_d):
#     """ replaces a segment of the RNA with a newly generated segment """
#     # First swap in a newly generated segment
#     frag_len = {}
#     for x in RNA_seq.seg_order:
#         frag_len[x] = len(RNA_seq.RNA_frag[x])

#     RNA_frag_temp, seg_order_temp = generate_simple_seq(seg_order_start = RNA_seq.seg_order,
#         frag_len_start = frag_len, **args_d)

#     length = len(RNA_seq.RNA_frag[seg_id])
#     new_frag = match_seq_len(RNA_frag_temp[seg_id], length, DNA=args_d['dna']) # Turn the sequence into the proper length
    
#     RNA_seq.RNA_frag[seg_id] = new_frag # Turn the sequence into the proper length
#     assert(length == len(new_frag))

#     # Next mutate it but maybe unnecessarily
#     #RNA_seq.RNA_frag[seg_id] = mutate_edge(RNA_seq.RNA_frag[seg_id])
#     return RNA_seq, 'Swap:{}'.format(seg_id)