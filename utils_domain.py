import utils_sequence as us
import pickle
import numpy as np

def _generate_sequence(sequence, length, category, options, parent_sequence=None, allow_nucleotide_repeat=False):
    dna = category == 'dna'
    freq = 'f' in options
    if category == 'aa': # amino acid
        if length != len(sequence)*3:
                raise ValueError(f"Incorrect length for mRNA. Easy solution is to use 'F'.")
        return _aa2rna(sequence, allow_nucleotide_repeat, freq=freq)

    elif category == 'mrna': # if mRNA
        if length != len(sequence):
            raise ValueError(f"Incorrect length for mRNA. Easy solution is to use 'F'.")
        return _aa2rna(None, allow_nucleotide_repeat, freq=freq, starting_seq=sequence)

    elif category == 'rna' or category == 'dna':
        if 'g' in options:
            return _generate_filled_sequence(length, sequence, category, DNA=dna)
        if 'c' in options:
            return _generate_complement(length, category, options, parent_sequence, DNA=dna)
        elif 'i' in options:
            return us.match_seq_len(parent_sequence, length, DNA=dna)
        else:
            return us.match_seq_len(sequence, length, DNA=dna)
    else:
        raise ValueError(f"Unknown category {category}.")

def _generate_filled_sequence(length, sequence, category, DNA=False):
    if sequence == '*H':
        filled_sequence = us.gen_random_hairpin_sequence(length, DNA=DNA)
    elif sequence == '*N':
        filled_sequence = us.gen_random_sequence(length, DNA=DNA)
    elif 'N' in sequence:
        filled_sequence = us.gen_N_sequence(sequence, DNA=DNA)
    else:
        raise ValueError(f"Unkown sequenced to be filled: {sequence}")
    return filled_sequence

def _generate_complement(length, category, options, parent_sequence, DNA=False):

    gu_level = 0
    return us.gen_complement(parent_sequence, length, gu_level, DNA=DNA)

def _process_length(sequence, length, category, options):
    if length == 'F':
        if 'c' in options:
            raise ValueError("Must specify actual length is domain is a complement of another domain. 'F' is not valid.")
        if category == 'aa':
            length = len(sequence)*3
        else:
            length = len(sequence)
    else:
        length = length
    return length

def _aa2rna(aa_seq, allow_nucleotide_repeat, freq=False, starting_seq=""):
    """ Converts amino acid sequence to RNA """
    # Load dictionary of amino acid to codons
    #AA2C = pickle.load(open( "reference/aa2c.pickle", "rb" ))
    #C2AA = pickle.load(open( "reference/c2aa.pickle", "rb" ))
    aa2c = {'F': [('UUU', 0.46), ('UUC', 0.54)], 'L': [('UUA', 0.08),  ('UUG', 0.13),  ('CUU', 0.13),  ('CUC', 0.2),  ('CUA', 0.07),  ('CUG', 0.4)], 
    'S': [('UCU', 0.19),  ('UCC', 0.22),  ('UCA', 0.15),  ('UCG', 0.05),  ('AGU', 0.15),  ('AGC', 0.24)], 'Y': [('UAU', 0.44), ('UAC', 0.56)], 
    '*': [('UAA', 0.3), ('UAG', 0.24), ('UGA', 0.47)], 'C': [('UGU', 0.46), ('UGC', 0.54)], 'W': [('UGG', 1.0)], 
    'P': [('CCU', 0.29), ('CCC', 0.32), ('CCA', 0.28), ('CCG', 0.11)], 'H': [('CAU', 0.42), ('CAC', 0.58)], 'Q': [('CAA', 0.27), ('CAG', 0.73)], 
    'R': [('CGU', 0.08),  ('CGC', 0.18),  ('CGA', 0.11),  ('CGG', 0.2),  ('AGA', 0.21),  ('AGG', 0.21)], 'I': [('AUU', 0.36), ('AUC', 0.47), ('AUA', 0.17)], 
    'M': [('AUG', 1.0)], 'T': [('ACU', 0.25), ('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.11)], 'N': [('AAU', 0.47), ('AAC', 0.53)], 
    'K': [('AAA', 0.43), ('AAG', 0.57)], 'V': [('GUU', 0.18), ('GUC', 0.24), ('GUA', 0.12), ('GUG', 0.46)], 
    'A': [('GCU', 0.27), ('GCC', 0.4), ('GCA', 0.23), ('GCG', 0.11)], 'D': [('GAU', 0.46), ('GAC', 0.54)], 'E': [('GAA', 0.42), ('GAG', 0.58)], 
    'G': [('GGU', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)]}

    c2aa = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 
    'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGU': 'C', 'UGC': 'C', 
    'UGA': '*', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 
    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 
    'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 
    'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 
    'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 
    'GGG': 'G'}

    # Creating starting RNA sequence from amino acids
    if len(starting_seq) == 0:
        rna = []
        for aa in aa_seq:
            codon_list = [x[0] for x in aa2c[aa]]
            if freq: # Choose codon based on natural frequency
                p_list = np.array([x[1] for x in aa2c[aa]])
                codon = np.random.choice(codon_list, 1, p=p_list/sum(p_list))[0]
            else: # Choose randomly
                codon = np.random.choice(codon_list, 1)[0]
            rna.append(codon)
    else:
        assert len(starting_seq) % 3 == 0, "Starting mRNA sequence must have length divisible by 3."
        rna = [starting_seq[i:i+3] for i in range(0, len(starting_seq), 3)]
        aa_seq = ''.join([c2aa[c] for c in rna])
    
    # Fix nucleotide repeats (no stretches of A,C,G,U)
    repeat_list = us.get_nucleotide_repeat_list()
    counter = 0
    rna_seq = ''.join(rna)
    while us.check_nucleotide_repeat(rna_seq) and not allow_nucleotide_repeat:
        counter += 1        
        # find the codons that need to be swapped
        for nr in repeat_list:
            idx = rna_seq.find(nr)
            if idx != -1: # found a repeat
                idx = idx//3 + np.random.randint(0,2)
                
                # swap codon
                aa = aa_seq[idx]
                codon_list = [x[0] for x in aa2c[aa]]
                cur_codon_idx = codon_list.index(rna[idx])
                if freq: # Choose codon based on natural frequency
                    p_list = np.array([x[1] for x in aa2c[aa]])
                    if len(p_list) > 1:
                        p_list[cur_codon_idx] = 0 # Remove starting codon
                    codon = np.random.choice(codon_list, 1, p=p_list/sum(p_list))[0]
                else: # Choose randomly
                    if len(codon_list) > 1:
                        del codon_list[cur_codon_idx]
                    codon = np.random.choice(codon_list, 1)[0]
                rna[idx] = codon
        rna_seq = ''.join(rna)
        
        if counter >= 10000:
            raise ValueError(f"Unable to generate an mRNA sequence for {rna_seq}.")
    return rna_seq