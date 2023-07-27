import numpy as np
import random
from Bio import pairwise2


RNA_BASES = set('ACGU')
DNA_BASES = set('ACGT')


def _base_complement(base, gu_fraction, DNA=False):
    """ get complement of single base.
        Args:
            - base: single letter (string, uppercase)
            - gu_fraction: probability of GU/UG pairing (float)
            - DNA: does the sequence use T instead of U? (boolean)
        Returns:
            - single letter (string)
    """
    if base == 'A':
        return 'T' if DNA else 'U'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'U' if np.random.uniform(0,1) < gu_fraction else 'C'
    elif base == 'U':
        return 'G' if np.random.uniform(0,1) < gu_fraction else 'A'
    elif base == 'T':
        return 'A'
    else:
        raise ValueError('Unknown base is used. Only A,C,G,U, and T bases are accepted.')


def reverse_complement(seq, gu_fraction, DNA):
    """ get reverse complement of a sequence string.
        Args:
            - seq: the sequence to reverse-complement (string)
            - gu_fraction: fraction of G/U bases where the rev. comp. will have U/G (float)
            - DNA: does the sequence use T instead of U? (boolean)
        Returns:
            - the input sequence, reverse complemented
    """
    seq_reverse = seq[::-1]
    seq_rev_comp = [_base_complement(base, gu_fraction, DNA=DNA) for base in seq_reverse]
    return ''.join(seq_rev_comp)


def match_seq_len(seq, length, DNA=False):
    """ given a sequence it will add/remove bases until its length is met """
    if DNA:
        base_set = DNA_BASES
    else:
        base_set = RNA_BASES

    cur_length = len(seq)
    if cur_length == length:
        return seq
    if cur_length < length: # gotta add bases
        seq = list(seq)
        for i in range(length - cur_length):
            base = random.sample(base_set,1)[0]
            idx = np.random.choice([0,len(seq)])
            seq.insert(idx,base)
        seq = ''.join(seq)
    else: # gotta remove bases
        idx = np.random.randint(cur_length - length)
        seq = seq[idx:idx+length]
    return seq


def gen_complement(seq, target_len, gu_fraction, DNA=False):
    """ generate complement to sequence of length target_len """
    """ get reverse complement of a sequence string. If the sequence
        is shorter than target_len, increase by appending random bases
        on both ends. If the sequence is longer than target_len,
        randomly remove bases, with some chance of removing a base in
        the middle of the sequence if can_delete_middle is True.
        Args:
            - seq: the sequence to reverse-complement (string)
            - gu_fraction: fraction of G/U bases where the rev. comp. will have U/G (float)
            - DNA: does the sequence use T instead of U? (boolean)
        Returns:
            - the input sequence, reverse complemented (string)
    """
    
    seq_matched_len = match_seq_len(seq, target_len, DNA)
    revcomp = reverse_complement(seq_matched_len, gu_fraction, DNA)

    assert len(revcomp) == target_len, (len(revcomp), target_len)
    return revcomp

def gen_random_sequence(length, DNA = False):
    if DNA:
        return ''.join(np.random.choice(list(DNA_BASES), length))
    return ''.join(np.random.choice(list(RNA_BASES), length))

def gen_random_hairpin_sequence(total_length, DNA = False):
    hairpin_first_half = gen_random_sequence(total_length//2 + total_length%2)
    hairpin_second_half = gen_complement(hairpin_first_half, total_length//2, 0, DNA=DNA)
    return hairpin_first_half + hairpin_second_half

def gen_N_sequence(sequence, DNA = False):
    sequence = list(sequence)

    for idx, c in enumerate(sequence):
        if c == 'N':
            if DNA:
                sequence[idx] = np.random.choice(list(DNA_BASES), 1)[0]
            else:
                sequence[idx] = np.random.choice(list(RNA_BASES), 1)[0]
    return ''.join(sequence)

def combine_sequence(seq1, seq2):
    """ Computes an approximate consensus sequence
        for two sequences of DNA/RNA by aligning them together.
        Args:
            - seq1, seq2: the two sequences to combine (string)
        Returns:
            - the consensus sequence
    """

    # Set up match vs. mismatch scoring matrix,
    # and specifically allow for GU pairing by
    # considering A and G (or C and U) equivalent.
    matrix = {}
    for base1 in 'ACGUT':
        for base2 in 'ACGUT':
            bases = base1 + base2
            if base1 == base2:
                score = 1
            elif bases == 'AG' or bases == 'GA':  # A and G pair with U
                score = 1
            elif bases == 'CU' or bases == 'UC':  # C and U pair with G
                score = 1
            else:
                score = 0
            matrix[(base1, base2)] = score
    
    # Align
    alignments = pairwise2.align.globalds(seq1,seq2, matrix, -0.5, -0.5)

    # Retrieve consensus sequence from the top alignment
    combined_seq = ''
    for base1, base2 in zip(alignments[0][0], alignments[0][1]):
        bases = base1 + base2
        if base1 == base2:
            base_to_add = base1
        elif bases == 'AG' or bases == 'GA':
            base_to_add = 'A'
        elif bases == 'CU' or bases == 'UC':
            base_to_add = 'C'
        elif base1 == '-':
            base_to_add = base2
        elif base2 == '-':
            base_to_add = base1
        else:
            base_to_add = np.random.choice([base1,base2], 1)[0]
        combined_seq += base_to_add
    return combined_seq

def get_nucleotide_repeat_list():
    return ['AAAAA', 'CCCC', 'GGGG', 'UUUU', 'TTTT']

def check_restriction_sites(seq, restriction_sites):
    """ Checks if the sequence contains any of the prohibited restriction sites.
    If restriction site in seq: returns seq location of that site
    else: returns -1 """
    for x in restriction_sites:
        if x in seq:
            return seq.find(x)

    return -1

def check_nucleotide_repeat(seq):
    """ Checks if there are any stretches of nucleotide repeats
        of A, C, G and U """

    nucleotide_repeat = get_nucleotide_repeat_list()

    for repeat in nucleotide_repeat:
        if repeat in seq:
            return seq.find(repeat)

    return -1