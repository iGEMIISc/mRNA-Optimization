import subprocess
import os
import shutil
import tempfile
from Bio.Seq import Seq
from arnie.mfe import mfe
import numpy as np


# GLOBAL

def uses_other_chars(s):
    """ Returns True if non-base characters are in the string """
    base_set = set('ACGU')
    return not set(s) <= base_set

def CDSFold(sequence, CDSFold_path, w=None, start_idx=None, end_idx=None, avoid_list=[], verbose=False):
    call_list = [CDSFold_path]
    
    if w is not None:
        call_list.append('-w')
        call_list.append(f'{w}')

    if start_idx is not None and end_idx is not None:
        call_list.append('-r')
        call_list.append('-f')
        call_list.append(f'{start_idx}')
        call_list.append('-t')
        call_list.append(f'{end_idx}')

    if len(avoid_list) != 0:
        call_list.append('-e')
        call_list.append(','.join(avoid_list))


    # Create temp folder safely
    fd, path = tempfile.mkstemp()
    try:
        with os.fdopen(fd, 'w') as tmp:
            tmp.write('>temp_seq\n')
            tmp.write(sequence + '\n')
            
        call_list.append(path)
        proc = subprocess.Popen(call_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.stdout.read()
    finally:
        os.remove(path)
        
    # Get raw output and extract sequence
    try:
        output_list = output.decode().split('\n')[::-1] # reverse the output

        # Get the index of MFE output since that is always 2 indicies ahead of sequence
        idx = next((idx for idx, s in enumerate(output_list) if "MFE" in s), None)
        sequence_new = output_list[idx+2]
    except:
        print(f"cmd {call_list}")
        print(f"input aa sequence: {sequence}")
        print(output)
        raise ValueError("Unable to get output from CDSFold")

    # Check if the sequence from CDSFold is valid
    if uses_other_chars(sequence_new) or len(sequence)*3 != len(sequence_new):
        print(f"cmd {call_list}")
        print(f"input aa sequence: {sequence}")
        print(f"output sequence: {sequence_new}")
        print(output_list)
        raise ValueError("Bad output from CDSFold")

    if verbose:
        return sequence_new, output_list
    return sequence_new

def rna2aa(seq):
    return str(Seq(seq).translate())

def expand_same_selection_mod(start_idx, end_idx, ss):
    """ Expands the indices to encompass unpaired bases 
    e.g.    ss = .....(((.......))).....(((....)))...
    initial idx:      |<--------->|
    final   idx: |<------------------->|
    index is inclusive [start_idx, end_idx]
    """
    while (start_idx > 0) and (ss[start_idx-1] == '.'):
        start_idx -= 1
    while (end_idx < len(ss)-1) and (ss[end_idx+1] == '.'):
        end_idx += 1
    return start_idx, end_idx

def get_selection_index_mod(ss, min_selection_length=10):
    """ From a starting origin_idx, selects the largest stem """
    selection_length = 0
    counter = 0
    selection_length_list = []
    while selection_length <= min_selection_length:
        # Start at a location of a basepair
        origin_idx = np.random.randint(len(ss))
        while ss[origin_idx] == '.':
            origin_idx = np.random.randint(len(ss))
        final_idx = origin_idx
        ss_origin_idx = ss[origin_idx]

        # Check which direction the basepair partner is
        if ss_origin_idx == '(':
            delta = 1
        else:
            delta = -1

        # Now found the index that basepairs with origin_idx
        open_count = 1
        while open_count != 0:
            final_idx += delta
            if ss[final_idx] == ss_origin_idx:
                open_count += 1
            elif ss[final_idx] == '.':
                open_count += 0
            else:
                open_count -= 1
        start_idx, end_idx = sorted([origin_idx, final_idx])

        # Expand selection to include unpaired 5' and 3' indices
        start_idx, end_idx = expand_same_selection_mod(start_idx, end_idx, ss)
        
        counter += 1
        selection_length = end_idx - start_idx + 1
        selection_length_list.append(selection_length)
        if counter > 1000:
            print(ss)
            raise ValueError(f"Unable to find suitable region for selection.")
    return start_idx, end_idx

def mutate_fold(seq, CDSFold_path, unstructured_5_prime_length=15, minimum_length=10, maximum_length=450):
    """ Run CDSFold to mutate the sequence
    *Only works for mRNA """
    new_seq = seq
    # Get structure
    ss = mfe(seq, package='eternafold', linear=True)
    counter = 0
    while seq == new_seq:

        # Select random region
        # selection_length = np.random.randint(minimum_length,min(maximum_length,len(ss))+1)
        # start_idx = np.random.randint(len(seq)-selection_length+1)
        # end_idx = start_idx+selection_length

        # select random stem
        start_idx, end_idx = get_selection_index_mod(ss, min_selection_length=minimum_length)

        # Expand selection to fit codon
        start_idx -= start_idx%3
        end_idx += 3-(end_idx+1)%3
        selection_length = end_idx-start_idx+1

        # check if selection is too large        
        if selection_length > maximum_length:
            delta = maximum_length - selection_length
            d = (selection_length - maximum_length)//3 # extra codons to remove
            d = d//2*3 # extra bases to remove on both sides
            start_idx += d
            end_idx -= d
            selection_length = end_idx-start_idx+1

        # Run CDSFold
        aa = rna2aa(seq[start_idx:end_idx+1])
        
        # Random w
        w = np.random.randint(minimum_length, selection_length+1)

        if start_idx < unstructured_5_prime_length:
            unstructured_length = int(np.ceil((15-start_idx)/3))
            rna = CDSFold(aa, CDSFold_path, w=w, start_idx=1, end_idx=unstructured_length)
        else:
            rna = CDSFold(aa, CDSFold_path, w=w)
        
        # update new_seq
        new_seq = list(new_seq)
        new_seq[start_idx:end_idx+1] = list(rna)
        new_seq = ''.join(new_seq)

        counter += 1
        if counter > 100:
            raise ValueError(f"Failed to mutate using CDSFold {seg_id}")

    return new_seq