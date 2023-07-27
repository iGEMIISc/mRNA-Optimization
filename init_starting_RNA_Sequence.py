import numpy as np
import random
import re
from RNA_Sequence import *
from utils_sequence import *


RNA_BASES = set('ACGU')
DNA_BASES = set('ACGT')

def _get_normally_distributed_length(len_input_str):
    len_inputs = len_input_str.split('*')
    avg, std_dev = len_inputs[0], len_inputs[1]
    # Must have avg, std_dev is optional
    avg = float(avg)
    try:
        std_dev = float(std_dev)
    except:
        print('Default standard devation for segment_length. 25% of mean')
        std_dev = avg / 4
    return _sample_from_normal_dist(avg, std_dev)
        
    
def _sample_from_normal_dist(avg, std_dev):
    num = np.round(np.random.normal(avg, std_dev)).astype(int)
    while num <= 0:
        num = np.round(np.random.normal(avg, std_dev)).astype(int)
    return num


def _parse_input_segment_lengths(segment_IDs, input_segment_lengths):
    # Fill in segment length dictionary using lengths given in segment_lengths_list,
    # when possible. Add segment IDs to segments_without_lengths if no length info
    # ("*", "*A", etc.) was input for that segment.
    given_segment_lengths = {}
    segments_with_norm_dist_len = {}
    segments_with_random_len = set()
    for segment_ID, length_input in zip(segment_IDs, input_segment_lengths):
        try:
            segment_len = int(length_input)
            if segment_len <= 0:
                raise ValueError("Input segment length <= 0: " + str(length_input))
            given_segment_lengths[segment_ID] = segment_len
            
        except:
            if '*' not in length_input:
                raise ValueError("Input fragment length cannot be parsed: " + str(length_input))             
            if length_input.startswith('*'):
                # fragment length will be determined by the number of bases leftover
                segments_with_random_len.add(segment_ID)
            else:
                # length_input must be in "avg*std_dev" format
                # fragment length will be determined by sampling from normal dist.
                segments_with_norm_dist_len[segment_ID] = length_input
                
    return given_segment_lengths, segments_with_norm_dist_len, segments_with_random_len


def _get_random_segment_lengths(segments_with_norm_dist_len, segments_with_random_len, target_length):
    # Try to assign lengths to the rest of the segments, using random sampling
    # Stop when goal RNA_length is reached, or after 1000 attempts
    count = 0
    while True:
        random_segment_lengths = {}
        for segment_ID in segments_with_norm_dist_len.keys():
            length_input = segments_with_norm_dist_len[segment_ID]
            random_segment_lengths[segment_ID] = _get_normally_distributed_length(length_input)
                    
        length_so_far = sum(random_segment_lengths.values())
        
        if len(segments_with_random_len) > 0:
            avg_base_left = target_length / len(segments_with_random_len)
            std_dev = avg_base_left*0.1
            for segment_ID in segments_with_random_len:
                random_segment_lengths[segment_ID] = _sample_from_normal_dist(avg_base_left, std_dev)

        total_length = sum(random_segment_lengths.values())

        if total_length == target_length or target_length == 0:
            break
        
        count += 1
        if count == 1000:
            raise ValueError("Cannot generate sequence of length specified due to incompatible length constraints")
   
    return random_segment_lengths

def determine_segment_lengths(segment_IDs, target_length, input_segment_lengths = []):
    """ Parses input segment lengths, when given, and randomly samples to get lengths
        for the remaining segments.

        First, lengths for all segments where an integer was given as the input length
        are assigned.
        Second, lengths for all segments where a mean (and optional standard deviation)
        length were specified are determined by sampling from a normal distribution.
        Third, lengths for all segments left are determined by sampling from a normal
        distribution, with mean equal to the number of bases needed to reach the target
        length, divided by the number of segments remaining.

        If an ampty list is input for input_segment_lengths, all segment lengths are
        randomly sampled (third step above).

        Args:
            - segment_IDs: list of IDs for each segment, in the same order as 
                input_segment_lengths
            - target_length: the desired length of the entire RNA sequence (sum of the
                lengths of all the segments)
            - input_segment_lengths: optional, list of input lengths for each segment;
                can be a number or a string in various formats (TODO)
        Returns:
            - segment_lengths: dictionary mapping segment IDs to their lengths
    """
    if len(input_segment_lengths) == 0:
        # TODO: figure out if this ever really happens
        segments_with_norm_dist_len = {}  # empty --> will be skipped in func
        segments_with_random_len = set(segment_IDs)
        return _get_random_segment_lengths(segments_with_norm_dist_len,
                                           segments_with_random_len,
                                           target_length)
            
    (segment_lengths,
     segments_with_norm_dist_len,
     segments_with_random_len) = _parse_input_segment_lengths(segment_IDs,
                                                              input_segment_lengths)
    
    length_before_random_sampling = sum(segment_lengths.values())

    if target_length > 0: # if 0 then any length is acceptable
        target_length_left = target_length - length_before_random_sampling
    else:
        target_length_left = 0

    random_segment_lengths = _get_random_segment_lengths(segments_with_norm_dist_len,
                                                         segments_with_random_len,
                                                         target_length_left)
     
    segment_lengths.update(random_segment_lengths)
    return segment_lengths


def generate_segment_sequences(domain_list, domain_dict, gu_level, allow_nucleotide_repeat):
    # TODO: implement gu_level
    # 0. Reset all domains in case they already have sequences generated
    for key, domain in domain_dict.items():
        domain.reset()

    # 1. Generate all domains that do not rely on other domains
    domain_set = set([domain_dict[key] for key in domain_list])
    for k, v in domain_dict.items():
        if not v.check_dependent():
            v.generate_init_sequence(allow_nucleotide_repeat=allow_nucleotide_repeat)
            if v in domain_set:
                domain_set.remove(v)

    # 2. Generate the rest of the domains
    counter = 0
    while len(domain_set) > 0:
        counter += 1

        # Iterate through domain list
        for key, domain in domain_dict.items():
            if domain.check_processed(): continue

            parent_domain_name_list = domain.get_og_sequence().split('+')
            parent_domain_list = [domain_dict[name] for name in parent_domain_name_list]

            if all([parent_domain.check_processed() for parent_domain in parent_domain_list]):
                parent_sequence_list = [domain.get_sequence() for domain in parent_domain_list]

                if len(parent_sequence_list) == 1:
                    parent_sequence = parent_sequence_list[0]
                elif len(parent_sequence_list) == 2:
                    parent_sequence = combine_sequence(parent_sequence_list[0], parent_sequence_list[1])
                else:
                    raise ValueError("Cannot combine more than 2 sequences. A+B is valid but A+B+C is not.")

                domain.generate_init_sequence(parent_sequence=parent_sequence, allow_nucleotide_repeat=allow_nucleotide_repeat)
                domain_set.remove(domain)
        if counter >= 1000:
            raise ValueError("Unable to generate sequences for all domains.")
    return

def generate_simple_seq(domain_list, domain, length, seg_order_start=[],
        frag_len_start={}, gu_level=0, domain_order=False, dna=False, allow_nucleotide_repeat=False, **kwargs):

    # Get the segment order (and replace all *)
    if len(seg_order_start) == 0:        
        if domain_order:
            d_order = domain_list
        else:
            d_order = np.copy(domain_list)
            np.random.shuffle(d_order)
    else:
        d_order = seg_order_start

    # Get each domain length
    if len(frag_len_start) == 0:
        domain_length_list = [domain[key].get_og_length() for key in domain_list]
        frag_len = determine_segment_lengths(domain_list, length, domain_length_list)
    else:
        frag_len = frag_len_start

    # Save the domain length
    for k,v in frag_len.items():
        domain[k].set_length(v)
    
    # Generate actual sequence for each domain
    generate_segment_sequences(domain_list, domain, gu_level, allow_nucleotide_repeat)

    RNA_frag = {}
    for d in d_order:
        RNA_frag[d] = domain[d].get_sequence()
    
    return RNA_frag, d_order


def gen_RNA_start_simple(args_d):
    """ generates a starting sequence for starting a monte carlo search

    Args:
        seed: RNG seed
        puzzle: DMSSD puzzle object for sequences to bind
        seg_order: Determines the order of the RNA segments
        frag_len: Determines the length of each RNA segment

    Returns:
    RNA_frag: Dictionary where key is segment id (e.g. 'A') and value is
        the corresponding RNA sequence (e.g. 'ACGUGUCAUGGC')
    seg_order: Order of the segments of RNA
        e.g. ['A', 'D', 'B', 'C']
        seq: The full sequence of RNA from RNA_frag ordered by seg_order
    """

    # Choose a random seed if seed doesn't exist
    if args_d['seed'] is None:
        args_d['seed'] = np.random.randint(2**32-1) # 32 bit unsigned integer
    np.random.seed(args_d['seed'])
    print(f"Starting seed is {args_d['seed']}.")

    # 5' and 3' constant regions
    cprime = [args_d['constant_5_prime'], args_d['constant_3_prime']]

    RNA_frag, seg_order = generate_simple_seq(**args_d)
    RNA_seq = RNA_Sequence(RNA_frag, seg_order, cprime)

    args_d['segment_list'] = seg_order # save original segment order for future reference
    return RNA_seq