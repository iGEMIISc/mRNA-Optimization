import argparse
import re
from evaluate import ensemble_score_func
from utils_parser import *

def print_useful_info(args):
    """ Prints useful info to terminal """

    # Print out what mode is being run
    if args['c_const'] <= 0:
        mode = 'linear Monte Carlo'
    else:
        if args['beam']:
            mode = 'beam MCTS'
        else:
            mode = 'MCTS'
    print(f'Running {mode}.\n')

    # Print all args
    print('Printing all arguments:')
    for key, value in args.items():
        print(f'{key}:{value}')

def organize_input(args):
    """ Fixes arguments that depends on other arguments """

    # Fix --domain_list
    # If domain_list is missing assume all domains are in the sequence
    if args['domain_list'] is None:
        args['domain_list'] = list(args['domain'].keys())

    # Fix --mrna
    if any([v.get_category() in ['mrna', 'aa'] for k,v in args['domain'].items()]):
        args['mrna'] = True

    # Fix --condition if the value is blank
    for k,v in args['condition'].items():
        if type(v) == str:
            args['condition'][k] = [[v, None]]
        elif type(v[0]) == str:
            args['condition'][k] = [[v[0], None]]
    return args

def check_valid(args):
    """ Checks if the arguments are internally consistent """
    # 1. Check if all RNA or DNA. Also checks if consistent with --dna flag
    seq = ''
    for k, v in args['domain'].items():
        if v.get_category() not in ['aa']: 
            seq += v.get_og_sequence()

    if 'T' in seq:
        raise ValueError("Must use RNA bases for ribotree-mrna.")

    if args['linearfold'] == True:
        if args['package'] == 'nupack':
            raise ValueError("LinearFold only works with contrafold, vienna, and eternafold. Does not work with nupack.")

    return True

def get_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--file', type=open, action=LoadFromFile,
        help="Text file where flags are stored.")
    parser.add_argument('--input_file', type=str,
        help=argparse.SUPPRESS) # Hidden flag. Only for code use. Stores file name from --file
    parser.add_argument('--output', type=parse_output_flag, default='',
        help="Specifies the output location to save results.")

    # Parameters to define sequence/domains
    parser.add_argument('--domain_list', type=parse_domain_list_flag,# Required
        help="""Specifies domains inside the RNA using the domain name. E.g. \
        -domain_list G H""")
    parser.add_argument('--preserve_seq', type=parse_flag_brackets, default=[],
        help='List of segments to never mutate. E.g. -perserve_seq A B C')
    parser.add_argument('--preserve_order', type=parse_flag_brackets, default=None,
        help="""List of segment orders preserver during segment swaps. This is \
        important if you have an RNA structure made up of more than 2 pieces. \
        E.g. If you have a an aptamers that made of three RNA pieces then use: \
        --preserve_order [A B C] [D E F]""")
    parser.add_argument('--domain_order', type=str2bool, nargs='?', const=True, default=False,
        help="""Whether to use the order in segment_list as the initial RNA order.""")
    parser.add_argument('--constant_5_prime', type=parse_constant_prime_flag, default='',
        help="""Appends the RNA sequence to the 5 prime end. Becareful of four \
        nucleotide repeats since the default mode is to not allow those during the \
        MCTS. E.g. If you wanted to add a barcode to the 5 prime end \
        -constant_5_prime GGAUUCUA""")
    parser.add_argument('--constant_3_prime', type=parse_constant_prime_flag, default='',
        help="""Appends the RNA sequence to the 5 prime end. Becareful of four \
        nucleotide repeats since the default mode is to not allow those during the \
        MCTS. E.g. If you wanted to add a ribosome binding side to the 3 prime end \
        -constant_3_prime AGGAGG""")
    parser.add_argument('--domain', type=parse_domain_flag, action=dictionary_append, nargs='+', \
        help="""Domains that make up the RNA. See README for details.""")

    # Parameters for MCTS
    parser.add_argument('--length', type=int, default=0,
        help="Length of RNA to generate. Leave as 0 for mRNA design.")
    parser.add_argument('--n_iter', type=int, default=300,
        help="Number of iterations to run the MCTS. Default is 300.")
    parser.add_argument('--stride', type=int, default=100,
        help="Number of iterations to save full statistics on the run. Default is 100.")
    parser.add_argument('--n_children', type=int, default=3,
        help="Number of children to generate from each parent node. Default is 3")
    parser.add_argument('--c_const', type=float, default=-1,
        help="""Exploration constant in Upper Confidence Bound 1. Large c favors \
        exploration. c is a float from [0, inf]. If c is < 0 then the code will \
        instead run a linear MCT. Default is -1 (defaults to linear MCTS)""")
    parser.add_argument('--beam', type=str2bool, nargs='?', const=True, default=False,
        help="Whether to do a beam MCTS or not. Beam is set to n_children + 1.")
    
    # Parameters for MCTS moves/criteria
    parser.add_argument('--num_shuffle', type=int, default=1,
        help="""Number of segments to shuffle for each shuffle move. During a shuffle \
        move N segments are randomy removed from the segment_list. They are randomly \
        added back and checked to ensure the segment order has changed. E.g. If \
        shuffle is 1 then [A B C] -> [A B] (remove C) -> [A C B] or [C A B] (add C)""")
    parser.add_argument('--num_mutate', type=int, default=1,
        help="""Number of mutations to make at once. May be useful for long sequences.""")
    parser.add_argument('--shuffle_prob', type=restricted_float, default=0.33,
        help="Probability of shuffling segments. Float from 0 to 1. Default is 0.33")
    parser.add_argument('--constant_length', type=str2bool, nargs='?', const=True, default=False,
        help="Whether to keep the segments constant length during the MCTS. There is \
        a move that will remove an edge base and insert a random base into another \
        segment. By default this move is allowed. Turning this flag on will stop that \
        move from occuring.")

    # CDSFold parameters
    parser.add_argument('--CDSFold_path', type=str, default=None,
        help="""Probability of performing CDSFold as a move. Default is 0.""")
    parser.add_argument('--CDSFold_prob', type=restricted_float, default=0,
        help="""Probability of performing CDSFold as a move. Default is 0. Must be a \
            float from [0 1].""")

    #Sequence parameters
    parser.add_argument('--gu_level', type=restricted_float, default=0, 
        help="""Probability of GU wobble base pair when generating complementary \
        RNA. This only matters when generating the initial sequence for running MCTS. \
        Must be a float from [0 1] where 0 indicates NO GU base pairs.""")
    parser.add_argument('--allow_nucleotide_repeat', type=str2bool, nargs='?', const=True, default=False,
        help="""Whether to restrict the following nucleotide repeats AAAAA, CCCC, \
        UUUU, and GGGG. By default repeats are not allowed since sequences may \
        be difficult to synthesize.""")

    # Evaluation parameters
    parser.add_argument('--condition', type=parse_flag_brackets, action=dictionary_append, nargs='+', \
        help="""Conditions for evaluating switch performance.""")
    parser.add_argument('--package', type=parse_package_flag, default='contrafold',
        help="""Package to use to calculate metrics. Choose between 'contrafold', 'vienna', and 'eternafold'""")
    parser.add_argument('--linearfold', type=str2bool, nargs='?', const=True, default=False,
        help="""Whether to use LinearFold mode when using the folding package.""")
    parser.add_argument('--sequence', type=str, default=None,
        help="RNA/DNA sequence for testing against all conditions. DOES NOT run MCTS. \
            ONLY EVALUATES the sequence and ends.")
    parser.add_argument('--mod_U', type=str2bool, nargs='?', const=True, default=False,
        help="""Whether to mask uracils when running degscore.""")

    parser.add_argument('--average', type=parse_average_flag, default=None,
        help="""Specifies if you want to average the results of high and low signal. \
        E.g. -average DR:1""")
    parser.add_argument('--Temp', type=float, default=37,
        help="Temperature for predictions in celcius. Default is 37 C. Only affects output of folding packages that accepts temperature.")

    parser.add_argument('--T', type=float, default=330.15,
        help="Temperature for probability of accepting moves in Kelvin. Will be multiplied by gas constant R.")
    parser.add_argument('--scale', type=parse_scale_flag, default=1,
        help="Scaling factor for evaluation. Max is 400.")
    parser.add_argument('--mrna', type=str2bool, nargs='?', const=True, default=False,
        help="Whether to do mRNA mode for minimizing delta G.")
    parser.add_argument('--ratio', type=str2bool, nargs='?', const=True, default=False,
        help="Whether to do ratio mode. Only use when specifying two concentrations in --condition.")
    parser.add_argument('--restriction_sites', type=parse_restriction_sites,
        help='Restriction sites to avoid. Will also avoid reverse complement of sites listed here.')

    # Visualization/print/etc parameters
    parser.add_argument('--verbose', type=str2bool, nargs='?', const=True, default=False,
        help="Boolean to print everything during run.")
    parser.add_argument('--plot', type=str2bool, nargs='?', const=True, default=False,
        help="Boolean to display plots during run.")
    parser.add_argument('--save_ss', type=str2bool, nargs='?', const=True, default=False,
        help="Boolean to save the results of nupack for each sequence.")
    parser.add_argument('--seed', type=int, default=None,
        help="Seed for random number generator.")

    return parser

def get_args(argv=None):
    parser = get_parser()

    # Create dictionary of arguments from parser
    if argv is not None:
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()
    args = vars(args)
    
    print('NOTE: Must use "," as delimiter for values inside [...]. E.g. [A,B,C]')

    # Correct some values of flags depending on other flags
    args = organize_input(args)

    check_valid(args)
    return args