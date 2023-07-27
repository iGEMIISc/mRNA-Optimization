from functools import partial
import networkx as nx
import numpy as np
import pandas as pd
import random
from RNA_sequence_mutations import * 
from evaluate import ensemble_score_func, ensemble_score_comparison
from Tree import Tree
from plot_tree import plot_MC
from copy import copy
from save_results import save_results, save_node

RNA_BASES = set('ACGU')
DNA_BASES = set('ACGT')
#-----------------------------------------------------------------------------------

# functions for the Monte Carlo Tree Search Algorithm

def MCTS(args_d):
    """ Runs a MCTS to search for sequences

    Args:
        seed: RNG seed for starting point
        inc: Whether the RNA is to be increasing or decreasing w.r.t the
            reporter (R)
        n_iter: Number of iterations to perform
            DOES NOT EQUAL NUMBER OF NODES
            ~n_iter + 1 nodes if n_iter % n_children == 0
        n_children: Number of children for branching out
        constant: Balance between exploring depth vs breath
            0 -> depth
            inf -> breathd
        r_round: Specified which DMSSD round
        num_shuffle: Number of segments to shuffle when shufflings the segment 
            order
        shuffle_prob: The probability that the segments will be shuffled
            Units in decimal (0.5 == 50%)
        constant_length: Whether to preserve the length of each segment
        verbose: To print everything for debugging    
    """

    n_children = args_d['n_children']
    C = args_d['c_const']

    # generate starting point
    tree = Tree(args_d)
    root = tree.create_root()
    child, solution = tree.add_children()

    # save starting point
    MC_results = save_node(root, MC_results = [])
    save_results(args_d, root, root, MC_results, save_full=False)

    best_leaf = root # Keep track of best solution so far
    for i in range(args_d['n_iter']):
        # 1: Traversal/Selection
        # select leaf
        leaf = tree.select_leaf()

        # 2: Node expansion
        # if leaf has never been visited before rollout/simulation
        # else create children and then pick a random child for rollout/simulation
        leaf, solution = tree.expand_leaf(leaf)
        MC_results = save_node(leaf, MC_results)

        if solution:
            print('found solution on iteration {}'.format(i+1))
            return root, leaf

        # 3: Rollout/simulation
        tree.simulate_leaf(leaf)

        # 4: backpropagation based on results
        tree.back_propagate(leaf, C=C)

        # 5: beam prune
        if C > 0 and args_d['beam']:
            tree.beam_prune(beam = n_children+1)  # TODO: maybe allow for option of different beam sizes

        # Update best leaf found thus far
        # Lower dG is better
        if ensemble_score_comparison(best_leaf.get_score(), leaf.get_score(), **args_d)[0] > 1.0:
            best_leaf = leaf
            print('Current best node:', best_leaf.get_node_status())

        # Plot result and print result
        if args_d['plot'] and (i % n_children)==0:
            plot_MC(root)
        if (i+1) % args_d['stride'] == 0:
            print(f"Iteration: {i+1}. Saving.")
            save_results(args_d, root, best_leaf, MC_results, save_full=False)

    print("Max iterations reached.")
    save_results(args_d, root, best_leaf, MC_results, save_full=False)
    save_results(args_d, root, best_leaf, MC_results, save_full=True)
    return

def MC_move(RNA_seq_new, cur_id, args_d, cur_id_category, cur_id_options, preserve_list):
    # Mutate the segment using a strategy

    # Choose move
    move_list = ['shuffle', 'mutate']
    p_list = np.array([args_d['shuffle_prob'], 1 - args_d['shuffle_prob']])

    if (args_d['num_shuffle'] == 0) or (len(RNA_seq_new.seg_order) == 1):
        p_list[0] = 0
    if cur_id in preserve_list:
        p_list[1] = 0
    move = np.random.choice(move_list, p=p_list/np.sum(p_list))
    
    if move == 'shuffle':
        seg_order_new, move = shuffle(args_d['num_shuffle'], RNA_seq_new)
        if args_d['preserve_order'] is not None: # check that the shuffle doesn't mess up preserved_order
            check = False
            while not check:
                check = True
                for order_list in args_d['preserve_order']:
                    temp = [np.where(seg_order_new == x)[0] for x in order_list]
                    if np.any([len(x) == 0 for x in temp]):
                        raise ValueError('Missing segments specified in -preserve_order')
                    dx_list = np.diff([x[0] for x in temp])
                    if not (np.all(dx_list > 0) or np.all(dx_list < 0)):
                        check = False
                seg_order_new, move = shuffle(args_d['num_shuffle'], RNA_seq_new)

        RNA_seq_new.seg_order = seg_order_new

    elif move == 'mutate': # if random pick either strategy
        if np.random.randint(3) == 2 and cur_id_category not in  ['aa', 'mrna']: # 1/3 chance swap, 2/3 chance mutate (since 2 different mutations)
            RNA_seq_new, move = swap_fragment(RNA_seq_new, cur_id, args_d)

        else:
            RNA_seq_new, move = mutate_fragment(RNA_seq_new, cur_id, \
                args_d['constant_length'], args_d['preserve_seq'], \
                num_mutate=args_d['num_mutate'], DNA=args_d['dna'], \
                category=cur_id_category, options=cur_id_options, \
                CDSFold_path=args_d['CDSFold_path'], \
                CDSFold_prob=args_d['CDSFold_prob'])
    return RNA_seq_new, move

#------------------------------------------------------------------------------
# Functions for template/base based monte carlo search

def MC_base_segment(RNA_seq):
    """ returns a random segment to modify in the MC search 
    Will avoid choosing and segments in preserve_list
    """
    cur_id = np.random.choice(RNA_seq.seg_order)
    # while cur_id in preserve_list:
    #     cur_id = np.random.choice(RNA_seq.seg_order)
    cur_id_category = RNA_seq.get_category(cur_id)
    cur_id_options = RNA_seq.get_options(cur_id)
    return cur_id, cur_id_category, cur_id_options

def MC_search(RNA_seq, args_d, verbose=False):
    """ performs a monte carlo search around RNA_seq

    Args:
        RNA_seq: An instance of RNA_Sequence class (MUST GIVE A COPY RNA_seq.get_copy())
        verbose: To print everything for debugging
    Returns:
        RNA_seq_new: A new RNA_Sequence found by MC.
    """
    search_num = 1 # Keeps track of how many searches it takes
    start_seq = RNA_seq.get_seq()

    targeted_repeat_restriction_search = False

    while True:
        # Start with a copy so that you can keep the original
        RNA_seq_new = RNA_seq.get_copy()

        # Pick segment id to modify (e.g. 'R')
        if not targeted_repeat_restriction_search:
            cur_id, cur_id_category, cur_id_options = MC_base_segment(RNA_seq_new)

        if verbose:
            print('Chose segment: {}\n'.format(cur_id), 'Original: {}'.format(RNA_seq_new.RNA_frag[cur_id]))

        # Mutate the segment using a strategy
        RNA_seq_new, move = MC_move(RNA_seq_new, cur_id, args_d, cur_id_category, cur_id_options, args_d['preserve_seq'])
        if verbose:
            print(move)

        # Ensure the sequence is not the same
        seq = RNA_seq_new.get_seq()
        if seq == start_seq:
            continue

        # Check for stretches of nucleotide repeats or restriction sites
        if not args_d['allow_nucleotide_repeat']: # we want to prohibit nucleotide repeats
            repeat_location = RNA_seq_new.check_repeat()
            if repeat_location==-1: # There wasn't a repeat
                break #we're good!
            else:
                if verbose: print('Found a repeat')
                targeted_repeat_restriction_search=True
                cur_id_options = copy(repeat_location)

        else:
            break # we didn't care

        if args_d['restriction_sites'] is not None: # we have restriction sites to avoid
            restrict_location = RNA_seq_new.check_restriction_sites()
            if restrict_location==-1: # No taboo restriction sites
                break # we're good!
            else:
                if verbose: print('Found a restriction site')
                targeted_repeat_restriction_search=True
                cur_id_options = copy(repeat_location)

        else:
            break
 
        if verbose: print('searching again... %d repeat %s restriction %s'%\
         (search_num, RNA_seq.check_repeat(), RNA_seq.check_restriction_sites()))
        search_num += 1

    if verbose:
        print([RNA_seq_new.RNA_frag[x] for x in RNA_seq_new.seg_order], \
            RNA_seq_new.seg_order, 'chosen segment: {}'.format(cur_id))
    RNA_seq_new.reset_evaluation()
    return RNA_seq_new, move

def compile_MC_result(RNA_seq, ensemble_score, accept, move = 'start', MC_result = []):
    """ compiles the MC result numpy array as described in MC_iteration """
    MC_result_temp = [accept] + [RNA_seq.RNA_frag[x] for x in RNA_seq.seg_order] + list(RNA_seq.seg_order)
    for score in ensemble_score: # iterate through binary, float, and prob scores
        for condition in score: # iterate through each condition
            for value in condition: # append each value in each condition
                MC_result_temp.append(value)
    MC_result_temp.append(move)
    MC_result_temp = [MC_result_temp]

    if len(MC_result) > 0:
        return MC_result + MC_result_temp
    else:
        return MC_result_temp

def MC_iteration(RNA_seq, n_iter, verbose=False):
    """ performs a Monte Carlo (MC) search around the given sequence 
    
    Args:
        RNA_seq: An instance of RNA_Sequence class
            Stores all the Monte Carlo parameters as well
        ensemble_score: Ensemble score of the RNA_seq sequence that the MC search
            is starting wtih
        n_ter: Number of MC iterations to perform
        verbose: To print everything for debugging

    Returns:
        MC_result: The trajectory of the MC simulation as a numpy array
            each row comtains the following:
            [RNA_frag (str), seg_order (str), constraint_binary, constraint_float, raw_constraint]
            note: raw_constraint will be flattened to a 1-d array
        enesemble_score: The most recent accept score
    """
    MC_result = np.array([], dtype=object)
    args_d = RNA_seq.get_args_d()
    
    T = args_d['T']
    counter = 0
    num_failures = 0 # Keep track of failures

    seq_set = set() # Keep track of all sequences generated to prevent repeating
    # Start MC searching
    while counter < n_iter:

        # Perform a MC search from the RNA_seq
        RNA_seq_new, move = MC_search(RNA_seq, args_d, verbose=verbose)

        if RNA_seq_new.get_seq() in seq_set:
            # TODO: if n_iter is set to 1, what happens here???
            if verbose: print('Already checked this sequence.')
            continue
        else:
            seq_set.add(RNA_seq_new.get_seq())

        # Calculate the score of the new sequence
        ensemble_score_new = ensemble_score_func(args_d, RNA_seq = RNA_seq_new)
        RNA_seq_new.store_score(ensemble_score_new)

        # Accept or Reject new result
        p, solution_found = ensemble_score_comparison(RNA_seq.get_score(), RNA_seq_new.get_score(), T)
        if verbose: print('Accept Probability: {0:.2f}%'.format(p*100))

        accept = 0
        if (np.random.rand() <= p): # Even if worse, may accept
            accept = 1
            RNA_seq = RNA_seq_new
            counter += 1
            # reset values
            num_failures = 0
            T = args_d['T']
            if verbose: print('Accept')

        else:
            num_failures += 1 
            # Maybe add simulated annealing to T.
            if verbose: print('Reject')

        # Concatenate results
        MC_result = compile_MC_result(RNA_seq_new, RNA_seq_new.get_score(), accept, move = move, MC_result = MC_result)
        #MC_result_df = compile_MC_result_v2(RNA_seq_new, RNA_seq_new.get_score(), accept, move = move, MC_result = MC_result_df)

        if accept:
            if verbose: 
                print("{:<5}".format(counter), RNA_seq_new.get_seq(), RNA_seq_new.get_score()[2][0][0])

        else:
            if verbose:
                print("\n{:<5}".format('Rejected'), RNA_seq_new.get_seq(),\
                    RNA_seq_new.seg_order, RNA_seq_new.get_score()[0], RNA_seq_new.get_score()[1], RNA_seq_new.get_score()[2], move)

        # If solution is found break early
        if solution_found:
            RNA_seq.set_solution()
            break
    return RNA_seq, MC_result