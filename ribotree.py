from MCTS_class import MCTS
from evaluate import ensemble_score_func
import ribotree_parser
import sys
import time
import os

def main():
	# Gets arguments
	args_d = ribotree_parser.get_args(argv=sys.argv[1:])

	# Print out information to inform the user what type of run is going on
	ribotree_parser.print_useful_info(args_d)

	# Create a unique ID for saving output results
	timestamp = time.strftime("%Y%m%d-%H%M%S")
	unique_id = f"{timestamp}-{os.getpid()}"
	args_d['unique_id'] = unique_id
	args_d['dna'] = False # check later

	# Check if sequence mode should be activated
	if args_d['sequence'] is not None:
		print('Activating sequence mode. Will not run MCTS.')
		constraint_binary, constraint_float, prob_list = ensemble_score_func(args_d, seq=args_d['seq'])
		print(args_d['sequence'], constraint_binary, constraint_float, prob_list)
	else:
		MCTS(args_d)
	print('Done')

if __name__ == '__main__':
	main()
