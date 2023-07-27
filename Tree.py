from Node import *
from collections import defaultdict
from save_results import *
from init_starting_RNA_Sequence import gen_RNA_start_simple
from utils_sequence import match_seq_len
import MCTS_class
import multiprocessing as mp
from functools import partial


class Tree(object):
    def __init__(self, args_d):
        self.n_children = args_d['n_children']
        self.args_d = args_d
        self.root = None
        self.layer_dict = defaultdict(lambda : set())


    # Step 0
    def create_root(self):
        # TODO: move some of this where it belongs (RNA_sequence.py?)
        RNA_seq = gen_RNA_start_simple(self.args_d)
        RNA_seq.store_seed(self.args_d['seed'])
        RNA_seq.store_args_d(self.args_d)

        ensemble_score = ensemble_score_func(self.args_d, RNA_seq = RNA_seq)

        # store score and parameters
        RNA_seq.store_score(ensemble_score)
        MC_result = MCTS_class.compile_MC_result(RNA_seq, ensemble_score, 1)
    
        print("Creating root")
        # print('root ensemble score', ensemble_score)
        # print('root info', "\n{:<5}({:3})".format(0, 1), RNA_seq.get_seq(),RNA_seq.seg_order, ensemble_score[0], ensemble_score[1], ensemble_score[2])

        root = Node(RNA_seq, MC_result, layer = 0)
        root.set_id("R")
        self.root = root
        return root


    # Step 1
    def select_leaf(self):
        """ Recursively selects the child with the max value until
            a leaf is reached.
            Returns:
                - a leaf node with a locally maximal value.
        """
        cur_child = self.root
        cur_children = cur_child.get_children()
        while len(cur_children) > 0:
            cur_child = _select_max_child(cur_children)
            cur_children = cur_child.get_children()
        return cur_child


    # Step 2
    def expand_leaf(self, node):
        """ If the leaf node has been visited before, generate children
            for the node and return one of them; else, return input node.
            Args:
                - node: leaf node to either add children to or return
            Returns:
                - either the input node or one of its new children
                - boolean, was a solution found
        """
        if node.get_n_visits() > 0:
            node, solution = self.add_children(node)
            return node, solution
        else:
            return node, False


    # Step 3
    def simulate_leaf(self, node):
        node.run_simulation()


    # Step 4
    def back_propagate(self, node, C = 2):
        if node.get_parent() is None:
            raise ValueError('Somehow parent node was lost')
        
        cur_node = node
        value = node.get_total_value()
        cur_node.visit(value)

        cur_parent = cur_node.get_parent()
        node_traversed = []
        while cur_parent is not None:
            node_traversed.append(cur_node)
            cur_parent.visit(value) # Add to the total value & visit of the parent
            cur_node = cur_parent
            cur_parent = cur_node.get_parent()

        # Update going down since the visits are increased going up  # TODO: ???
        for node in node_traversed:
            node.update_UCT1(C=C)


    def beam_prune(self, beam = np.inf):
        # TODO: check correctness of implementation

        # Check each layer to ensure that there are only beam # alive
        layer_list = sorted(self.layer_dict.keys())
        last_layer = layer_list[-1]
        for layer in layer_list:
            node_list = self.layer_dict[layer]
            if layer != last_layer:
                for node in node_list:
                    if not node.get_dead() and len(node.get_children()) == 0 and len(self.layer_dict[layer+1]) > beam:
                        node.set_dead()
            num_alive = sum([not node.get_dead() for node in node_list])
            UCT1_list = [node.get_UCT1(check_children=True, dG=True) for node in node_list]

            # Sort node_list and remove if not top N
            node_list = [node for _,node in sorted(zip(UCT1_list,node_list), key=lambda x:x[0])]

            if num_alive > beam:
                node_dead = node_list[:-beam]
                for node in node_dead:
                    node.set_dead()


    def add_children(self, parent_node = None):
        """ Adds n_children to the input leaf parent node, returns one child.
            Args:
                - parent_node: the node to create children for.
                - verbose: boolean, whether or not to print extra info.
            Returns:
                - a solution child node, if one is found; else, a child node
                    selected at random.
                - boolean, was a solution node found
            Note:
                - Bottleneck in MCTS due to calculating RNA folding energy
        """
        if parent_node is None:
            parent_node = self.root

        new_children = _create_children(parent_node, self.n_children, self.args_d['verbose'])
        self._add_children_to_tree(new_children, parent_node)

        for child_node in new_children:
            if child_node.check_solution_found():
                return child_node, True
        return np.random.choice(new_children), False


    def _add_children_to_tree(self, children_nodes, parent_node):
        for child_node in children_nodes:
            child_num = len(parent_node.get_children()) + 1
            child_node.set_id(parent_node.get_id() + '.' + str(child_num + 1))
            child_node.set_parent(parent_node)
            parent_node.set_children(child_node)
            self.layer_dict[parent_node.get_layer() + 1].add(child_node)


def _create_children(parent_node, n_children, verbose = True):
	# first, check that this node does not already have children
    assert(len(parent_node.get_children()) == 0)

    # TODO: use apply_async instead?
    pool = mp.Pool()

    # Do not send in a tree objects into multiprocess.
    parent_sequence = parent_node.get_RNA_seq()
    parent_layer = parent_node.get_layer()
    new_children = pool.map(partial(_create_child, parent_sequence=parent_sequence, \
        parent_layer=parent_layer, verbose = verbose), range(n_children))
    pool.close()
    pool.join()

    # check that the multiprocessing has correctly finished
    assert len(new_children) == n_children
    return new_children


def _create_child(child_num, parent_sequence, parent_layer, verbose = True):
    np.random.seed()  # TODO: ???
    RNA_seq, MC_result = MCTS_class.MC_iteration(parent_sequence, 1, verbose=verbose) # 1 iteration
    child = Node(RNA_seq, MC_result, parent_layer + 1)
    return child


def _select_max_child(nodes):
    """ Selects the node with greatest UCT1 score from a list of nodes. """
    live_nodes = [node for node in nodes if not node.get_dead()]
    live_nodes_UCT1s = [node.get_UCT1() for node in live_nodes]
    max_UCT1_index = live_nodes_UCT1s.index(max(live_nodes_UCT1s))
    return live_nodes[max_UCT1_index]