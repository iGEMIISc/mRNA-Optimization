import numpy as np
import pandas as pd
from evaluate import *


rt = 1.987e-3*310 # for calcualting probabilities from delta G
#-----------------------------------------------------------------------------------

class Node(object):
    def __init__(self, RNA_seq, MC_result, layer = -1):
        self.RNA_seq = RNA_seq
        self.MC_result = MC_result
        self.ensemble_score = RNA_seq.get_score()
        self.parent = None
        self.children = []
        self.total_value = 0
        self.n_visits = 0
        self.UCT1 = np.inf
        self.dead = False
        self.solution = RNA_seq.get_solution()
        self.layer = layer

    def check_solution_found(self):
        return self.solution

    def run_simulation(self):
		# TODO: ????
        # total value is -dG since I want to maxmize this
        # (-inf == all conditions true)
        dG = ensemble_score_dG(self)
        self.total_value = np.exp(-dG/rt)
        return self.check_solution_found()

    def set_parent(self, node):
        self.parent = node

    def set_id(self, node_id):
        self.id = node_id

    def get_id(self):
        return self.id

    def get_layer(self):
        return self.layer

    def set_children(self, node):
        self.children.append(node)

    def set_dead(self):
        self.dead = True
        parent = self.parent
        while parent is not None:
            # If all are dead then set parent to dead and move up, else
            if all([child.get_dead() for child in parent.get_children()]):
                parent.dead = True
                parent = parent.get_parent()
            else:
                break

    def set_UCT1(self, value):
        self.UCT1 = value

    def get_children(self):
        return self.children

    def get_dead(self):
        return self.dead

    def get_score(self):
        return self.ensemble_score

    def get_parent(self):
        return self.parent

    def get_MC_result(self):
        return self.MC_result

    def get_RNA_seq(self):
        return self.RNA_seq

    def get_UCT1(self, check_children = False, dG = False):
		# TODO: move the actual founctionality to Tree, simplify to getter
        if check_children:
            UCT1 = -np.inf
            for child in self.children:
                if child.get_dead(): continue # ignore dead children
                grandchildren = child.get_children()
                if len(grandchildren) == 0:
                    UCT1 = max(child.get_UCT1(dG=dG), UCT1)
                else:
                    UCT1 = max(child.get_UCT1(dG=dG), UCT1)
                    child.get_UCT1(check_children=True)
            return UCT1
        if dG:
            dG_temp = ensemble_score_dG(self)
            return np.exp(dG_temp/rt)
        else:
            return self.UCT1

    def get_n_visits(self):
        return self.n_visits

    def get_total_value(self):
        return self.total_value

    def get_node_status(self):
        return self.RNA_seq.get_status()

    def update_UCT1(self, C = 2):
		# TODO: combine with get_UCT1?
        n = self.n_visits
        N = self.parent.n_visits
        if C >= 0:
            UCT1 = self.total_value/float(n) + C*np.sqrt(np.log(N)/n)
        else: # If negative C, do a linear MCTS
            UCT1 = self.total_value
        self.UCT1 = UCT1

    def visit(self, value):
        self.n_visits += 1
        self.total_value += value


