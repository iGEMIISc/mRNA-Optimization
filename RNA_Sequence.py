import copy
import numpy as np
import utils_sequence as us
from arnie.bpps import bpps
from arnie.mfe import mfe
from RiboGraphViz import RGV
from utils_features import calc_cai

class RNA_Sequence(object):
    """ object to represent the RNA sequence for monte carlo searching

    Args:
    RNA_frag: Dictionary where key is segment id (e.g. 'A') and value is
        the corresponding RNA sequence (e.g. 'ACGUGUCAUGGC')
       seg_order: Order of the segments of RNA
           e.g. ['A', 'D', 'B', 'C']
    """
    def __init__(self, RNA_frag, seg_order, cprime):
        self.RNA_frag = copy.copy(RNA_frag)
        self.seg_order = np.copy(seg_order)

        # Store constant 5' and 3' constraints
        self.c5prime = cprime[0]
        self.c3prime = cprime[1]

        # Store original values for future reference
        self.seg_order_og = np.copy(self.seg_order)
        self.seg_len_og = {}
        for x in self.seg_order:
            self.seg_len_og[x] = len(RNA_frag[x])
        self.solution = False

        self.args_d = None
        self.score = None
        self.bpm = None
        self.mfe = None
        self.dG_MFE = None
        self.mfe_ss_constraint = None
        self.RGV = None

    def check_ratio_mode(self):
        return self.args_d['ratio']

    def check_mrna_mode(self):
        return self.args_d['mrna']

    def check_repeat(self):
        """ Boolean if a 5 nucleotide repeat exists IGNORING constant regions """
        seq = ''.join([self.RNA_frag[x] for x in self.seg_order])
        return us.check_nucleotide_repeat(seq)

    def check_restriction_sites(self):
        """ Boolean if a 5 nucleotide repeat exists IGNORING constant regions """
        seq = ''.join([self.RNA_frag[x] for x in self.seg_order])
        return us.check_restriction_sites(seq, self.args_d['restriction_sites'])

    def get_args_d(self):
        return self.args_d

    def get_bpm(self):
        # Assumes package, param_file, linear doesn't change
        if self.bpm is None:
            dna = self.args_d['dna']
            if dna:
                param_file = 'DNA'
            else:
                param_file = None
            self.bpm = bpps(self.get_seq(), package=self.args_d['package'], 
                param_file=param_file, 
                linear=self.args_d['linearfold'])
        return self.bpm

    def get_category(self, domain_name):
        return self.args_d['domain'][domain_name].get_category()

    def get_RGV(self):
        #if self.RGV is None:
        #    self.RGV = RGV(self.get_mfe())
        self.RGV = RGV(self.get_mfe())
        self.RGV.run_structure_properties()

        return self.RGV

    def get_dG_MFE(self):
        if self.dG_MFE is None:
            self.mfe, self.dG_MFE = mfe(self.get_seq(), package=self.args_d['package'], linear=True, return_dG_MFE = True)
        return self.dG_MFE

    def get_mfe(self, ss_constraint=None):
        if self.mfe is None or (self.mfe_ss_constraint != ss_constraint and ss_constraint is not None):

            self.mfe, self.dG_MFE = mfe(self.get_seq(), package=self.args_d['package'], linear=True, return_dG_MFE = True)

            self.mfe_ss_constraint = ss_constraint
        return self.mfe

    def get_options(self, domain_name):
        return self.args_d['domain'][domain_name].get_options()

    def get_scale(self):
        return self.args_d['scale']

    def get_copy(self):
        """ returns a deep copy of class instance """
        return copy.deepcopy(self)

    def get_score(self):
        return self.score

    def get_segment(self, seg_id):
        return self.RNA_frag[seg_id]

    def get_segment_bounds(self, seg_id):
        """ Returns the index bounds for the segment in the full sequence """
        length_list = [len(self.RNA_frag[x]) for x in self.seg_order]
        idx = np.where(self.seg_order == seg_id)[0][0]
        start_idx = len(self.c5prime)
        start_idx += int(sum(length_list[:idx]))
        end_idx = start_idx + len(self.RNA_frag[seg_id])
        return start_idx, end_idx

    def get_ss_open(self, rbs_length=14):
        """ Returns the secondary structure of an open mRNA """
        total_length = len(self.get_seq())
        c5p_length = len(self.c5prime)
        c3p_length = len(self.c3prime)

        ss = c5p_length*['x']
        ss += min(rbs_length, total_length - c5p_length - c3p_length)*['x']
        ss += max(0, total_length - c5p_length - c3p_length - rbs_length)*['.']
        ss += c3p_length*['.']
        ss = ''.join(ss)
        return ss

    def get_ss_con(self):
        """ Returns secondary structure of constrained structure """
        ss = '.'*len(self.c5prime)
        for domain in self.seg_order:
            ss += self.args_d[domain].get_ss()
        ss += '.'*len(self.c3prime)
        return ss

    def get_seq(self):
        """ returns the RNA sequence """
        seq = ''.join([self.RNA_frag[x] for x in self.seg_order])
        seq = self.c5prime + seq + self.c3prime
        return seq

    # def get_precompiled_metrics(self):
    #     dct = {'sequence': self.get_seq(), 'CAI': calc_cai(self.get_seq())}

    #     if self.mfe is not None:
    #         dct.update({'mfe': self.mfe, 'dG_MFE': self.dG_MFE})

    #     if self.bpm is not None:
    #         punp_vec = 1 - np.sum(self.bpm, axis=0)
    #         dct.update({'AUP': np.mean(punp_vec), 'AUP_init14': np.mean(punp_vec[:14])})

    #     return dct

    def get_solution(self):
        return self.solution

    def get_status(self):
        return "%s %.3f" % (self.get_seq(), self.get_score()[-1][0][0])

    def store_seed(self, seed):
        self.seed = seed

    def store_args_d(self, args_d):
        self.args_d = args_d

    def store_score(self, score):
        self.score = score

    def set_solution(self):
        self.solution = True

    def reset_evaluation(self):
        self.score = None
        self.bpm = None
        self.mfe = None
        self.mfe_ss_constraint = None
