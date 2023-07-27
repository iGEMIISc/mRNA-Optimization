import copy
import utils_domain as ud

RNA_BASES = set('ACGU')
DNA_BASES = set('ACGT')

class Domain(object):
    """ object to the domains in the RNA sequence """
    def __init__(self, name, sequence, category, options):
        # Clean inputs and check if valid
        sequence = sequence.upper()
        category = category.lower()
        options = options.lower()

        valid_categories = ['aa','mrna','rna','dna']
        assert category in valid_categories, f"Domain category must be on of the following: {valid_categories}" 

        # Save values
        self.name = name
        self.og_sequence = sequence # store original sequence for future reference
        self.sequence = None # Not set until generate_new_sequence() is run
        self.og_length = ud._process_length(sequence, "F", category, options) # assumes FULL LENGTH
        self.length = self.og_length
        self.ss = None

        self.options = options
        self.category = category

    def __repr__(self):
        return f"<Domain name:{self.name}, sequence:{self.og_sequence}, ss:{self.ss}, length:{self.og_length}, category:{self.category}, options:{self.options}>"

    def check_dependent(self):
        """ Checks if the domain depends on another domain """
        if self.category == 'mrna':
            return False
        elif self.category == 'aa':
            return False
        elif 'g' in self.options: # g is for generate
            return False
        elif len(self.options) == 0:
            return False
        else:
            return True

    def check_processed(self):
        """ Checks if domain has been fully processed """
        if self.sequence is None:
            return False
        else:
            return True

    def generate_init_sequence(self, parent_sequence=None, allow_nucleotide_repeat=False):
        self.sequence = ud._generate_sequence(self.og_sequence, self.length, \
            self.category, self.options, parent_sequence=parent_sequence, allow_nucleotide_repeat=allow_nucleotide_repeat)

    def get_copy(self):
        """ returns a deep copy of class instance """
        return copy.deepcopy(self)

    def get_name(self):
        return self.name

    def get_sequence(self):
        return self.sequence

    def get_og_sequence(self):
        return self.og_sequence

    def get_ss(self):
        if self.ss is None:
            return '.'*len(self.sequence)
        else:
            return self.ss

    def get_og_length(self):
        return self.og_length

    def get_length(self):
        return self.length

    def get_category(self):
        return self.category

    def get_options(self):
        return self.options

    def get_moves(self):
        raise ValueError("Not coded yet.")

    def set_name(self, name):
        self.name = name

    def set_length(self, length):
        self.length = length

    def set_sequence(self, sequence):
        self.sequence = sequence

    def reset(self):
        self.sequence = None