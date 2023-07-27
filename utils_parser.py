import argparse
import numpy as np
import pickle
import re
import utils_sequence as us
from domain import Domain

#-----------------------------------------------------------------------------------
# Argparse custom actions

class LoadFromFile (argparse.Action):
    """ Read flags from a file """
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # Copy the filename for storage as 'input_file'
            parser.parse_args(f'--input_file {f.name}'.split(), namespace)

            text = f.read().split('\n')
            for line in text:
                line = line.strip() # Remove excess whitespace
                line = line.split('#')[0]
                # Skip empty lines or lines without a flag
                if len(line) == 0 or '-' not in line: continue
                
                # convert all -flags to --flags
                if '--' not in line: line = '-' + line
                
                # split line into flag and argument
                # only split by first space
                line = [x.strip() for x in line.split(' ', 1)] 

                args, unknown_args = parser.parse_known_args(line, namespace)
                if len(unknown_args) > 0: # Unkown flag
                    raise argparse.ArgumentTypeError(f"{unknown_args} is not a valid flag.")

def values2kv(x):
    """ splits input into a key and value for dictionary storage later 
    Helper function for class dictionary_append
    """
    k = x[0]
    v = x[1:]
    
    # Clean up k
    if isinstance(k, list):
        if len(k) != 1:
            raise argparse.ArgumentTypeError(f"{self.dest} flag is incorrect. Check brackets.")
        k = k[0]

    # Clean up v
    if isinstance(v, list):
        if len(v) == 1:
            v = v[0]
    return k,v

class dictionary_append(argparse.Action):
    """ Converts the flag value into a dictionary """
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(dictionary_append, self).__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        previous_d = getattr(namespace, self.dest)
        k, v = values2kv(values)
        if previous_d is not None: # Append to dictionary
            previous_d[k] = v
            setattr(namespace, self.dest, previous_d)
        else: # Create new dictionary if it does not exist
            args = {k:v}
            setattr(namespace, self.dest, args)

#-----------------------------------------------------------------------------------
# Parser helper functions

def str2bool(x):
    """ Converts a string to a boolean """
    if isinstance(x, bool):
        return x
    if x.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif x.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f"{x} is not a valid boolean. True or False only.")

def restricted_float(x):
    """ Converts to a float between 0 and 1.0 """
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} is not a floating-point literal")
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} is not in range [0.0, 1.0]")
    return x

def to_int(x):
    """ Converts to an integer """
    try:
        x = int(x)
    except ValueError:
        pass
    return x

def to_float(x):
    """ Converts to a float """
    try:
        x = float(x)
    except ValueError:
        pass 
    return x

def split_value(x):
    """ Splits a string containing ':' """
    x_list = x.split(':')
    if len(x_list) != 2:
        raise argparse.ArgumentTypeError(f"{x} is not a valid use of ':'")
    x_list[1] = to_float(x_list[1])
    return x_list

def parse_flag_brackets(flag, allow_nested=False):
    """ Parses a string containing brackets and maintains organization """
    # Clean up flag
    flag = re.sub(r'[;]',',',flag) # Get rid of semicolons
    flag = re.sub(r'[\t ]+:[\t ]+', r':',flag) # Fixes spacing around colons
    flag = re.sub(r'\s+', ',', flag) # Replace all whitespace with ','
    flag = re.sub(r'([,])\1+', r'\1', flag) # Remove contiguous stretches of ','

    # Now convert to a list
    flag = flag.replace('],[','][')
    flag = re.split('\[|\]',flag)
    flag = [x for x in flag if len(x)]

    # Iterate and parse
    for idx, x in enumerate(flag):
        x = x.split(',')
        x = [v for v in x if isinstance(v, str) and len(v)]

        # Convert to float or split if ':' inside value
        flag[idx] = [ to_float(v) if ':' not in v else split_value(v) for v in x]
    
    # Get rid of nested lists
    if not allow_nested:
        if isinstance(flag, list) and len(flag) == 1:
            flag = flag[0]
    return flag

def parse_input_flag(x):
    x = parse_flag_brackets(x)
    return np.array(x)

def parse_domain_list_flag(x):
    """ Parses the segment list into a of domains and lengths and categories """
    d_list = parse_flag_brackets(x)

    if len(d_list) != len(set(d_list)):
        raise ValueError("{x} is not a valid domain_list. Each item in list must be unique.")
    return d_list

def parse_constant_prime_flag(x):
    """ Returns uppercase string of the bases """
    x = x.upper()
    regex_search_bases = re.compile(r'[^ACGTU]').search(x)
    valid = not bool(regex_search_bases)
    if valid and not ('U' in x and 'T' in x):
        return x.upper()
    else:
        raise argparse.ArgumentTypeError(f"""{x} contains invalid nucleotides. Must be \
            RNA or DNA bases.""")

def parse_domain_flag(x):
    """ Formats the domain inputs if using amino acids """

    x = parse_flag_brackets(x)
    if len(x) < 3 or len(x) > 4:
        raise ValueError("Domain flag consists of name, sequence, category, options. {x} is invalid")
    # Parse to create Domain object
    name = x[0]
    
    sequence = x[1].upper()
    
    category = x[2]
    
    if len(x) == 4:
        options = x[3]
    else:
        options = ""

    # convert DNA to RNA
    if category == 'mrna':
        sequence = sequence.replace('T','U')

    x = [name, Domain(name, sequence, category, options)]
    return x

def parse_package_flag(x):
    """ Returns lowercase string """
    x = x.lower()
    if x in ['contrafold', 'vienna', 'nupack', 'eternafold']:
        return x
    else:
        raise argparse.ArgumentTypeError(f"{x} is not a valid RNA/DNA package.")

def parse_output_flag(folder):
    """ Formats a folder with '/' if missing """
    if len(folder) > 0:
        try:
            if folder[-1] != '/':
                folder = folder + '/'
        except ValueError:
            raise argparse.ArgumentTypeError(f"{folder} is not a valid folder")
    return folder

def parse_average_flag(x):
    """ Parses the average flag """
    x =  "".join(x.split()) # Remove all spaces in the string
    if ':' in x:
        temp = x.split(':')
        temp[1] = float(temp[1])

        # check if valid --avarage flag
        # must be (w)AR/DR/ARDR
        metric = temp[0]
        if metric[0] == 'w': # remove 'w' at the start
            metric = metric[1:]
        if metric not in ['AR', 'DR', 'ARDR']:
            raise argparse.ArgumentTypeError(f"--average must use (w)AR/DR/ARDR")
        return temp
    raise argparse.ArgumentTypeError(f"--average flag is missing <metric> ':' <value> ")

def parse_scale_flag(x):
    x = to_float(x)
    if x <= 0 or x > 400:
        raise ValueError(f"--scale {x} is not valid. x must be (0, 400].")
    return x

def parse_restriction_sites(x):
    sites = parse_flag_brackets(x)

    for site in sites:
        for char in site:
            if char.upper() not in list('ACTGU'):
                raise ValueError('non-nucleotide in restriction site %s'% site)

    rcs = [us.reverse_complement(x,0,False) for x in sites]
    sites.extend(rcs)

    sites = [x.replace('T','U') for x in sites] # convert to RNA
    return sites


