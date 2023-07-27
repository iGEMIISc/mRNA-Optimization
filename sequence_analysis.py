from ribotree_parser import *
from evaluate import ligand_score_func, ligand_polynomial
import itertools
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
import sympy


def add_extra_condition(args_d, conc_list):
    # now add all conditions for plotting/analysis
    condition_list = args_d['condition_list']
    conc_list_og = []
    for condition in condition_list:
        cond = args_d[condition]
        conc_list_og.append(cond[0])
        if len(cond) == 3: 
            conc_list_og.append(cond[1])
    
    condition_list = []
    
    for idx, conc in enumerate(conc_list):
        key = f'condition_{idx+1}'
        condition_list.append(key)
        args_d[key] = [conc, [None]]
    args_d['condition_list'] = condition_list
    return args_d



def signal_prediction(seq, conc_list, filename, package=None, T=37, args_d=None):
    """ Function for predicting sequence response 
    
    args:
        seq: sequence of riboswitch/dna switch
        conc_list: numpy array of all the conditions to be tested (from get_conc_list)
        filename: filename of input file (file used to run RiboTree)
        package: package to use (vienna or contrafold)
        T: temperature (C)
        args_d: dictionary detailing all the aptamers, etc.

    return:
        signal_list: signal response for each concentration
        binary_list: binary response for each concentration (lowest energy microstate Z)
        kd_list: kd response for each concentration
    """

    if args_d is None:
        # Get the ards_d dictionary detailing all the aptamers
        # Also adds the concentration conditions
        args_d = get_args(conc_list, filename)

    if package is not None:
        args_d['package'] = package.lower()
    args_d['Temp'] = T

    # Calculate the signal at all concentrations
    signal_list, binary_list, kd_list = ligand_score_func(seq, args_d, verbose=False)
    signal_list = np.array(signal_list)
    binary_list = np.array(binary_list) # only looks at the lowest energy structure at each condition
    kd_list = np.array(kd_list)
    return signal_list, binary_list, kd_list

def get_NZk(seq, filename, package=None, T=37, args_d=None):
    """ Returns N, Z, kd list for a given sequence and file """
    if args_d is None:
        # Get the ards_d dictionary detailing all the aptamers
        args_d = process_file(filename)

    if package is not None:
        args_d['package'] = package.lower()
    args_d['Temp'] = T

    Z_arr, N_arr, kd_list = ligand_polynomial(seq, args_d, verbose=False)
    return Z_arr, N_arr, kd_list

def calc_polynomial(N_arr, Z_arr, symbol_list=None, output_idx=-1):
    """ Returns the numerator and denominator for the polynial given N and Z """
    num_symbol = N_arr.shape[1]

    # organizes the symbols
    if symbol_list is None:
        symbol_list = ' '.join([f'x_{idx}' for idx in range(1, num_symbol+1)])
    else:
        if num_symbol != len(symbol_list.split(' ')):
            raise ValueError('Too many or too few symbols')

    # Convert symbol to variables
    var_list = sympy.symbols(symbol_list)

    # Get all terms
    term_list = []
    for N, Z in zip(N_arr, Z_arr):
        temp = 1
        for idx in range(len(N)):
            temp *= var_list[idx]**int(N[idx])
        term_list.append(Z*temp)

    # Get ratio expression
    num = sum([term for term in term_list if var_list[output_idx]  in term.free_symbols])
    den = sum(term_list)
    return num, den, var_list

def get_polynomial(seq, filename, package=None, T=37, args_d=None, symbol_list=None, output_idx=-1):
    Z_arr, N_arr, kd_list = get_NZk(seq, filename, package=package, T=T, args_d=args_d)
    num, den, var_list = calc_polynomial(N_arr, Z_arr, symbol_list=symbol_list, output_idx=output_idx)
    return num, den, var_list

def get_conc_list(num_input, lower_bound, upper_bound, output_idx=-1, output_conc=1, num_point=50, log=True):
    """ Returns a list of concentrations from lower_bound to upper_bound 

    args:
        num_input: number of inputs to the riboswitch (do not count the output aptamer)
        lower_bound: lower_bound concentration [M] (log space)
            e.g. -9 is 1e-9
        upper_bound: upper_bound concentration [M] (log space)
        output_idx: specifies which ligand is the output ligand
            Default is the last ligand
        output_conc: concentration of the output ligand (Default = 1)
        num_point: number of points to generate between the bounds (Default = 50)

    return:
        conc_list: numpy array of concentration where each row is a different condition
        x: numpy array of the concentrations iterated over (1-D array)
    """
    if not log:
        lower_bound = np.log10(lower_bound)
        upper_bound = np.log10(upper_bound)
    conc_list = []
    x = np.logspace(lower_bound, upper_bound, num_point)
    conc_list = np.array(list(itertools.product(x, repeat=num_input)))
    if output_idx == -1:
        output_idx = num_input
    conc_list = np.insert(conc_list, output_idx, output_conc*np.ones(conc_list.shape[0]), axis=1)
    return conc_list, x

def plot_signal_2D(x, y, z, label=[], title=''):
    """ Contour plot """

    X, Y = np.meshgrid(x,y)
    Z = z.reshape(X.shape)
    
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes()
    ax.set_aspect('equal')
    n_level = 15
    plt.contourf(X, Y, Z, n_level, cmap='viridis')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar()
    if len(label) == 2:
        plt.xlabel(label[0])
        plt.ylabel(label[1])
    else:
        plt.xlabel('Input Ligand 1 [M]')
        plt.ylabel('Input Ligand 2 [M]')
    AR = np.amax(z)/np.amin(z)
    if len(title): title = title + ' : '
    plt.title(f'{title}AR {AR:.2f}')
    plt.show()
    return AR

def plot_signal_1D(x, y, c='k', scatter=False, alpha=1, s=20, vline=None, title=''):
    """ 1D plot 

    """

    # sort
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes()
    ax.set_xscale('log')
    if scatter:
        if c is not None:
            plt.scatter(x, y, c=c, cmap='viridis', alpha=alpha, s=s)
            plt.colorbar()
        else:
            plt.scatter(x, y, c=c, alpha=alpha, s=s)
    else:
        plt.plot(x,y)
    if vline is not None:
        plt.axvline(x=vline, c='k')
    plt.ylim([0,np.ceil(np.amax(y))])
    plt.xlabel('Input Ligand [M]')
    plt.ylabel('Output Bound [%]')
    AR = np.amax(y)/np.amin(y)
    if len(title): title = title + ' : '
    plt.title(f'{title}AR {AR:.2f}')
    plt.show()
    return AR

def get_ROC(signal, calc_t, cut_off):
    # Compute ROC
    y = (calc_t > cut_off).astype(int)
    if len(signal.shape) == 1:
        X = signal[:,None]
    else:
        X = signal

    # Calculate population difference
    # Delta from Above_Thresh - Below_Thresh
    dp = np.mean(X[y==1]) - np.mean(X[y==0])

    # Fit the classifier
    clf = linear_model.LogisticRegression()
    clf.fit(X, y)

    # Compute ROC curve and ROC area
    y_score = clf.decision_function(X)
    y_hat = clf.predict(X)
    fpr, tpr, _ = roc_curve(y, y_score)
    y_cutoff = -clf.intercept_/clf.coef_
    roc_auc = auc(fpr, tpr)
    cm  = confusion_matrix(y, y_hat)
    return roc_auc, cm, dp