import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from queue import PriorityQueue
from matplotlib.patches import ConnectionPatch
import itertools
import re
from scipy.spatial.distance import cdist

def neutralize_atoms(smiles):
    #pulled from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return Chem.MolToSmiles(mol)


def standardize_reaction_smiles(reaction_smiles):
    try:
        reactants, products = [neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(s))) for s in reaction_smiles.split('>>')]
        return '>>'.join([reactants, products])
    except:
        return reaction_smiles
    
    
def standardize_smiles(smiles):
    smiles = smiles.replace("'","")
    if ' ' in smiles or 'R' in smiles or '*' in smiles:
        return smiles
    try:
        smiles = neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        return smiles
    except:
        if Chem.MolFromSmiles(smiles)==None:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
            fix_c = AllChem.ReactionFromSmarts('[#6-:1]>>[C;+0:1]')
            if m:
                f = fix_c.RunReactants([m])
                if len (f)>0:
                    f = f[0][0]
                    print ('fixed smiles :', smiles,Chem.MolToSmiles(f))
                    return neutralize_atoms(Chem.MolToSmiles(f))
                else:
                    return smiles
            else:
                return smiles
        else:
            return smiles


def metabolic_dijkstra (g, starting_nodes, path_length_field="path_length", 
                        shortest_path_field="shortest_pathway"):
    """
    Implementation based on Gallo et al, 1993
    """
    g = g.copy()
    pq = PriorityQueue()  
    for node in g.nodes:
        g.nodes[node][path_length_field] = np.inf
        g.nodes[node][shortest_path_field] = []
        g.nodes[node]['in_queue'] = False
        if '>>' in node:
            g.nodes[node]['visit_counter'] = 0
            g.nodes[node]['pred_count'] = len(list(g.predecessors(node)))

    for node in starting_nodes:
        g.nodes[node][path_length_field] = 0
        if '>>' not in node:
            pq.put((g.nodes[node][path_length_field], node))
            g.nodes[node]['in_queue'] = True

    
    while not pq.empty():
        _, curr_node = pq.get()
        g.nodes[curr_node]['in_queue'] = False
        for reaction_node in list(g.successors(curr_node)):
            g.nodes[reaction_node]['visit_counter'] += 1
            if g.nodes[reaction_node]['visit_counter'] == g.nodes[reaction_node]['pred_count']:
                f = sum([g.nodes[r][path_length_field] for r in g.predecessors(reaction_node)])
                f_path = list(itertools.chain(*[g.nodes[r][shortest_path_field] for r in g.predecessors(reaction_node)]))
                g.nodes[reaction_node][path_length_field] = f
                g.nodes[reaction_node][shortest_path_field] = f_path
                for p in g.successors(reaction_node):
                    orig_path_length = g.nodes[p][path_length_field]
                    if  orig_path_length > f + 1: # currently all edge weights = 1
                        g.nodes[p][path_length_field] = f + 1
                        g.nodes[p][shortest_path_field] = f_path + [reaction_node]
                        if not g.nodes[p]['in_queue']:
                            g.nodes[p]['in_queue'] = True
                            pq.put((g.nodes[p][path_length_field], p))
                            if orig_path_length < np.inf:
                                for e in list(g.successors(p)):
                                    g.nodes[e]['visit_counter'] -= 1

    return g


def flip_reaction (reaction_smiles):
    reactants, products = reaction_smiles.split('>>')
    return '>>'.join([products, reactants])

def clean_name(string: str) -> str:
    new = re.sub("[\s\.,]", "_", string)
    new = re.sub("[\[\]\(\)']", "", new)
    new = re.sub('&rarr;', '->', new)
    new = re.sub('<[a-zA-z]*>|<\/[a-zA-z]*>|;|&|^[Aa]n |^[0-9]* ','', new)
    new = re.sub('\+', 'plus', new)
    new = re.sub('^-', 'minus', new)
    new = re.sub(',|\/|\\\\', '-', new)
    return new

def remove_consecutive_ints (int_ls, method='keep_first'):
    """
    Args:
        int_ls (list[int]): list of ascending integer values
        method (str): `keep_first` or `midpoint`. If `keep_firt`
            keep the lowest value integer in the list of consecutive
            values. If `midpoint`, keep the midpoint value in the list 
            of consecutive integers.
    Returns:
        a list of integers without consecutive values
    """
    res_ls = []
    if method=='keep_first':
        last_number = np.inf
        for i in int_ls:
            if (i - last_number) == 1:
                pass
            else:
                res_ls.append(i)

            last_number = i
    
    elif method=='midpoint':
        last_number = np.inf
        consecutive_lists = []
        for i in int_ls:
            if (i - last_number) == 1:
                consecutive_lists[-1].append(i)
            else:
                consecutive_lists.append([i])
            last_number = i
        res_ls = [ls[int(len(ls)/2)] for ls in consecutive_lists] # find midpoint
    else:
        res_ls = None
        print ('No such method: {}!'.format(method))
    return res_ls

def parse_spectrum (spectrum, intensity_threshold=0.05, derivative_threshold=0.0001, method='midpoint',
                   threshold_subspectra = True, validate=False):
    """
    Split a spectrum at valley points into a list of spectra 
    Args:
        spectrum (array-like): spectrum to parse
        intensity_threshold (float): relative intensity value under which a value in the spectrum
            is considered a valley (separating peaks)
        derivative_threshold (float): value for the derivative of a point with a positive second 
            derivative under which the point is considered to be a valley
        method (string): method to use for picking which point in a valley to split the spectra on
        threshold_subspectra (bool): if true, remove subspectra that do not have any points above the 
            `intensity_threshold`
        validate(bool): if true, raise an error if there is a discrepancy greater than 5% between
            the sum of the area under the subspectra and the area under the original spectrum
    Returns:
        A list of subspectra that sum to `spectrum`. Each subspectrum contains a single 
        identified peak
    """
    
    spectrum = np.array(spectrum)
    norm_spectrum = spectrum / np.max(spectrum)
    ERROR_MARGIN = 0.05*np.sum(norm_spectrum)

    dx = np.gradient(norm_spectrum)
    ddx = np.gradient(dx)
    
    min_idxs = np.where((norm_spectrum<intensity_threshold) | (np.abs(dx) <= derivative_threshold) * (ddx > 0))[0]
    split_points = remove_consecutive_ints(min_idxs, method)
    if len(split_points):
        split_points = [0] + split_points + [-1]
        basis_spectra = []
        
        for i in range(len(split_points)-1):
            partial_spec = spectrum.copy()
            partial_spec[:split_points[i]] = 0
            partial_spec[split_points[i+1]:] = 0
            
            # don't include fake spectra 
            if threshold_subspectra and any(norm_spectrum[split_points[i]:split_points[i+1]]>intensity_threshold):
                basis_spectra.append(partial_spec)
            elif not threshold_subspectra:
                basis_spectra.append(partial_spec)
                
    else:
        return [spectrum]
#     print (np.sum(np.abs(np.sum(basis_spectra, axis=0)/np.max(spectrum) - norm_spectrum)))
    if validate:
        assert 1-np.sum(np.abs(np.sum(basis_spectra, axis=0)/np.sum(spectrum))) < ERROR_MARGIN
    return basis_spectra
    
    
def filter_df (df, lower_bound=350, upper_bound=1000, int_threshold=0.005):
    """
    Args:
    df: dataframe to filter
    lower_bound: lower bound of the wavelength range to keep
    upper_bound: upper bound of the wavelength range to keep
    int_threshold: minimum intensity threhold. Spectra with a sum of values below this value are filtered out
    """
    df = df.copy()
    print ("DOCLSN",df.columns)
    df.columns = df.columns.astype(float)
    c_idxs = [c for c in df.columns if c > lower_bound and c < upper_bound]
    df = df.loc[:, c_idxs]
    df = df[(df.sum(axis=1))>int_threshold]
    return df


def parse_spectra(df):
    original_spectra = df.values
    # parse the spectra by peak -- can be sensitive to parsing parameters
    failed_to_parse = []
    INTENSITY_THRESHOLD=0.05
    DERIVATIVE_THRESHOLD=0.0005
    parsed_spectra = []
    for i in range(original_spectra.shape[0]):
        parsed_spectra.append(parse_spectrum(original_spectra[i,:],INTENSITY_THRESHOLD, DERIVATIVE_THRESHOLD,'midpoint', validate=True))
        
    assert len(failed_to_parse) == 0
            
    counter = 0
    #track which spectra correspond to which subspectra
    original_to_parsed_idxs_dict = {}
    for i in range(len(original_spectra)):
        parsed_idxs = []
        for j in parsed_spectra[i]:
            parsed_idxs.append(counter)
            counter += 1
        original_to_parsed_idxs_dict[i] = parsed_idxs
    flat_parsed_spectra = [x for y in parsed_spectra for x in y]
    return flat_parsed_spectra, original_to_parsed_idxs_dict


def calculate_uniquness (background_specs, spectra_to_score, metric="cosine"):
    parsed_specs, original_to_parsed_idxs_dict = parse_spectra(spectra_to_score)
    distances = cdist(background_specs, parsed_specs, metric=metric)
    scores = []
    for i in range(len(original_to_parsed_idxs_dict.keys())):
        idxs = original_to_parsed_idxs_dict[i]
        score = np.sum(np.nanmean(distances[:,idxs], axis=0))
        scores.append(score)
    return scores
