import streamlit as st
from scipy import signal
import pandas as pd
import pickle
import json
from utils.spectranalysis_utils import *
import hashlib

def hash_smiles(smiles):
    """
    Generate modified md5 hash of input SMILES string.
    """
    return str(int(hashlib.md5(smiles.encode("utf-8")).hexdigest(), 16))

@st.cache_data
def load_csv(path_name, transpose=False, *args, **kwargs):
  if transpose:
    return pd.read_csv(path_name, *args, **kwargs).transpose()
  else:
    return pd.read_csv(path_name, *args, **kwargs)

@st.cache_data
def load_json(path_name):
    with open(path_name, 'r') as f:
        return json.load(f)

@st.cache_data
def load_nx_graph(path_name):
    with open(path_name, 'rb') as f:
        g = pickle.load(f)
    return g

def load_building_blocks(path_name, **kwargs):
    building_blocks = load_csv(path_name, **kwargs)
    building_blocks = set([standardize_smiles(str(s)) for s in building_blocks['smiles'].values if '*' not in str(s)])
    return building_blocks

def load_spec_df(spec_df_path, transpose=True):
    spec_df = load_csv(spec_df_path, transpose=transpose, index_col=0)
    return spec_df
    
def get_family_labels(df, smiles_col='smiles'):
    quin = Chem.MolFromSmarts('[#6]:1:[#6](~O):[#6]:[#6]:[#6](~O):[#6]:1')
    carot = Chem.MolFromSmarts('C=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC=C')
    porph = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~1')
    bil = Chem.MolFromSmarts('[#6]~[#6]~[#7;R1]~[#6]~[#6]~[#6]~[#7;R1]~[#6]~[#6]~[#6]~[#7;R1]~[#6]~[#6]~[#6]~[#7;R1]')
    flav = Chem.MolFromSmarts('[#8]~C~1~C~C(~C~2~C~C~C~C~C~2)~[#8]~C3~C~1~C~C~C~C~3'.replace('C','[#6]'))
    cyan =  Chem.MolFromSmarts('C~1~C~C(~C~2~C~C~C~C~C~2)~[#8+]~C3~C~1~C~C~C~C~3'.replace('C','[#6]'))
    phenazine =  Chem.MolFromSmarts('C:1:C:C:C:2:C(:C:1):N:C:3:C:C:C:C:C:3:N2  '.replace('C','[#6]').replace('N','[#7]'))
    
    df['mols'] = df[smiles_col].map(lambda x : Chem.MolFromSmiles(x))
    df['Quinone'] = df['mols'].map(lambda x : x.HasSubstructMatch(quin) if x else False) 
    df['Carotenoid'] = df['mols'].map(lambda x : x.HasSubstructMatch(carot) if x else False)
    df['Porphyrin'] = df['mols'].map(lambda x : x.HasSubstructMatch(porph) if x else False)
    df['Bilin'] = df['mols'].map(lambda x : x.HasSubstructMatch(bil) if x else False)
    df['Flavonoid'] = df['mols'].map(lambda x : x.HasSubstructMatch(flav) if x else False)
    df['Cyanidin'] = df['mols'].map(lambda x : x.HasSubstructMatch(cyan) if x else False)
    df['Phenazine'] = df['mols'].map(lambda x : x.HasSubstructMatch(phenazine) if x else False)
    
    family_ls = []
    families = np.array(['Quinone', 'Carotenoid', 'Porphyrin', 'Bilin', 'Flavonoid', 'Cyanidin', 'Phenazine'])
    for i in df.index:
        family = families[df.loc[i, families].values.astype(bool)]
        if len(family):
            family_ls.append(family[0])
        else:
            family_ls.append("Other")   
    df['Family'] = family_ls
    return df.drop(columns=['mols'])

def get_peaks(spectrum, wl_index, **kwargs):
    peak_idxs = signal.find_peaks(spectrum, **kwargs)
    return wl_index[peak_idxs[0]]
