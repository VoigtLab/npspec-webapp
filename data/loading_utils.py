import streamlit as st
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
    building_blocks = set([standardize_smiles(str(s)) for s in building_blocks['smiles'].values if '*' not in s])
    return building_blocks

def load_spec_df(spec_df_path, transpose=True):
    spec_df = load_csv(spec_df_path, transpose=transpose, index_col=0)
    return spec_df
    