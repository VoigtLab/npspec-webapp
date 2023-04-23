import streamlit as st
import pandas as pd
import pickle
import json

@st.cache_data
def load_csv(path_name, *args, **kwargs):
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
