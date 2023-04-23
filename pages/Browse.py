import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import os
from data.loaded_data import uniqueness_df, pred_spec_bars
from utils.utils import *
from components import display
DISPLAY_STRUCTURES = False


st.title('Browse Predicted Metabolite Spectra Data')


df = uniqueness_df
pred_spec = pred_spec_bars

for hsmi, smi in zip(df['hashed_smiles'].values, df['smiles'].values):
  #check if file already exists
  mol = Chem.MolFromSmiles(smi)
  if not os.path.isfile('imgs/'+hsmi+'.png') and mol is not None:
    img = Draw.MolToImage(mol, size=(600,600), clearBackground=True)
    img.thumbnail((300,300))
    img.save('imgs/'+hsmi+".png")

df_to_display = df.loc[:, ['name', 'smiles', 'lambda_max']]
df_to_display = df_to_display.rename(columns={'smiles': 'SMILES', 'lambda_max': 'Lambda Max'})

if DISPLAY_STRUCTURES:
  df_to_display['structure'] = ['imgs/'+smi+'.png' for smi in df['hashed_smiles'].values]
  html = convert_df(df_to_display)
  st.markdown(
      html,
      unsafe_allow_html=True
  )

else:
  st.dataframe(df_to_display, use_container_width=True)
  display.show_structures_and_spectra(df_to_display, df, pred_spec)
