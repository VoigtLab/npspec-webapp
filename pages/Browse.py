import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import os
from data.loaded_data import uniqueness_df, pred_spec_bars
from utils.utils import *
from components import display
import re
DISPLAY_STRUCTURES = False


st.title('Browse Predicted Metabolite Spectra Data')

df = uniqueness_df
pred_spec = pred_spec_bars

for hsmi, smi in zip(df['hashed_smiles'].values, df['smiles'].values):
  #check if file already exists
  mol = Chem.MolFromSmiles(re.sub('R[0-9]*', '*', smi))
  if not os.path.isfile('imgs/'+hsmi+'.png') and mol is not None:
    img = Draw.MolToImage(mol, size=(600,600), clearBackground=True)
    img.thumbnail((300,300))
    img.save('imgs/'+hsmi+".png")
st.markdown('\n')
wls_to_show = st.slider(
    r'Select a range of $\lambda_{max}$ (nm) values to filter by:',
    0.0, 2000.0, (0.0, 1000.0))
wl_filt = (df['lambda_max']>wls_to_show[0]) & (df['lambda_max']<wls_to_show[1])
st.markdown('`{}` metabolites predicted to have peak of maximal absorbance in the desired range'.format(wl_filt.sum()))
st.markdown('\n')
df_to_display = df.loc[:, ['display_name', 'lambda_max', 'smiles']].drop_duplicates(subset=['display_name'])
df_to_display = df_to_display.rename(columns={'display_name':'Metabolite Name', 'smiles': 'SMILES', 'lambda_max': 'Lambda Max'})
# df_to_display['name'] = convert_names_for_display(df_to_display)

if DISPLAY_STRUCTURES:
  df_to_display['structure'] = ['imgs/'+smi+'.png' for smi in df['hashed_smiles'].values]
  html = convert_df(df_to_display[wl_filt])
  st.markdown(
      html,
      unsafe_allow_html=True
  )

else:
  container = st.container()
  display_smiles = st.checkbox("Display molecule SMILES strings", value=False)
  if display_smiles:
    container.dataframe(df_to_display[wl_filt], use_container_width=True)
  else:
    container.dataframe(df_to_display[wl_filt].drop(columns=['SMILES']), use_container_width=True)
  display.show_structures_and_spectra(df_to_display, df, pred_spec, name_col='Metabolite Name')
