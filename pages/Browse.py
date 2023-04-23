import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import os
from data.loading_utils import *
from utils.utils import *
import matplotlib.pyplot as plt
import plotly.express as px
DISPLAY_STRUCTURES = False


st.title('Browse Predicted Metabolite Spectra Data')


df = load_csv('../spectranalysis/22Apr2023_computed_uniqueness_scores.csv')
df['hashed_smiles'] = [str(hash_smiles(s)) for s in df['smiles']]

pred_spec = load_json('../biospectral/data/spectra/all_accessible_metabolites/7Apr2023_all_accesible_bars.json')
for k in list(pred_spec.keys()):
  pred_spec[k.replace('_b3lyp','').split('/')[-1]] = pred_spec[k]

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
  selected_names = st.multiselect('#### Select molecules:', df_to_display['name'], default=None, max_selections=10, 
  help = 'Enter molecule names to display structure and predicted absorbance peaks')
  # st.write('Enter molecule names above to display structure and predicted absorbance peaks', use_container_width=True)
  selected_indices = df_to_display[df_to_display['name'].map(lambda x: x in selected_names)].index
  selected_rows = df_to_display.loc[selected_indices, :]
  
  if len(selected_rows)==0:
    st.write('### No molecules selected', use_container_width=True)

  else:
    st.write('### Selected molecules', selected_rows, use_container_width=True)
    st.write('Displaying molecule structure and predicted absorbance peaks', use_container_width=True)
    cmap = plt.get_cmap('Set2')

    col1, col2 = st.columns(2)
    mols_to_show = df.loc[selected_indices, 'hashed_smiles'].values
    names_to_show = df.loc[selected_indices, 'name'].values
    for i, (mol, name) in enumerate(zip(mols_to_show, names_to_show)):
      if os.path.isfile('imgs/'+mol+'.png'):
        color= '{:02x}{:02x}{:02x}'.format(*[int(x*255) for x in cmap(i)])
        col1.markdown(f'<h1 style="color:#'+color+';font-size:18px;">'+str(i+1)+'. '+name+':</h1>', unsafe_allow_html=True)
        col1.image('imgs/'+mol+'.png', width=200)
      else:
        col1.write('Could not retrieve structure for {}'.format(name))
    
    fig, ax = plt.subplots()
    ax.set_facecolor((1,1,1,1))
    plt.rcParams["font.size"]=20
    
    for i, n in enumerate(selected_names):
      if n not in pred_spec.keys():
        st.write('Could not retrieve spectrum for ', n)
      else:
        results = list(zip(*pred_spec[n]))
        ax.bar(results[0], results[1], width=3, label = i+1, color = cmap(i), alpha=0.6)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Absorbance intensity (a.u.)')
    ax.legend()
    col2.pyplot(fig)