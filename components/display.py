import streamlit as st
import matplotlib.pyplot as plt
import os 

def show_structures_and_spectra(df_to_display, df, pred_spec, name_col='Name'):
  selected_names = st.multiselect('#### Select molecules:', df_to_display[name_col], default=None, max_selections=10, 
  help = 'Enter molecule names to display structure and predicted absorbance peaks')
  selected_indices = df_to_display[df_to_display[name_col].map(lambda x: x in selected_names)].index
  selected_rows = df_to_display.loc[selected_indices, :]
  
  if len(selected_rows)==0:
    st.write('### No molecules selected')

  else:
    st.write('### Select molecules to learn more', selected_rows)
    st.write('Displaying molecule structure and predicted absorbance peaks')
    cmap = plt.get_cmap('Set2')

    col1, col2 = st.columns(2)
    mols_to_show = df.loc[selected_indices, 'hashed_smiles'].values
    smiles_to_show = df.loc[selected_indices, 'smiles'].values
    names_to_show = df_to_display.loc[selected_indices, name_col].values
    for i, (mol, name) in enumerate(zip(mols_to_show, names_to_show)):
      color= '{:02x}{:02x}{:02x}'.format(*[int(x*255) for x in cmap(i)])
      pubchem_link = 'https://pubchem.ncbi.nlm.nih.gov/#query='+smiles_to_show[i]+'&tab=similarity&similaritythreshold=100'
      if os.path.isfile('imgs/'+mol+'.png'):
        # col1.markdown(f'<h1 style="color:#{color};font-size:18px;">{str(i+1)}. [{name}]({pubchem_link}):</h1>', unsafe_allow_html=True)
        # col1.markdown('<h1 style="color:#{};font-size:18px;">{}. [{}]({}):</h1>'.format(color, str(i+1), name, pubchem_link), unsafe_allow_html=True)
        col1.markdown('<a style="color:#{};font-size:18px;">{}. [{}]({}):</a>'.format(color, str(i+1), name, pubchem_link), unsafe_allow_html=True)
        col1.image('imgs/'+mol+'.png', width=200)
      else:
        col1.markdown(f'<a style="color:#{color};font-size:18px;">{str(i+1)}. {name}:</a>', unsafe_allow_html=True)
        col1.markdown(f'<h1 style="color:#{color};font-size:18px;">Cannot display structure</h1>', unsafe_allow_html=True)
    fig, ax = plt.subplots()
    ax.set_facecolor((1,1,1,1))
    plt.rcParams["font.size"]=20
    
    for i, n in enumerate(selected_indices):
      query_name = df.loc[n, 'name']
      if query_name not in pred_spec.keys():
        st.write('Could not retrieve spectrum for ', df_to_display.loc[n, name_col])
      else:
        results = list(zip(*pred_spec[query_name]))
        ax.bar(results[0], results[1], width=3, label = i+1, color = cmap(i), alpha=0.6)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Absorbance intensity (a.u.)')
    ax.legend()
    col2.pyplot(fig)
