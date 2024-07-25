import streamlit as st  
from data.loaded_data import uniqueness_df, reverse_name_converter, reaction_df, met_graph
from utils.utils import *
import re 

st.title("Explore the biosynthetic routes of molecules")

#select organism
st.markdown('### Chassis selection')
organism = st.selectbox(
        "## Select a chassis organism to define available starting metabolites",
        ("E. coli", "B. subtilis", "A. thaliana", "R. gelatinosus", "P. putida"),
    )
organsim_dict = {"E. coli":"EC",
                 "B. subtilis":"BS", 
                 "A. thaliana":"AT", 
                 "R. gelatinosus":"RG",
                 "P. putida":"PP"}

uniqueness_df['steps'] = uniqueness_df["{}_steps".format(organsim_dict[organism])]

disp_names = []
for idx in uniqueness_df.index:
  disp_names.append(str(uniqueness_df.loc[idx, 'steps'])  + ' step(s) - ' + str(uniqueness_df.loc[idx, 'display_name']).replace('inf', 'Unknown'))


df_to_display = uniqueness_df.loc[:, ['display_name', 'steps', 'smiles']]
df_to_display['display_name'] = disp_names
df_to_display = df_to_display.drop_duplicates(subset=['display_name'])

selected_name = st.selectbox('#### Select target molecule:', df_to_display['display_name'])


selected_name = ' - '.join(selected_name.split(' - ')[1:])
if selected_name in reverse_name_converter.keys():
  selected_name = reverse_name_converter[selected_name]

target_entry = uniqueness_df[uniqueness_df['name']==selected_name]
if len(target_entry):
  idx = target_entry.index

  smiles = target_entry.loc[idx, 'smiles'].values[0]
  col1, col2, _= st.columns(3)

  pubchem_link = 'https://pubchem.ncbi.nlm.nih.gov/#query='+smiles+'&tab=similarity&similaritythreshold=100'
  
  col1.markdown('<a style="font-size:18px;">[Target structure : ]({})</a>'.format(pubchem_link), unsafe_allow_html=True)
  
  fig, ax = plt.subplots(1)
  ax.imshow(Chem.Draw.MolToImage(Chem.MolFromSmiles(smiles)))
  plt.xticks([])
  plt.yticks([])
  plt.box(False)
  col2.pyplot(fig, use_container_width = True)

  rxn_list = met_graph.nodes[smiles]['shortest_path_{}'.format(organsim_dict[organism])]
  
  reaction_str_ls = []
  for step_num, rxn in enumerate(rxn_list):
    rd_obj = Chem.ReactionFromSmarts(rxn, useSmiles=True)
    entry = reaction_df[reaction_df['smiles'] == rxn].copy()
    if len(entry):
      entry['Reaction'] = entry['reaction_str'].fillna(entry['reaction_string']).fillna(entry['Reaction'])
      rhea_link_template = 'https://www.rhea-db.org/rhea/'
      kegg_link_template = 'https://www.genome.jp/dbget-bin/www_bget?rn:'
      metacyc_link_template = 'https://metacyc.org/reaction?id='
      
      rhea_link = entry['ID_rhea'].map(lambda x : rhea_link_template+re.findall('(?<=RHEA:)[0-9]+', x)[0] if not pd.isna(x) else None).values[0]
      metacyc_link = entry['UNIQUE-ID'].map(lambda x : metacyc_link_template+x if not pd.isna(x) else None).values[0]
      brenda_link = entry['Reaction_ID_BRENDA'].map(lambda x : get_brenda_link(x) if not pd.isna(x) else None).values[0]
      kegg_link = entry['Reaction_ID_KEGG'].map(lambda x : kegg_link_template+str(x).split(',')[0] if not pd.isna(x) else None).values[0]

    links_ls = []
    for db_name, link in zip(['Rhea', 'MetaCyc', 'BRENDA', 'KEGG'], 
                             [rhea_link, metacyc_link, brenda_link, kegg_link]):
      if link is not None:
        links_ls.append('<a>[{}]({})</a>'.format(db_name, link))
    
    links = ', '.join(links_ls)

    st.markdown(str(step_num+1) + '. ' + entry['Reaction'].values[0]+' ['+links+']', unsafe_allow_html=True)

    fig, ax = plt.subplots(1)
    ax.imshow(Chem.Draw.ReactionToImage(rd_obj, subImgSize=(300,500)))
    plt.xticks([])
    plt.yticks([])
    plt.box(False)
    st.pyplot(fig)
