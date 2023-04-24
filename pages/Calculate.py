import streamlit as st
import numpy as np
import plotly.express as px
from utils.utils import *
from utils.spectranalysis_utils import *
from data.loaded_data import met_graph, ec_building_blocks, uniqueness_df, pred_spec_bars, pred_spec_df
from data.loading_utils import load_building_blocks, load_spec_df
from components import display


st.title("Calculate spectral uniqueness and accessibility")

# Load in the graph

met = met_graph

# load in starting points
# Redefine buyables and step count for E Coli

col1, col2 = st.columns(2)
uploaded_metabolite_file = col1.file_uploader("Upload a file for metabolic starting points", type="csv", 
            help="The file should be a tsv with a column named 'smiles'")
col1.markdown("By default, the starting points are the metabolites naturally produced in E. Coli.")
uploaded_spectra_file = col2.file_uploader("Upload a file for spectral comparison", type="csv", 
            help="The file should be a csv with each row containing an absorbance spectrum and the columns representing wavelengths")
col2.markdown("By default, the spectral comparison is done with the predicted spectra for all metabolites")

# Calculate metabolic distance
if uploaded_metabolite_file is not None:
  building_blocks = load_building_blocks(uploaded_metabolite_file, sep='\t')
else:
  building_blocks = ec_building_blocks
starting_nodes = []
for node in met.nodes:
    if node in building_blocks:
        starting_nodes.append(node)
met = metabolic_dijkstra(met, starting_nodes, path_length_field="path_length", 
                         shortest_path_field="shortest_path")
uniqueness_df['steps'] = [met.nodes[s]['path_length'] if s in met.nodes else np.inf for s in uniqueness_df['smiles'] ]

# Calculate uniqueness
if uploaded_spectra_file is not None:
  uploaded_spec_df = load_spec_df(uploaded_spectra_file, transpose=True)
  background_spectra = filter_df(uploaded_spec_df).values
  spectra_to_score = filter_df(pred_spec_df.transpose())
  scores = calculate_uniquness(background_spectra, spectra_to_score)
  cols_to_keep = [c for c in uniqueness_df.columns if c != 'mean_dist']
  new_uniqueness_df = pd.DataFrame({'name': spectra_to_score.index, 'mean_dist': scores})

  # this is not the right df to merge onto, need to figure that out
  uniqueness_df = uniqueness_df.loc[:, cols_to_keep].merge(new_uniqueness_df, on='name', how='outer')
  uniqueness_df = uniqueness_df.dropna(subset='mean_dist')

# Plot distance vs. uniqueness
to_plot = uniqueness_df.copy()
max_val = int(np.nanmax(uniqueness_df[uniqueness_df['steps']<np.inf]['steps'])+10)
to_plot['steps'] = to_plot['steps'].map(lambda x: max_val if x == np.inf else x)
fig = px.scatter(to_plot, y='steps', x='mean_dist', hover_name='name', 
      labels={'steps': 'Met. Steps', 'mean_dist': 'Uniqueness'}) 
fig.update_layout(
    yaxis = dict(
        tickmode = 'array',
        tickvals = list(range(0, max_val+1, 2)),
        ticktext = [str(x) for x in range(0, max_val, 2)] + ["Unknown"]
    )
)
st.plotly_chart(fig)

# Display additional information
df_to_display = uniqueness_df.loc[:, ['display_name','mean_dist', 'steps', 'smiles']].drop_duplicates(subset=['display_name'])
df_to_display = df_to_display.rename(columns={'display_name':'Name','smiles': 'SMILES', 'mean_dist': 'Uniqueness', 'steps': 'Met. Steps'})
display.show_structures_and_spectra(df_to_display, uniqueness_df, pred_spec_bars)
