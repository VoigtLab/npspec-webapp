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
uniqueness_df = uniqueness_df.dropna(subset='ws_dist')
met = met_graph
# load in starting points
# Redefine buyables and step count for E Coli

col1, col2 = st.columns(2)
col1.markdown('### 1. Chassis selection')
col2.markdown('### 2. Spectral background selection')
uploaded_metabolite_file = col1.file_uploader("(Optional) Upload a file for metabolic starting points or select an organism below", type="csv", 
            help="The file should be a tsv with a column named 'smiles'")
col2.markdown("By default, molecules are organized by thieir 'global' uniqueness score")

metric_map = metrics = {c.replace('_', ' '):c for c in uniqueness_df}
metrics = [c.replace('_', ' ') for c in uniqueness_df \
              if pd.api.types.is_numeric_dtype(uniqueness_df[c])\
              and not pd.api.types.is_bool_dtype(uniqueness_df[c])\
              and all(uniqueness_df[c]<np.inf) and c != 'color_dist_to_query']
uniqueness_metric = col2.selectbox(
        "## Select a metric ",
        metrics
    )
selected_metric = metric_map[uniqueness_metric]

drop_down_disabled=False
if uploaded_metabolite_file is not None:
  drop_down_disabled = True
organism = col1.selectbox(
        "## Select a chassis organism or upload a custom file above",
        ("E. coli", "B. subtilis", "A. thaliana", "R. gelatinosus", "R. capsulatus", "P. putida"),
        disabled=drop_down_disabled
    )
organsim_dict = {"E. coli":"EC",
                 "B. subtilis":"BS", 
                 "A. thaliana":"AT", 
                 "R. gelatinosus":"RG", 
                 "R. capsulatus":"RC",
                 "P. putida":"PP"}

# Calculate metabolic distance

if uploaded_metabolite_file is not None:
  drop_down_disabled=True
  building_blocks = load_building_blocks(uploaded_metabolite_file, sep='\t')
  starting_nodes = []
  for node in met.nodes:
      if node in building_blocks:
          starting_nodes.append(node)
  met = metabolic_dijkstra(met, starting_nodes, path_length_field="path_length", 
                          shortest_path_field="shortest_path")
  uniqueness_df['steps'] = [met.nodes[s]['path_length'] if s in met.nodes else np.inf for s in uniqueness_df['smiles'] ]
else:
  uniqueness_df['steps'] = uniqueness_df["{}_steps".format(organsim_dict[organism])]


# Plot distance vs. uniqueness
uniqueness_to_show = st.slider(
    r'Select a range for uniqueness {} to filter by:'.format(uniqueness_metric),
    0.0, float(int(np.max(uniqueness_df[selected_metric]))+1), 
    (0.0, float(int(np.max(uniqueness_df[selected_metric]))+1)))
u_filt = (uniqueness_df[selected_metric]>=uniqueness_to_show[0]) & (uniqueness_df[selected_metric]<=uniqueness_to_show[1])

filt_steps = st.checkbox('Filter by number of steps')
s_filt = uniqueness_df['steps']>=0
if filt_steps:
  non_inf = uniqueness_df[uniqueness_df['steps'] < np.inf]
  steps_to_show = st.slider(
      r'Select a range for number of heterelogous enzymatic steps to filter by:',
      0, int(np.max(non_inf['steps']))+1, 
      (0, int(np.max(non_inf['steps']))+1))
  s_filt = (uniqueness_df['steps']>=steps_to_show[0]) & \
            (uniqueness_df['steps']<=steps_to_show[1]) 

to_plot = uniqueness_df.copy()[u_filt & s_filt]
max_val = int(np.nanmax(uniqueness_df[uniqueness_df['steps']<np.inf]['steps'])+10)
to_plot['steps'] = to_plot['steps'].map(lambda x: max_val if x == np.inf else x)
fig = px.scatter(to_plot, y='steps', x=selected_metric, hover_name='display_name', color='Family', 
      labels={'steps': 'Met. Steps', selected_metric: 'Uniqueness ({})'.format(uniqueness_metric)}) 
fig.update_layout(
    yaxis = dict(
        tickmode = 'array',
        tickvals = list(range(0, max_val+1, 2)),
        ticktext = [str(x) for x in range(0, max_val, 2)] + ["Unknown"]
    )
)
st.plotly_chart(fig)

# Display additional information
df_to_display = uniqueness_df.loc[:, ['display_name',selected_metric, 'steps', 'smiles']].drop_duplicates(subset=['display_name'])
df_to_display = df_to_display.rename(columns={'display_name':'Name','smiles': 'SMILES', selected_metric: 'Uniqueness ({})'.format(selected_metric), 'steps': 'Met. Steps'})

st.markdown("`{}` metabolites match filters".format(len(df_to_display[u_filt & s_filt])))
st.dataframe(df_to_display[u_filt & s_filt].sort_values('Uniqueness ({})'.format(selected_metric), ascending=False), use_container_width=True)

display.show_structures_and_spectra(df_to_display, uniqueness_df, pred_spec_bars)
