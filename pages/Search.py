import streamlit as st  
from data.loaded_data import uniqueness_df, name_converter
from utils.utils import *

st.title("Search")



st.markdown("<h1 style='font-size:18px;'>Search molecules by predicted absorbance color using the selector tool in the sidebar</h1>", unsafe_allow_html=True)
st.sidebar.markdown("<h1 style='font-size:20px;'>Color selector:</h1>", unsafe_allow_html=True)
_, _, col, _, _ = st.sidebar.columns(5)
query_color = col.color_picker("Color selection", value="#00FF00", help="Select a color", label_visibility="collapsed")

#convert hex to rgb
query_color = tuple(int(query_color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
st.sidebar.markdown("<h1 style='font-size:20px;'>Number of hits to display:</h1>", unsafe_allow_html=True)
num_hits = st.sidebar.slider("Number of hits to display", min_value=1, max_value=50, value=10, label_visibility="collapsed") 

st.markdown("### Top hits")
st.markdown("More information below")
uniqueness_df['rgb'] = [convert_wl_to_rgb(w) for w in uniqueness_df['lambda_max']]
uniqueness_df['color_dist_to_query'] = [l2_dist(query_color, rgb) for rgb in uniqueness_df['rgb']]
to_display = uniqueness_df.sort_values(by='color_dist_to_query').iloc[:num_hits].loc[:, ['display_name', 'lambda_max', 'peak_info']]
to_display.rename(columns={'display_name':'Name', 'lambda_max': 'Wavelength of max abs. (nm)', 'peak_info':'Peak positions (nm)'}, inplace=True)

col_1, col_2, col_3 = st.columns(3)
for i, idx in enumerate(uniqueness_df.sort_values(by='color_dist_to_query').index[:num_hits]):
  mol = uniqueness_df.loc[idx, 'hashed_smiles']
  name = uniqueness_df.loc[idx, 'name']
  if name in name_converter.keys():
    name = name_converter[name]
  if len(name) > 25:
    name = name[:22] + '... '
  color = '{:02x}{:02x}{:02x}'.format(*[int(x) for x in uniqueness_df.loc[idx, 'rgb']])
  display_text = f'<h1 style="color:#{color};text-align: left;font-size:14px;">{i+1}. {name}:</h1>'
  
  # get the reaminder
  if i%3==0:
    col = col_1
  elif i%3==1:
    col = col_2
  else:
    col = col_3
  col.markdown(display_text, unsafe_allow_html=True)
  col.image('imgs/'+mol+'.png', width=200)
  
  

st.dataframe(to_display, use_container_width=True)