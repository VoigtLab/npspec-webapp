import streamlit as st  
from data.loaded_data import uniqueness_df
from utils.utils import *

st.title("Search")



st.markdown("<h1 style='font-size:18px;'>Select a color you would like your molecule to absorb using the square below</h1>", unsafe_allow_html=True)
_,col2,_ = st.columns(3)
query_color = col2.color_picker("", value="#FF0000", help="Select a color")


#convert hex to rgb
query_color = tuple(int(query_color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
num_hits = st.slider("Number of hits to display", min_value=1, max_value=50, value=10) 
# st.write("Selected color: {} (R, G, B)".format(query_color))

st.markdown("### Top hits")
st.markdown("More information below")
uniqueness_df['rgb'] = [convert_wl_to_rgb(w) for w in uniqueness_df['lambda_max']]
uniqueness_df['color_dist_to_query'] = [l2_dist(query_color, rgb) for rgb in uniqueness_df['rgb']]
to_display = uniqueness_df.sort_values(by='color_dist_to_query').iloc[:num_hits].loc[:, ['name', 'lambda_max', 'peak_info']]
to_display.rename(columns={'lambda_max': 'Wavelength of max abs. (nm)', 'peak_info':'Peak positions (nm)'}, inplace=True)

col_1, col_2, col_3 = st.columns(3)
for i, idx in enumerate(uniqueness_df.sort_values(by='color_dist_to_query').index[:num_hits]):
  mol = uniqueness_df.loc[idx, 'hashed_smiles']
  name = uniqueness_df.loc[idx, 'name']
  color = '{:02x}{:02x}{:02x}'.format(*[int(x) for x in uniqueness_df.loc[idx, 'rgb']])
  display_text = f'<h1 style="color:#{color};text-align: center;font-size:18px;">{i+1}. {name}:</h1>'
  
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