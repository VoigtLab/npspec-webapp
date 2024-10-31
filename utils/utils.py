from PIL import Image
from io import BytesIO
import base64
import hashlib
import streamlit as st
import pandas as pd
import json
import pickle
import numpy as np
from data.loaded_data import name_converter
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt

def convert_wl_to_rgb(w):
    if (w >= 380 and w < 440):
        red   = -(w - 440) / (440 - 380)
        green = 0.0
        blue  = 1.0
    elif(w >= 440 and w < 490):
        red   = 0.0
        green = (w - 440) / (490 - 440)
        blue  = 1.0
    elif (w >= 490 and w < 510):
        red   = 0.0
        green = 1.0
        blue  = -(w - 510) / (510 - 490)
    elif (w >= 510 and w < 580):
        red   = (w - 510) / (580 - 510)
        green = 1.0
        blue  = 0.0
    elif (w >= 580 and w < 645):
        red   = 1.0
        green = -(w - 645) / (645 - 580)
        blue  = 0.0
    elif (w >= 645 and w < 781):
        red   = 1.0
        green = 0.0
        blue  = 0.0
    elif w <380:
        red = 1
        green = 0
        blue = 1
    elif w > 781:
        red = 1
        green = 0
        blue = 0
    else:
        red   = 0.0
        green = 0.0
        blue  = 0.0
    #Let the intensity fall off near the vision limits
    if (w >= 200 and w < 420):
        factor = 0.3 + 0.7*(w - 200) / (420 - 200)
    elif (w >= 420 and w < 701):
        factor = 1.0
    elif (w >= 701 and w < 1000):
        factor = 0.3 + 0.7*(1000 - w) / (1000 - 700)
    else:
        factor = 0.0
    
    gamma = 0.80
    R = 255*(red * factor)** gamma if red > 0 else 0
    G = 255*(green * factor)** gamma if green > 0 else 0
    B = 255*(blue * factor)** gamma if blue > 0 else 0
    return R, G, B

def l2_dist(arr1, arr2):
    return np.sqrt(np.sum(np.square([x1-x2 for x1,x2 in zip(arr1, arr2)])))

def get_thumbnail(path: str) -> Image:
    img = Image.open(path)
    img.thumbnail((200, 200))
    return img

def image_to_base64(img_path: str) -> str:
    try:
      img = get_thumbnail(img_path)
      with BytesIO() as buffer:
          img.save(buffer, 'png') # or 'jpeg'
          return base64.b64encode(buffer.getvalue()).decode()
    except:
      return None

def convert_df(input_df):
     return input_df.to_html(escape=False, formatters=dict(structure=image_formatter))

def image_formatter(img_path: str) -> str:
    return f'<img src="data:image/png;base64,{image_to_base64(img_path)}">'

def convert_names_for_display(df):
    return [name_converter[n] if n in name_converter.keys() else n for n in df['name'].values]

def get_brenda_link(brenda_id):
    # some code is copied from clean_new_brenda_smiles.ipynb
    # should have just replaced the id in the dataframe maybe
    if pd.isna(brenda_id):
        return None
    else:
        brids = brenda_id.split(',')
        if len(brids) == 0:
            return none
        elif len(brids) == 1:
            selected_brid = brids[0]
    
        else:
            #check if there are BR and BS
            if not all([x.startswith('BR') for x in brids]):
                brids = [x for x in brids if not x.startswith('BR')]
            selected_brid = brids[0]

        if selected_brid.startswith('BR'):
            return 'https://www.brenda-enzymes.org/structure.php?show=reaction&id={}&type=I&displayType=marvin'.format(selected_brid.replace('BR',''))
        elif selected_brid.startswith('BS'):
            return 'https://www.brenda-enzymes.org/structure.php?show=reaction&id={}&type=S&displayType=marvin'.format(selected_brid.replace('BS',''))
        else:
            print ('WARNING: Unknown Brenda ID type')
            return None
