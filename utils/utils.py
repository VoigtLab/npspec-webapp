from PIL import Image
from io import BytesIO
import base64
import hashlib
import streamlit as st
import pandas as pd
import json
import pickle

def hash_smiles(smiles):
    """
    Generate modified md5 hash of input SMILES string.
    """
    return str(int(hashlib.md5(smiles.encode("utf-8")).hexdigest(), 16))

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

