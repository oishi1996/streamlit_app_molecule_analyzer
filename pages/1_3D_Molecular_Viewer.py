import streamlit as st
from stmol import showmol
from streamlit_ketcher import st_ketcher
from utils import molecule_handler

# Streamlit UI
st.title(":material/Science: 3D Molecular Viewer")
st.markdown("""
            - This app allows you to draw chemical structures and generate 3D structure
            - Draw your molecule below...
            """)

# --- SMILES 入力方法の選択 ---
input_mode = st.radio("Choose input method:", ["Text Input", "Draw Molecule"], horizontal=True)

if input_mode == "Text Input":
    smiles = st.text_input("Enter a SMILES string:", 
        "COc3nc(OCc2ccc(C#N)c(c1ccc(C(=O)O)cc1)c2P(=O)(O)O)ccc3C[NH2+]")
else:
    smiles = st_ketcher()


# 表示スタイル選択
view_style = st.selectbox("Select a view style:", ["stick", "line", "sphere"])

if smiles:
    mol, view = molecule_handler.get_optimized_3Dmol(smiles, view_style)
    if mol:
        st.success("3D structure generated successfully.")
        showmol(view, height=500, width=800)
    else:
        st.error("3D structure generation failed. Please check the SMILES string.")

