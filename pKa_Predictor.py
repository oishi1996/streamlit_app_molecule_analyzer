import os
import streamlit as st
from streamlit_ketcher import st_ketcher
from dotenv import load_dotenv
from utils import model_utils
load_dotenv(os.path.join(os.path.dirname(__file__), ".", ".env"))

# モデルのロード
model_path = os.getenv("MODEL_PATH")
if not model_path:
    raise ValueError("MODEL_PATH is not set in the environment variables.")
model = model_utils.load_model(model_path)

# Streamlit UI
st.title(":material/Science: pKa Predictor")
st.markdown("""
            - This app allows you to draw chemical structures and predict their pKa values
            - Draw your molecule below...
            """)


with st.container():
    smiles = st_ketcher()

    if smiles:
        predicted_pka = model_utils.predict_pka(smiles, model)
        st.success(f"**Generated SMILES:** `{smiles}`")
        st.success(f"Predicted pKa: {predicted_pka:.2f}")
    else:
        st.warning("No molecule detected. Please draw a structure.")
