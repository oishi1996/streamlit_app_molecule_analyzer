import os
import streamlit as st
from streamlit_ketcher import st_ketcher

from utils import model_utils
from services.initialize_setting import initialize_setting
from database.database import get_db_session
from database.models import PkaData

initialize_setting()

# モデルのロード
model_path = os.getenv("MODEL_PATH")
if not model_path:
    raise ValueError("MODEL_PATH is not set in the environment variables.")
model = model_utils.load_model(model_path)

# Streamlit UI
st.title(":material/Science: pKa Predictor")
st.markdown("""
            - This app allows you to draw chemical structures and predict their pKa values
            - (1). **Draw** molecule below
            - (2). **Click** the **Apply** button to predict pKa value
            """)


with st.container():
    smiles = st_ketcher()

    if smiles:
        with st.spinner("Predicting pKa..."):
            predicted_pka = model_utils.predict_pka(smiles, model)

        if predicted_pka:
            st.success(f"**Generated SMILES:** {smiles}  \n" + f"**Predicted pKa:** {predicted_pka:.2f}")

            #データベースに保存
            with get_db_session() as db:
                new_pkadata = PkaData(smiles=smiles, predicted_pka=predicted_pka)
                db.add(new_pkadata)
                db.commit()
                db.refresh(new_pkadata)
                st.toast("successfuly saved pKa", icon="🎉")

        else:
            st.error("Invalid SMILES. Please try again.")
    
    else:
        st.warning("No molecule detected. Please draw a structure.")


# データベース取得
with get_db_session() as db:
    recent_entries = db.query(PkaData).order_by(PkaData.created_at.desc()).limit(5).all()

    if recent_entries:
        st.subheader(":material/Database: Recent pKa data (top 5)")
        
        # 表形式に変換
        data = [
            {
                "SMILES": entry.smiles,
                "Predicted pKa": round(entry.predicted_pka, 2),
                "Saved At": entry.created_at.strftime('%Y-%m-%d %H:%M:%S')
            }
            for entry in recent_entries
        ]

        # 表示
        st.dataframe(data, use_container_width=True)
        
    else:
        st.info("No molecule detected.")

# Streamlit cloud上では分子の2D描画が難しいため、画像生成はあきらめた
# with get_db_session() as db:
#     recent_entries = db.query(PkaData).order_by(PkaData.created_at.desc()).limit(5).all()

#     if recent_entries:
#         st.subheader(":material/Database: Recent pKa data (top 5)")

#         for entry in recent_entries:
#             cols = st.columns([2, 3, 2, 3])

#             mol_img = molecule_handler.smiles_to_image(entry.smiles)
#             if mol_img:
#                 cols[0].image(mol_img, use_container_width =True)
#             else:
#                 cols[0].warning("Invalid SMILES")

#             cols[1].markdown(f"**SMILES:** `{entry.smiles}`")
#             cols[2].markdown(f"**pKa:** {round(entry.predicted_pka, 2)}")
#             cols[3].markdown(f"**Saved At:** {entry.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
#             st.markdown("---")
#     else:
#         st.info("No molecule detected.")
