import os
import streamlit as st
import logging
from dotenv import load_dotenv

from database.database import create_tables


def initialize_setting():
    """
    環境設定を初期化する関数。

    .envファイルを読み込み、必要なディレクトリを作成し、
    Streamlitのsession_stateにパスを設定します。

    Raises:
        OSError: ディレクトリの作成に失敗した場合に発生します。
    """
    load_dotenv(os.path.join(os.path.dirname(__file__), "..", ".env"))

    if "is_initialized" not in st.session_state:
        st.session_state.is_initialized = False

    if not st.session_state.is_initialized:
            # ログディレクトリ作成
        log_path = os.getenv("LOG_PATH")
        os.makedirs(log_path, exist_ok=True)
        log_file = os.path.join(log_path, "app.log")
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )

        logging.info("Log file created at %s", log_file)
        
        # DB初期化
        create_tables()

        st.session_state.is_initialized = True

    # CSSでst.toastのスタイルをカスタマイズ
    st.markdown(
        """
        <style>
        .stToast {
            background-color: #64748b !important;
            color: white !important;
            border-radius: 10px;
            padding: 18px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )