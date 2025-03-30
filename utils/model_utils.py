import joblib
import pandas as pd
import logging
from typing import Union
from sklearn.base import BaseEstimator

from utils.molecule_handler import (
    descriptor_names,
    canonicalize_smiles,
    compute_all_descriptors
)


def load_model(model_path: str) -> BaseEstimator:
    try:
        model = joblib.load(model_path)
        logging.info(f"Model loaded from: {model_path}")
        return model
    except Exception as e:
        logging.error(f"Failed to load model from {model_path}: {e}")
        return None


def predict_pka(smiles: str, model: BaseEstimator) -> Union[float, str]:
    try:
        smiles = canonicalize_smiles(smiles)
        descriptors = compute_all_descriptors(smiles)
        descriptors_df = pd.DataFrame([descriptors], columns=descriptor_names)
        prediction = model.predict(descriptors_df)[0]
        logging.info(f"pKa predicted for SMILES {smiles}: {prediction:.2f}")
        return prediction

    except Exception as e:
        logging.error(f"Error during pKa prediction for SMILES '{smiles}': {e}")
        return None
