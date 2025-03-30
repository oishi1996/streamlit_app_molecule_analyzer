from typing import Optional, List
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from PIL.Image import Image
import logging

descriptor_names: List[str] = [desc[0] for desc in Descriptors.descList]

def canonicalize_smiles(smiles: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Invalid SMILES: {smiles}")
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logging.error(f"Failed to canonicalize SMILES: {smiles}, Error: {e}")
        return None


def compute_all_descriptors(smiles: str) -> Optional[List[float]]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Invalid SMILES for descriptor calculation: {smiles}")
            return None
        return [getattr(Descriptors, desc)(mol) for desc in descriptor_names]
    except Exception as e:
        logging.error(f"Descriptor calculation failed for SMILES: {smiles}, Error: {e}")
        return None


def smiles_to_image(smiles: str, size: tuple = (100, 100)) -> Optional[Image]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Invalid SMILES for image generation: {smiles}")
            return None
        img = Draw.MolToImage(mol, size=size)
        return img
    except Exception as e:
        logging.error(f"Failed to generate image for SMILES: {smiles}, Error: {e}")
        return None
