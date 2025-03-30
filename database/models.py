import pandas as pd
from sqlalchemy import Column, Integer, String, DateTime, Text, Float, Boolean

from .database import Base

class PkaData(Base):
    __tablename__ = "pkadata"

    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(String, index=True)
    predicted_pka = Column(Float, index=True)
    created_at = Column(
        DateTime, default=pd.Timestamp.now(tz="UTC").tz_convert("Asia/Tokyo").floor("s")
    )
