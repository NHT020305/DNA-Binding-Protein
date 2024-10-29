import numpy as np
import pandas as pd
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler


def normalize(data):
    min_max_scaler = MinMaxScaler()
    normalized_data = min_max_scaler.fit_transform(data)
    normalized_df = pd.DataFrame(normalized_data, columns=data.columns)
    return normalized_df

    
def standardize(data):
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(data)
    standardized_df = pd.DataFrame(standardized_data, columns=data.columns)
    return standardized_df

def robust_scaling(data):
    scaler = RobustScaler()
    robust_scaled_data = scaler.fit_transform(data)
    robus_scaled_df = pd.DataFrame(robust_scaled_data, columns=data.columns)
    return robus_scaled_df