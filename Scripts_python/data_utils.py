# data_utils.py

import os
import numpy as np
import pandas as pd
from scipy.io import loadmat
import torch
from sklearn.preprocessing import StandardScaler

def vectorize_fc(fc_mat):
    """Return upper triangle of FC matrix (excluding diagonal)."""
    return fc_mat[np.triu_indices(fc_mat.shape[0], k=1)]

def load_all_fc_data(sub_cond_path, base_nifti_folder):
    """Load and vectorize FC data, returning data, labels, and IDs."""
    SubInfo = pd.read_excel(sub_cond_path)
    Subs = SubInfo[SubInfo['Include'] == 1]['SubID'].tolist()
    sessions = ['D0', 'S1D1', 'S1D2', 'S2D1', 'S2D2', 'S3D1', 'S3D2']
    
    # Map StimOrder to TMS type sequence
    order_map = {
        123: ['N', 'C', 'S', 'S', 'C', 'S', 'S'],
        132: ['N', 'C', 'S', 'S', 'S', 'S', 'C'],
        213: ['N', 'S', 'C', 'C', 'S', 'S', 'S'],
        231: ['N', 'S', 'C', 'S', 'S', 'C', 'S'],
        312: ['N', 'S', 'S', 'C', 'S', 'S', 'C'],
        321: ['N', 'S', 'S', 'S', 'C', 'C', 'S'],
    }

    all_corr_data = []
    all_tms_type = []
    all_subject_id = []

    for _, row in SubInfo.iterrows():
        subject_id = row['SubID']
        if row['Include'] != 1:
            continue
        tms_types = order_map.get(row['StimOrder'], ['N'] * 7)
        for j, session in enumerate(sessions):
            mat_file = os.path.join(base_nifti_folder, subject_id, session, 'conn_matrix.mat')
            if os.path.exists(mat_file):
                matdat = loadmat(mat_file)
                dat_corr = matdat['FC_matrix']
                dat_vec = vectorize_fc(dat_corr)
                all_corr_data.append(dat_vec)
                all_tms_type.append(tms_types[j])
                all_subject_id.append(subject_id)
            else:
                print(f"[WARN] File not found: {mat_file}")
    
    return np.array(all_corr_data), all_tms_type, all_subject_id

def preprocess_for_torch(X):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # remove columns with any NaN
    nan_cols = np.isnan(X_scaled).any(axis=0)
    if nan_cols.any():
        print(f"Columns with NaN: {nan_cols.sum()} / {X_scaled.shape[1]}")
        X_scaled = X_scaled[:, ~nan_cols]
    else:
        print("No NaNs found in features.")

    # convert to torch tensor
    X_tensor = torch.tensor(X_scaled, dtype=torch.float32)
    return X_tensor
