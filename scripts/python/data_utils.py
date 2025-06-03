# data_utils.py

import os
import numpy as np
import pandas as pd
from scipy.io import loadmat
import torch
from sklearn.preprocessing import StandardScaler

def vectorize_fc(fc_mat, row_idx=None):
    """
    Vectorize FC matrix.

    Parameters:
        fc_mat (np.ndarray): Functional connectivity matrix (square or rectangular).
        row_idx (int or None): If provided, only return the column at this index 

    Returns:
        np.ndarray: 1D vector of connectivity values.
    """
    if row_idx is not None:
        # Return the full row (exclude self-connection if square)
        row = fc_mat[:, row_idx]
        return row.flatten()
    else:
        if fc_mat.shape[0] == fc_mat.shape[1]:
            return fc_mat[np.triu_indices(fc_mat.shape[0], k=1)]
        else:
            return fc_mat.flatten()

def load_all_fc_data(sub_cond_path, base_nifti_folder, mat_filename='conn_matrix.mat',
                     key_name='FC_matrix', row_idx=None):
    """
    Load and vectorize FC data from .mat files for all included subjects.

    Parameters:
        sub_cond_path (str): Path to subject condition Excel file.
        base_nifti_folder (str): Base directory containing subject/session folders.
        mat_filename (str): Name of .mat file inside each session folder.
        key_name (str): Key to access FC matrix inside .mat file.
        row_idx (int or None): If specified, only vectorize the row at this index in the FC matrix.

    Returns:
        Tuple of np.ndarray: (all_corr_data, all_tms_type, all_subject_id, all_stimloc, all_session)
    """
    SubInfo = pd.read_excel(sub_cond_path)
    Subs = SubInfo[SubInfo['Include'] == 1]['SubID'].tolist()
    sessions = ['D0', 'S1D1', 'S1D2', 'S2D1', 'S2D2', 'S3D1', 'S3D2']
    
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
    all_stimloc = []
    all_session = []

    for _, row in SubInfo.iterrows():
        subject_id = row['SubID']
        if row['Include'] != 1:
            continue
        stimloc = row['StimLoc']
        tms_types = order_map.get(row['StimOrder'], ['N'] * 7)
        for j, session in enumerate(sessions):
            mat_file = os.path.join(base_nifti_folder, subject_id, session, mat_filename)
            if os.path.exists(mat_file):
                matdat = loadmat(mat_file)
                dat_corr = matdat[key_name]
                dat_vec = vectorize_fc(dat_corr, row_idx=row_idx)
                all_corr_data.append(dat_vec)
                all_tms_type.append(tms_types[j])
                all_subject_id.append(subject_id)
                all_stimloc.append(stimloc)
                all_session.append(session)
            else:
                print(f"[WARN] File not found: {mat_file}")
    
    return np.array(all_corr_data), np.array(all_tms_type), np.array(all_subject_id), np.array(all_stimloc), np.array(all_session)


def preprocess_for_torch(X):
    # Remove columns with any NaN FIRST
    nan_cols = np.isnan(X).any(axis=0)
    if nan_cols.any():
        print(f"Columns with NaN: {nan_cols.sum()} / {X.shape[1]}")
        X = X[:, ~nan_cols]
    else:
        print("No NaNs found in features.")

    # Then scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Convert to torch tensor
    X_tensor = torch.tensor(X_scaled, dtype=torch.float32)
    print(f"Tensor shape after preprocessing: {X_tensor.shape}")
    return X_tensor

