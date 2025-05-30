import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from tqdm import tqdm
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from scipy.spatial.distance import euclidean

def loss_function(recon_x, x, mu, logvar):
    # Reconstruction loss + KL divergence
    recon_loss = F.mse_loss(recon_x, x, reduction='sum')
    kl_div = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    total_loss = recon_loss + kl_div
    return total_loss, recon_loss, kl_div

def train_vae(model, dataloader, optimizer, device, epochs, label_encoder=LabelEncoder(), log_batch_info=False):
    model.to(device)
    model.train()

    total_loss_hist = []
    recon_loss_hist = []
    kld_loss_hist = []

    # Mapping to decode numeric condition to label
    inv_condition_map = {0: 'N', 1: 'S', 2: 'C'}

    for epoch in range(epochs):
        print(f"\n===== Epoch {epoch+1} =====")
        total_epoch_loss = 0
        recon_epoch_loss = 0
        kld_epoch_loss = 0

        for batch_idx, batch in enumerate(tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")):

            # Support both (x,) and (x, subj, cond) batches
            if len(batch) == 1:
                x = batch[0].to(device)
                subj_batch = None
                cond_batch = None
            else:
                x, subj_batch, cond_batch = batch
                x = x.to(device)
            
            # Optional: log batch info
            if log_batch_info and subj_batch is not None:
                subj_np = subj_batch.numpy()
                cond_np = cond_batch.numpy()
                for subj in np.unique(subj_np):
                    subj_mask = subj_np == subj
                    conds = cond_np[subj_mask]
                    cond_names = [inv_condition_map[int(c)] for c in conds]
                    subj_str = label_encoder.inverse_transform([subj])[0] if label_encoder else f"{subj:02d}"
                    print(f"  Subject {subj_str}: {cond_names}")

            optimizer.zero_grad()
            recon_batch, mu, logvar = model(x)
            loss, recon_loss, kl_div = loss_function(recon_batch, x, mu, logvar)
            loss.backward()
            optimizer.step()

            total_epoch_loss += loss.item()
            recon_epoch_loss += recon_loss.item()
            kld_epoch_loss += kl_div.item()

        total_loss_hist.append(total_epoch_loss / len(dataloader.dataset))
        recon_loss_hist.append(recon_epoch_loss / len(dataloader.dataset))
        kld_loss_hist.append(kld_epoch_loss / len(dataloader.dataset))

        print(f"Epoch {epoch+1}, Loss: {total_loss_hist[-1]:.4f}, Recon: {recon_loss_hist[-1]:.4f}, KL: {kld_loss_hist[-1]:.4f}")

    return total_loss_hist, recon_loss_hist, kld_loss_hist

def get_latent_mu(model, X_tensor, device='cpu', verbose=True):
    """
    Runs the VAE encoder on input data to extract latent mean vectors (mu).

    Parameters:
        model (nn.Module): Trained VAE model
        X_tensor (torch.Tensor): Input data of shape (n_samples, input_dim)
        device (str): 'cuda' or 'cpu'
        verbose (bool): Whether to print the output shape

    Returns:
        mu_all (np.ndarray): Latent means for all inputs (n_samples, latent_dim)
    """
    model.eval()
    mu_all = []

    with torch.no_grad():
        for i in range(X_tensor.shape[0]):
            x = X_tensor[i].unsqueeze(0).to(device)  # shape (1, input_dim)
            mu, _ = model.encode(x)
            mu_all.append(mu.cpu().numpy().flatten())

    mu_all = np.array(mu_all)

    if verbose:
        print("dim of latent mu:", mu_all.shape)

    return mu_all


def compute_condition_distances(mu_all, all_tms_type, all_subject_id, label_map={'N': 0, 'S': 1, 'C': 2}):
    """
    Computes per-subject latent vector distances between conditions using label mapping.

    Parameters:
    - mu_all: np.array, shape (n_samples, latent_dim)
    - all_tms_type: list or np.array, raw condition labels (e.g., 'N', 'S', 'C')
    - all_subject_id: list or np.array, subject identifiers
    - label_map: dict, e.g., {'N': 0, 'S': 1, 'C': 2}

    Returns:
    - pd.DataFrame with per-subject distances:
        ["subject", "d_null_sham", "d_null_real", "diff_real_minus_sham"]
    """
    y = np.array([label_map[t] for t in all_tms_type])
    results = []
    unique_subjects = np.unique(all_subject_id)

    for subj in unique_subjects:
        idx = all_subject_id == subj
        mu_subj = mu_all[idx]
        y_subj = y[idx]

        # Compute mean latent vectors per condition
        means = {}
        for cond in [0, 1, 2]:
            cond_mu = mu_subj[y_subj == cond]
            if len(cond_mu) > 0:
                means[cond] = np.mean(cond_mu, axis=0)

        # Only compute distances if all 3 conditions are present
        if all(c in means for c in [0, 1, 2]):
            d_0_1 = euclidean(means[0], means[1])
            d_0_2 = euclidean(means[0], means[2])
            results.append({
                "subject": subj,
                "d_null_sham": d_0_1,
                "d_null_real": d_0_2,
                "diff_real_minus_sham": d_0_2 - d_0_1
            })

    return pd.DataFrame(results)