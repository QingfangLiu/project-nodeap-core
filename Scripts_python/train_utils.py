import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from tqdm import tqdm

def loss_function(recon_x, x, mu, logvar):
    # Reconstruction loss + KL divergence
    recon_loss = F.mse_loss(recon_x, x, reduction='sum')
    kl_div = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    total_loss = recon_loss + kl_div
    return total_loss, recon_loss, kl_div

def train_vae(model, dataloader, optimizer, device='cpu', epochs=20):
    model.to(device)
    model.train()
    total_loss_hist = []
    recon_loss_hist = []
    kld_loss_hist = []

    for epoch in range(epochs):
        total_epoch_loss = 0
        recon_epoch_loss = 0
        kld_epoch_loss = 0

        for batch in tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}"):
            x = batch[0].to(device)
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