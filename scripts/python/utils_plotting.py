import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


def plot_vae_losses(total_loss_hist, recon_loss_hist, kld_loss_hist, figsize=(10, 6), save_path=None):
    """
    Plots VAE training loss curves over epochs.

    Parameters:
        total_loss_hist (list): Total loss per epoch.
        recon_loss_hist (list): Reconstruction loss per epoch.
        kld_loss_hist (list): KL divergence loss per epoch.
        figsize (tuple): Size of the figure (default: (10, 6)).
        save_path (str or None): If specified, saves the plot to this path.
    """
    plt.figure(figsize=figsize)
    plt.plot(total_loss_hist, label='Total Loss')
    plt.plot(recon_loss_hist, label='Reconstruction Loss')
    plt.plot(kld_loss_hist, label='KL Divergence')

    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('VAE Loss Components Over Epochs')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    else:
        plt.show()


def plot_latent_embedding_by_condition_and_subject(z_2d, all_tms_type, all_subject_id, label_map={'N': 0, 'S': 1, 'C': 2},
                                                    markers=['o', 's', '*'], figsize=(8, 8), title='VAE Latent Space: Subject Color, Condition Marker'):
    """
    Plots a 2D embedding (e.g., t-SNE or PCA) of latent space with:
    - Different markers for each condition (Null, Sham, Real)
    - Different colors for each subject

    Parameters:
        z_2d (ndarray): 2D array of shape (n_samples, 2), e.g. t-SNE or PCA output
        all_tms_type (list or array): List of TMS conditions, e.g., ['N', 'S', 'C', ...]
        all_subject_id (list or array): List of subject IDs corresponding to each row in z_2d
        label_map (dict): Mapping from condition labels to integers (default: {'N': 0, 'S': 1, 'C': 2})
        markers (list): List of marker styles for each condition (default: ['o', 's', '*'])
        figsize (tuple): Size of the figure
        title (str): Plot title
    """
    # Encode condition labels to integers
    y = np.array([label_map[t] for t in all_tms_type])
    unique_conds = np.unique(y)
    cond_to_marker = dict(zip(unique_conds, markers))

    # Encode subjects to unique color indices
    unique_subjects, subject_idx = np.unique(all_subject_id, return_inverse=True)
    palette = sns.color_palette("husl", len(unique_subjects))
    subject_colors = np.array([palette[i] for i in subject_idx])

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    for cond in unique_conds:
        mask = y == cond
        ax.scatter(
            z_2d[mask, 0],
            z_2d[mask, 1],
            marker=cond_to_marker[cond],
            s=30,
            c=subject_colors[mask],
            label=f'Condition {cond}',
            edgecolors='k',
            alpha=0.5
        )

    # Custom legend
    handles = [
        Line2D([0], [0], marker=cond_to_marker[c], color='gray', linestyle='',
               label={0: 'Null', 1: 'Sham', 2: 'Real'}[c], markersize=8)
        for c in unique_conds
    ]
    ax.legend(handles=handles, title='TMS Condition', loc='upper right')

    ax.set_title(title)
    ax.set_xlabel("t-SNE dim 1")
    ax.set_ylabel("t-SNE dim 2")
    ax.grid(True)
    plt.tight_layout()
    plt.show()





def plot_subject_distance_comparison(df):
    """
    Plots paired Euclidean distance comparisons per subject between 'd_null_sham' and 'd_null_real'.

    Parameters:
    - d: pd.DataFrame with columns ['subject', 'd_null_sham', 'd_null_real']
    """

    # Melt the DataFrame to long format
    df_long = pd.melt(
        df,
        id_vars=["subject"],
        value_vars=["d_null_sham", "d_null_real"],
        var_name="Comparison",
        value_name="Distance"
    )

    # Plot
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df_long, x="Comparison", y="Distance", color="lightgray", showfliers=False)

    # Overlay paired lines per subject
    for subj in df["subject"]:
        subj_data = df_long[df_long["subject"] == subj]
        plt.plot(
            subj_data["Comparison"],
            subj_data["Distance"],
            marker='o',
            linestyle='-',
            alpha=0.6,
            color='blue'
        )

    plt.title("Paired Distance Comparison per Subject")
    plt.ylabel("Euclidean Distance")
    plt.xlabel("")
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


