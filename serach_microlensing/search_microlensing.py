"""Automated Transient Detection and Ranking in Light Curve Data

This script processes a directory of astronomical light curve data files to 
identify statistically significant transient events. It computes a rolling 
geometric mean of absolute z-scores to isolate anomalous flux deviations while 
filtering out noisy data points.

Prerequisites & Environment Setup:
    1. Data Directory: The script expects a directory named 'light_curves_data/' 
       in the same path, containing light curve files in '.data' format.
    2. Data Format: Each '.data' file must have at least 3 columns:
       - Column 0: Time (Heliocentric Julian Date)
       - Column 1: Magnitude
       - Column 2: Magnitude Error

Usage:
    python rank_lightcurves.py

    The script will evaluate all files, rank them based on the maximum rolling 
    z-score, and automatically save the top 10 anomaly plots to a newly created 
    'top10_results' directory.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt


def plot_and_save_lightcurve(filename, t_good, f_good, f_err_good, weighted_zf, idx_max, save_path=None):
    """Generates an improved two-panel plot for the light curve and z-score.

    Args:
        filename (str): Original file path for title reference.
        t_good (np.ndarray): Filtered time array.
        f_good (np.ndarray): Filtered magnitude array.
        f_err_good (np.ndarray): Filtered magnitude error array.
        weighted_zf (np.ndarray): Calculated rolling z-score array.
        idx_max (int): Index of the maximum anomaly.
        save_path (str, optional): Path to save the figure. If None, displays it.
    """
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)
    
    # Top Panel: Light Curve (Magnitude)
    axs[0].errorbar(
        t_good, f_good, yerr=f_err_good, fmt='o', markersize=3, 
        color='black', ecolor='gray', linestyle='none', label='Observation'
    )
    axs[0].scatter(
        t_good[idx_max], f_good[idx_max], 
        color='red', s=60, zorder=5, label='Detected Anomaly'
    )
    axs[0].set_ylabel("Magnitude")
    axs[0].invert_yaxis()  # Astronomical magnitudes are inverted
    axs[0].set_title("Light Curve")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()
    
    # Bottom Panel: Rolling Z-Score
    axs[1].scatter(t_good, weighted_zf, s=10, color='blue')
    axs[1].axvline(t_good[idx_max], color='red', linestyle='--', alpha=0.5)
    axs[1].set_ylabel("Rolling Weighted Z-Score")
    axs[1].set_xlabel("Time (Heliocentric Julian Date)")
    axs[1].set_title("Anomaly Score")
    axs[1].grid(True, alpha=0.3)
    
    fig.suptitle(f"Transient Candidate: {os.path.basename(filename)}", fontsize=14)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


def process_lightcurve(filename, window=9):
    """Extracts data, applies quality filters, and calculates the rolling z-score.

    Args:
        filename (str): Path to the light curve data file.
        window (int): Size of the rolling window (must be odd).

    Returns:
        tuple: Contains valid data arrays and anomaly metrics, or None if invalid.
    """
    data = np.loadtxt(filename)
    t = data[:, 0]
    f = data[:, 1]
    f_err = data[:, 2]

    # Filter out NaNs, infinities, and invalid errors
    good = np.isfinite(f) & np.isfinite(f_err) & (f_err > 0)
    
    t_good = t[good]
    f_good = f[good]
    f_err_good = f_err[good]
    
    if len(f_good) < 300:
        return None

    w = 1.0 / (f_err_good**2)
    f0 = np.sum(w * f_good) / np.sum(w)
    
    # Calculate absolute z-score (corrected formula)
    zf = np.abs(f_good - f0) / f_err_good

    N = len(zf)
    weighted_zf = np.full(N, np.nan, dtype=float)
    
    k = window // 2
    for j in range(k, N - k):
        weighted_zf[j] = np.prod(zf[j-k : j+k+1]) ** (1 / window)
    
    std_zf = np.nanstd(weighted_zf)
    max_zscore = np.nanmax(weighted_zf)
    idx_max = np.nanargmax(weighted_zf)

    # Filtering condition
    if f_good[idx_max] < 0 and std_zf < 3:
        return t_good, f_good, f_err_good, weighted_zf, max_zscore, idx_max
    
    return None


if __name__ == "__main__":
    file_path = "light_curves_data/"
    output_dir = "top10_results/"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    files_full = glob.glob(os.path.join(file_path, '*.data'))
    
    results = []

    print("Analyzing light curves...")
    for i, file_full in enumerate(files_full):
        print(f"\rProcessing file {i+1}/{len(files_full)}", end="", flush=True)
        
        processed_data = process_lightcurve(file_full)
        if processed_data is not None:
            _, _, _, _, max_zscore, _ = processed_data
            results.append((file_full, max_zscore))
            
    print("\nAnalysis complete.")

    # Sort results by max_zscore descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    # Extract Top 10
    top_10 = results[:10]
    
    print("\n--- Generating Top 10 Plots ---")
    for rank, (file_full, score) in enumerate(top_10, start=1):
        print(f"Rank {rank}: {os.path.basename(file_full)} (Score: {score:.2f})")
        
        # Re-process to get plotting arrays (computationally cheap)
        t_good, f_good, f_err_good, weighted_zf, _, idx_max = process_lightcurve(file_full)
        
        save_name = os.path.join(output_dir, f"rank_{rank:02d}_{os.path.basename(file_full)}.png")
        plot_and_save_lightcurve(file_full, t_good, f_good, f_err_good, weighted_zf, idx_max, save_path=save_name)
        
    print(f"\nTop 10 plots successfully saved to '{output_dir}'.")