#!/usr/bin/env python3
import numpy as np
import struct
import matplotlib.pyplot as plt
import os
import argparse
from matplotlib.colors import LinearSegmentedColormap

def read_binary_matrix(filename):
    """
    Read a matrix from a binary file.
    
    Format:
    - rows (int64_t)
    - cols (int64_t)
    - matrix data in row-major order (float for each element)
    """
    with open(filename, 'rb') as f:
        # Read dimensions
        rows = struct.unpack('q', f.read(8))[0]  # int64_t for rows
        cols = struct.unpack('q', f.read(8))[0]  # int64_t for cols
        
        # Read matrix data
        matrix = np.zeros((rows, cols), dtype=np.float32)
        for i in range(rows):
            matrix[i, :] = struct.unpack(f'{cols}f', f.read(cols * 4))
        
        return matrix

def main():
    parser = argparse.ArgumentParser(description='Analyze drift in MAC outputs after multiple inferences')
    parser.add_argument('--input_dir', type=str, default='clean_xbar_model/input_data', help='Input directory containing output files')
    
    args = parser.parse_args()
    
    # Path to MAC output file
    mac_path = os.path.join(args.input_dir, 'mac_iterations.bin')
    
    # Read MAC output
    try:
        mac_data = read_binary_matrix(mac_path)
        print(f"MAC output shape: {mac_data.shape}")
        print(f"Number of iterations: {mac_data.shape[0]}")
        print(f"Number of columns: {mac_data.shape[1]}")
        
        # Calculate relative differences between consecutive iterations
        relative_diffs = np.zeros_like(mac_data)
        for i in range(1, mac_data.shape[0]):
            relative_diffs[i] = (mac_data[i] - mac_data[i-1]) / mac_data[0] * 100  # Percentage change
        
        # Calculate cumulative drift from first iteration
        cumulative_drift = np.zeros_like(mac_data)
        for i in range(1, mac_data.shape[0]):
            cumulative_drift[i] = (mac_data[i] - mac_data[0]) / mac_data[0] * 100  # Percentage change
        
        # Find columns with significant drift
        max_drift = np.max(np.abs(cumulative_drift[-1, :]))
        significant_cols = np.where(np.abs(cumulative_drift[-1, :]) > max_drift * 0.5)[0]
        
        print(f"Maximum cumulative drift: {max_drift:.4f}%")
        print(f"Columns with significant drift: {significant_cols}")
        
        # Visualize cumulative drift heatmap
        plt.figure(figsize=(12, 8))
        plt.imshow(cumulative_drift, cmap='coolwarm', aspect='auto', vmin=-max_drift, vmax=max_drift)
        plt.colorbar(label='Drift (%)')
        plt.title('Cumulative Drift in MAC Output Over Iterations')
        plt.xlabel('Column Index')
        plt.ylabel('Iteration')
        plt.savefig(os.path.join(args.input_dir, 'drift_heatmap.png'))
        print(f"Drift heatmap saved to {os.path.join(args.input_dir, 'drift_heatmap.png')}")
        
        # Visualize drift for significant columns
        plt.figure(figsize=(12, 8))
        for col in significant_cols[:5]:  # Plot up to 5 significant columns
            plt.plot(range(mac_data.shape[0]), cumulative_drift[:, col], label=f'Column {col}')
        
        plt.title('Drift in MAC Output for Significant Columns')
        plt.xlabel('Iteration')
        plt.ylabel('Drift (%)')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(args.input_dir, 'drift_significant_columns.png'))
        print(f"Drift plot saved to {os.path.join(args.input_dir, 'drift_significant_columns.png')}")
        
        # Visualize average drift across all columns
        plt.figure(figsize=(12, 6))
        avg_drift = np.mean(np.abs(cumulative_drift), axis=1)
        plt.plot(range(mac_data.shape[0]), avg_drift)
        plt.title('Average Absolute Drift Across All Columns')
        plt.xlabel('Iteration')
        plt.ylabel('Average Absolute Drift (%)')
        plt.grid(True)
        plt.savefig(os.path.join(args.input_dir, 'average_drift.png'))
        print(f"Average drift plot saved to {os.path.join(args.input_dir, 'average_drift.png')}")
        
        # Visualize distribution of final drift values
        plt.figure(figsize=(12, 6))
        plt.hist(cumulative_drift[-1, :], bins=30, alpha=0.7)
        plt.title('Distribution of Cumulative Drift After All Iterations')
        plt.xlabel('Drift (%)')
        plt.ylabel('Count')
        plt.grid(True)
        plt.savefig(os.path.join(args.input_dir, 'drift_distribution.png'))
        print(f"Drift distribution saved to {os.path.join(args.input_dir, 'drift_distribution.png')}")
        
    except Exception as e:
        print(f"Error analyzing MAC output: {e}") 

def plot_cell_conductance(data_dir="data", cell_idx=5, save_path=None):
    """
    Plot the conductance of a single cell over time.
    
    Args:
        data_dir: Directory containing the mac_iterations.bin file
        cell_idx: Index of the cell to plot
        save_path: Optional path to save the figure
    """
    # Path to MAC output file
    mac_path = os.path.join(data_dir, 'mac_iterations.bin')
    
    # Read MAC output data
    mac_data = read_binary_matrix(mac_path)
    
    print(f"MAC output shape: {mac_data.shape}")
    print(f"Number of iterations: {mac_data.shape[0]}")
    print(f"Number of columns: {mac_data.shape[1]}")
    
    # Check if cell_idx is valid
    if cell_idx >= mac_data.shape[1]:
        print(f"Error: Cell index {cell_idx} exceeds the number of columns {mac_data.shape[1]}")
        return
    
    # Extract conductance for the specified cell
    cell_conductance = mac_data[:, cell_idx]
    
    # Calculate relative change from the initial conductance
    initial_conductance = cell_conductance[0]
    relative_change = (cell_conductance - initial_conductance) / initial_conductance * 100
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Absolute conductance values
    ax1.plot(range(len(cell_conductance)), cell_conductance)
    ax1.set_title(f'Conductance of Cell {cell_idx} Over Time')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Conductance (S)')
    ax1.grid(True)
    
    # Plot 2: Relative change in conductance
    ax2.plot(range(len(relative_change)), relative_change)
    ax2.set_title(f'Relative Change in Conductance of Cell {cell_idx}')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Relative Change (%)')
    ax2.grid(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        print(f"Figure saved to {save_path}")
    
    # Find maximum drift
    max_drift = np.max(np.abs(relative_change))
    final_drift = relative_change[-1]
    
    print(f"Initial conductance of cell {cell_idx}: {initial_conductance:.6e} S")
    print(f"Final conductance of cell {cell_idx}: {cell_conductance[-1]:.6e} S")
    print(f"Maximum drift: {max_drift:.4f}%")
    print(f"Final drift: {final_drift:.4f}%")
    
    return fig

if __name__ == "__main__":
    main()
    
    # Plot conductance for a significant cell (cell with high drift)
    fig = plot_cell_conductance(data_dir="data", cell_idx=5, 
                          save_path="data/cell_conductance_drift.png")
    
    # Also plot another cell for comparison
    fig2 = plot_cell_conductance(data_dir="data", cell_idx=10,
                          save_path="data/cell_conductance_drift_2.png") 