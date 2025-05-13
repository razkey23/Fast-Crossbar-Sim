#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import os

try:
    import xbar_simulator
    print("Successfully imported xbar_simulator module")
except ImportError as e:
    print(f"Error importing xbar_simulator module: {e}")
    sys.exit(1)

def main():
    # Define crossbar dimensions
    M = 16  # Number of rows
    N = 16  # Number of columns
    
    print(f"Creating crossbar simulator with dimensions {M}x{N}")
    
    # Create the crossbar simulator
    xbar = xbar_simulator.CrossbarSimulator(M, N)
    
    # Generate a simple binary weight matrix (all cells set to 1)
    weights = np.ones((M, N), dtype=bool)
    print(f"Setting weights with shape {weights.shape}")
    
    # Set the weights in the crossbar
    xbar.set_weights(weights)
    
    # Generate a simple input vector (all 1s)
    input_vector = np.ones(M, dtype=float) * xbar_simulator.voltage_pulse_height
    print(f"Input vector shape: {input_vector.shape}")
    
    # Number of simulations to run
    num_simulations = 10000
    print(f"Running {num_simulations} simulations...")
    
    # Track the output current of a specific cell over all simulations
    # We'll track the cell at position (0, 0)
    cell_row = 0
    cell_col = 0
    
    # Store all output currents
    all_outputs = []
    
    # Run the simulations
    start_time = time.time()
    
    # Run multiple inferences
    outputs = xbar.run_multiple_inferences(input_vector, num_simulations)
    
    # Extract the current from the specific cell (column) for each simulation
    cell_currents = outputs[:, cell_col]
    
    end_time = time.time()
    total_time = end_time - start_time
    
    print(f"Total simulation time: {total_time:.2f} seconds")
    print(f"Average time per simulation: {total_time/num_simulations*1000:.2f} ms")
    
    # Calculate statistics
    avg_current = np.mean(cell_currents)
    std_current = np.std(cell_currents)
    min_current = np.min(cell_currents)
    max_current = np.max(cell_currents)
    
    print(f"\nCell ({cell_row}, {cell_col}) Current Statistics:")
    print(f"Average current: {avg_current:.6e} A")
    print(f"Standard deviation: {std_current:.6e} A")
    print(f"Minimum current: {min_current:.6e} A")
    print(f"Maximum current: {max_current:.6e} A")
    print(f"Drift percentage: {(max_current - min_current) / min_current * 100:.2f}%")
    
    # Plot the cell current over all simulations
    plt.figure(figsize=(10, 6))
    plt.plot(range(num_simulations), cell_currents)
    plt.title(f'Cell ({cell_row}, {cell_col}) Current Over {num_simulations} Simulations')
    plt.xlabel('Simulation Number')
    plt.ylabel('Current (A)')
    plt.grid(True)
    plt.savefig('cell_current_over_simulations.png')
    print("Current plot saved as 'cell_current_over_simulations.png'")
    
    # Plot the distribution of currents
    plt.figure(figsize=(10, 6))
    plt.hist(cell_currents, bins=50, alpha=0.7)
    plt.axvline(avg_current, color='r', linestyle='dashed', linewidth=2, label=f'Mean: {avg_current:.6e} A')
    plt.title(f'Distribution of Cell ({cell_row}, {cell_col}) Current Over {num_simulations} Simulations')
    plt.xlabel('Current (A)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    plt.savefig('cell_current_distribution.png')
    print("Current distribution saved as 'cell_current_distribution.png'")
    
    # Calculate moving average
    window_size = 100
    moving_avg = np.convolve(cell_currents, np.ones(window_size)/window_size, mode='valid')
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(moving_avg)), moving_avg)
    plt.title(f'Moving Average of Cell ({cell_row}, {cell_col}) Current (Window Size: {window_size})')
    plt.xlabel('Simulation Number')
    plt.ylabel('Moving Average Current (A)')
    plt.grid(True)
    plt.savefig('cell_current_moving_average.png')
    print("Moving average plot saved as 'cell_current_moving_average.png'")
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 