#!/usr/bin/env python3
import numpy as np
import time
import sys
import os
import multiprocessing
import xbar_simulator

# Set number of threads based on CPU cores
num_threads = multiprocessing.cpu_count()
os.environ['OMP_NUM_THREADS'] = str(num_threads)
xbar_simulator.set_eigen_threads(num_threads)

print(f"CPU cores available: {num_threads}")
print(f"OpenMP threads: {os.environ.get('OMP_NUM_THREADS')}")
print(f"Eigen threads: {xbar_simulator.get_eigen_threads()}")

def measure_inference_time(matrix_size, num_threads):
    # Set number of threads for OpenMP
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    
    # Set number of threads for Eigen
    xbar_simulator.set_eigen_threads(num_threads)
    
    # Create crossbar simulator
    xbar = xbar_simulator.CrossbarSimulator(matrix_size, matrix_size)
    
    # Generate random binary weights
    weights = np.random.choice([True, False], size=(matrix_size, matrix_size))
    
    # Set the weights in the crossbar
    xbar.set_weights(weights)
    
    # Generate input vector with voltage pulse height
    input_vector = np.ones(matrix_size, dtype=float) * xbar_simulator.voltage_pulse_height
    
    # Measure inference time
    start_time = time.time()
    output = xbar.run_inference(input_vector)
    end_time = time.time()
    
    inference_time = (end_time - start_time) * 1000  # Convert to milliseconds
    return inference_time


def load_vector(filename):
    """Load input vector from file with size header."""
    with open(filename, 'r') as f:
        size = int(f.readline().strip())
        vector = np.zeros(size)
        for i in range(size):
            vector[i] = float(f.readline().strip())
    return vector

def load_matrix(filename):
    """Load weight matrix from file with dimension header."""
    with open(filename, 'r') as f:
        # Read dimensions from first line
        rows, cols = map(int, f.readline().strip().split())
        # Read the matrix data
        matrix = np.zeros((rows, cols), dtype=bool)
        for i in range(rows):
            row = list(map(int, f.readline().strip().split()))
            matrix[i] = row
    return matrix

def generate_random_inputs(size):
    """Generate random binary input vector and weight matrix of given size."""
    input_vector = np.random.choice([0, 1], size=size)
    weight_matrix = np.random.choice([True, False], size=(size, size))
    return input_vector, weight_matrix

def run_benchmark(size, method="fixed-point", num_runs=1):
    # Set thread counts
    num_threads = multiprocessing.cpu_count()
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    xbar_simulator.set_eigen_threads(num_threads)
    
    # Create simulator
    xbar = xbar_simulator.CrossbarSimulator(size, size)
    xbar.set_parasitic_resistances(3, 1e20, 1e20, 5, 3, 2)
    times = []
    
    for _ in range(num_runs):
        input_vector = np.random.choice([0, 1], size=size)
        weight_matrix = np.random.choice([True, False], size=(size, size))
        xbar.set_weights(weight_matrix)
        input_scaled = input_vector * xbar_simulator.voltage_pulse_height
        
        start_time = time.time()
        _ = xbar.run_inference(input_scaled)
        times.append((time.time() - start_time) * 1000)
    
    return np.min(times), np.max(times), np.mean(times), np.std(times)

def main():
    methods = ["fixed-point"]
    sizes = [32, 64, 128,256]
    
    for method in methods:
        print(f"\n=== Testing {method} method ===")
        print(f"{'Size':>6} | {'Min (ms)':>10} | {'Max (ms)':>10} | {'Mean (ms)':>10} | {'Std (ms)':>10}")
        print("-" * 62)
        
        for size in sizes:
            min_t, max_t, mean_t, std_t = run_benchmark(size, method)
            print(f"{size:6d} | {min_t:10.2f} | {max_t:10.2f} | {mean_t:10.2f} | {std_t:10.2f}")

if __name__ == "__main__":
    main() 