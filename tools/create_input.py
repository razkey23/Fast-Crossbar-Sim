#!/usr/bin/env python3
import numpy as np
import struct
import os
import argparse

def write_binary_matrix(matrix, filename):
    """
    Write a matrix to a binary file in the format expected by the crossbar simulator.
    
    Format:
    - rows (int64_t)
    - cols (int64_t)
    - matrix data in row-major order (int64_t for each element)
    """
    rows, cols = matrix.shape
    with open(filename, 'wb') as f:
        # Write dimensions
        f.write(struct.pack('q', rows))  # int64_t for rows
        f.write(struct.pack('q', cols))  # int64_t for cols
        
        # Write matrix data
        for i in range(rows):
            f.write(struct.pack('q' * cols, *[int(x) for x in matrix[i, :]]))

def main():
    parser = argparse.ArgumentParser(description='Create input files for crossbar simulator')
    parser.add_argument('--input_rows', type=int, default=10, help='Number of input vectors')
    parser.add_argument('--weight_rows', type=int, default=32, help='Number of weight rows (crossbar height)')
    parser.add_argument('--weight_cols', type=int, default=32, help='Number of weight columns (crossbar width)')
    parser.add_argument('--output_dir', type=str, default='input_data', help='Output directory')
    parser.add_argument('--input_sparsity', type=float, default=0.5, help='Sparsity of input matrix (probability of zeros)')
    parser.add_argument('--weight_sparsity', type=float, default=0.5, help='Sparsity of weight matrix (probability of zeros)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate input matrix - binary matrix (0 or 1)
    # Each row is an input vector, columns match the number of rows in the weight matrix
    input_matrix = np.random.choice([0, 1], 
                                   size=(args.input_rows, args.weight_rows), 
                                   p=[args.input_sparsity, 1-args.input_sparsity])
    
    # Generate weight matrix - binary matrix (0 or 1)
    weight_matrix = np.random.choice([0, 1], 
                                    size=(args.weight_rows, args.weight_cols), 
                                    p=[args.weight_sparsity, 1-args.weight_sparsity])
    
    # Write matrices to binary files
    input_path = os.path.join(args.output_dir, 'input.bin')
    weight_path = os.path.join(args.output_dir, 'weight.bin')
    
    write_binary_matrix(input_matrix, input_path)
    write_binary_matrix(weight_matrix, weight_path)
    
    print(f"Generated input file: {input_path}")
    print(f"Input matrix shape: {input_matrix.shape}")
    print(f"Generated weight file: {weight_path}")
    print(f"Weight matrix shape: {weight_matrix.shape}")
    
    # Print a preview of the matrices
    print("\nInput matrix preview (first 5 rows, first 10 columns):")
    print(input_matrix[:min(5, args.input_rows), :min(10, args.weight_rows)])
    
    print("\nWeight matrix preview (first 5 rows, first 10 columns):")
    print(weight_matrix[:min(5, args.weight_rows), :min(10, args.weight_cols)])

if __name__ == "__main__":
    main() 