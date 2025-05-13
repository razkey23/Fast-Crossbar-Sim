# Python Wrapper for Memristor Crossbar Simulator

This directory contains a Python wrapper for the memristor crossbar simulator using PyBind11.

## Prerequisites

- CMake (version 3.10 or higher)
- Python 3 with development headers
- C++17 compatible compiler
- NumPy
- Matplotlib (for the example script)

## Building the Python Module

1. First, clone the PyBind11 repository in the project root:

```bash
git clone https://github.com/pybind/pybind11.git
```

2. Build the Python module:

```bash
mkdir build
cd build
cmake ..
make
```

This will create the `xbar_simulator` Python module in the `python` directory.

## Using the Python Module

The wrapper provides a simple interface to the C++ crossbar simulator:

```python
import numpy as np
import xbar_simulator

# Create a crossbar simulator with specified dimensions
xbar = xbar_simulator.CrossbarSimulator(M=32, N=32)

# Set weights (binary matrix)
weights = np.random.choice([0, 1], size=(32, 32)).astype(bool)
xbar.set_weights(weights)

# Run a single inference with an input vector
input_vector = np.random.choice([0, 1], size=32).astype(float)
output = xbar.run_inference(input_vector)

# Run multiple inferences to observe drift
outputs = xbar.run_multiple_inferences(input_vector, num_inferences=100)
```

## Example Script

An example script `example.py` is provided that demonstrates how to use the wrapper:

```bash
cd python
python example.py
```

This will:
1. Create a crossbar simulator with random weights
2. Run a single inference
3. Run multiple inferences to observe drift
4. Generate plots showing the drift in output currents

## API Reference

### `CrossbarSimulator(M, N)`

Creates a new crossbar simulator with M rows and N columns.

### `set_weights(weights)`

Sets the weights of the crossbar. `weights` should be a 2D NumPy array of booleans with shape (M, N).

### `run_inference(input, dt=1e-6)`

Runs a single inference with the given input vector. `input` should be a 1D NumPy array with length M.
Returns a 1D NumPy array with length N containing the output currents.

### `run_multiple_inferences(input, num_inferences, dt=1e-6)`

Runs multiple inferences with the same input vector. Returns a 2D NumPy array with shape (num_inferences, N)
containing the output currents for each inference.

### `set_parasitic_resistances(Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl)`

Sets the parasitic resistances of the crossbar. 