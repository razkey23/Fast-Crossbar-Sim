# Python Wrapper for Memristor Crossbar Simulator

This repository contains a Python wrapper for the memristor crossbar simulator using PyBind11. The wrapper allows you to use the C++ crossbar simulator from Python, making it easier to integrate with data analysis and visualization tools.

## Overview

The crossbar simulator models a memristor-based crossbar array for neuromorphic computing applications. It simulates the behavior of memristors, including nonlinear effects, parasitic resistances, and conductance drift over time.

With this Python wrapper, you can:
- Create crossbar simulators with custom dimensions
- Set binary weights in the crossbar
- Apply input vectors and get output currents
- Run multiple inferences to observe drift
- Analyze and visualize the results using Python's data science ecosystem

## Directory Structure

```
.
├── CMakeLists.txt           # CMake configuration for building the Python module
├── python/                  # Python module directory
│   ├── example.py           # Simple example script
│   ├── drift_analysis.py    # Advanced drift analysis script
│   ├── xbar_simulator_demo.ipynb  # Jupyter notebook demonstration
│   ├── README.md            # Python-specific documentation
│   └── setup.py             # Python setup script for installation
├── src/                     # Source code directory
│   ├── core/                # Core simulator code
│   ├── crossbar_model/      # Crossbar model implementation
│   ├── memristor_model/     # Memristor model implementation
│   └── python/              # Python wrapper code
│       └── xbar_wrapper.cpp # PyBind11 wrapper implementation
└── README_PYTHON.md         # This file
```

## Installation

### Prerequisites

- CMake (version 3.10 or higher)
- Python 3 with development headers
- C++17 compatible compiler
- NumPy
- Matplotlib (for visualization)

### Building from Source

1. Clone the PyBind11 repository in the project root:

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

### Installing via pip

You can also install the module using pip:

```bash
cd python
pip install -e .
```

## Usage

Here's a simple example of how to use the Python wrapper:

```python
import numpy as np
import xbar_simulator

# Create a crossbar simulator with specified dimensions
M, N = 32, 32
xbar = xbar_simulator.CrossbarSimulator(M, N)

# Set random binary weights
weights = np.random.choice([0, 1], size=(M, N)).astype(bool)
xbar.set_weights(weights)

# Create an input vector
input_vector = np.random.choice([0, 1], size=M).astype(float)
input_vector = input_vector * xbar_simulator.voltage_pulse_height

# Run a single inference
output = xbar.run_inference(input_vector)
print(f"Output shape: {output.shape}")
print(f"First few output values: {output[:5]}")

# Run multiple inferences to observe drift
num_inferences = 100
outputs = xbar.run_multiple_inferences(input_vector, num_inferences)
print(f"Outputs shape: {outputs.shape}")
```

For more advanced usage, see the example scripts and Jupyter notebook in the `python` directory.

## Example Scripts

- `example.py`: A simple example that demonstrates basic usage of the wrapper
- `drift_analysis.py`: A more advanced example that analyzes and visualizes drift in the crossbar
- `xbar_simulator_demo.ipynb`: A Jupyter notebook that provides an interactive demonstration

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

## Constants

- `voltage_pulse_height`: The height of the voltage pulse applied to the crossbar
- `simulation_time_step`: The time step used in the simulation

## License

This project is licensed under the same license as the original crossbar simulator. 