# Memristor Crossbar Simulator

A simulator for memristor crossbar arrays with non-linear behavior for neuromorphic computing applications.

## Overview

This simulator implements the JART VCM model for memristors and includes features like:
- Non-linear I-V characteristics
- Temperature effects
- Parasitic resistances
- State drift analysis

## Repository Structure

```
clean_xbar_model/
├── bin/                    # Compiled binaries
├── include/                # Header files
│   ├── core/               # Core simulator headers
│   ├── crossbar_model/     # Crossbar model headers
│   └── memristor_model/    # Memristor model headers
├── src/                    # Source files
│   ├── core/               # Core simulator source
│   ├── crossbar_model/     # Crossbar model implementation
│   └── memristor_model/    # Memristor model implementation
├── tools/                  # Utility scripts
│   ├── create_input.py     # Input file generator
│   └── analyze_drift.py    # Drift analysis tool
├── examples/               # Example configurations
├── data/                   # Input/output data directory
└── build/                  # Build artifacts
    └── obj/                # Object files
```

## Building the Simulator

```bash
make clean
make
```

## Usage

1. Generate input files:
```bash
python tools/create_input.py --input_rows 10 --weight_rows 32 --weight_cols 32 --output_dir data
```

2. Run the simulator:
```bash
./bin/xbar_simulator data
```

3. Analyze drift (for repeated inferences):
```bash
python tools/analyze_drift.py --input_dir data
```

## License

[License information] 