# Basic Inference Example

This example demonstrates how to perform a basic inference using the crossbar simulator.

## Steps

1. Generate input files:
```bash
python ../../tools/create_input.py --input_rows 10 --weight_rows 32 --weight_cols 32 --output_dir ../../data
```

2. Run the simulator:
```bash
../../bin/xbar_simulator ../../data
```

3. Analyze the results:
```bash
python ../../tools/analyze_drift.py --input_dir ../../data
```
