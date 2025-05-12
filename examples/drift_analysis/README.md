# Drift Analysis Example

This example demonstrates how to analyze state drift in memristors over multiple inferences.

## Steps

1. Generate input files:
```bash
python ../../tools/create_input.py --input_rows 1 --weight_rows 32 --weight_cols 32 --output_dir ../../data
```

2. Run the simulator with repeated inferences:
```bash
../../bin/xbar_simulator ../../data
```

3. Analyze the drift:
```bash
python ../../tools/analyze_drift.py --input_dir ../../data
```

## Expected Results

The analysis will show how memristor states drift over multiple inferences, which can affect the reliability of the crossbar array over time.
