#!/bin/bash

# Run the main program
./JCP_5_1_1

# Check if the output folder exists
if [ -d "Results" ]; then
    echo "Output folder 'Results/' exists."
else
    echo "Output folder 'Results/' is missing."
    exit 1
fi

# Optionally, compare outputs with expected values (simplified example)
if diff -q Results/output.txt expected_results/output.txt; then
    echo "Output matches the expected results."
else
    echo "Output does not match the expected results."
    exit 1
fi

echo "All tests passed!"
