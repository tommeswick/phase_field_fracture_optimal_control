#!/bin/bash

# Create a log file for the terminal output
terminal_output_log="terminal_output.log"

# Run the main program and redirect both stdout and stderr to the log file
echo "Running JCP_5_1_1 program..." | tee -a $terminal_output_log
./JCP_5_1_1 > $terminal_output_log 2>&1
status=$?
if [ $status -ne 0 ]; then
    echo "JCP_5_1_1 program failed with exit code $status." | tee -a $terminal_output_log
    exit 1
fi

# Debugging: Check if the Results directory exists
if [ -d "Results" ]; then
    echo "Results directory exists." | tee -a $terminal_output_log
else
    echo "Results directory is missing." | tee -a $terminal_output_log
    exit 1
fi

# Debugging: List the contents of the Results directory
echo "Listing contents of the Results directory:" | tee -a $terminal_output_log
ls -l Results | tee -a $terminal_output_log

# Check if the dope.log file exists
if [ -f "Results/dope.log" ]; then
    echo "Output file 'Results/dope.log' exists." | tee -a $terminal_output_log
else
    echo "Output file 'Results/dope.log' is missing." | tee -a $terminal_output_log
    exit 1
fi

echo "Test execution completed. Logs saved in terminal_output.log"
