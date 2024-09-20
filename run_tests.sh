#!/bin/bash

# Create a log file for the terminal output
terminal_output_log="terminal_output.log"

# Run the main program and redirect both stdout and stderr to the log file
./JCP_5_1_1 > $terminal_output_log 2>&1

# Check if the Results folder and dope.log file exist
if [ -d "Results" ] && [ -f "Results/dope.log" ]; then
    echo "Output file 'Results/dope.log' exists." | tee -a $terminal_output_log
else
    echo "Output file 'Results/dope.log' is missing." | tee -a $terminal_output_log
    exit 1
fi

echo "Test execution completed. Logs saved in terminal_output.log"
