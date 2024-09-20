/** Copyright (C) 2024 by Denis Khimin, Marc C. Steinbach, Thomas Wick   
 *
 *  This file is part of the phase_field_fracture_optimal_control library.
 *
 *  The phase_field_fracture_optimal_control library is free software; 
 *  you can use it, redistribute it, and/or modify it
 *  under the terms of the GNU Lesser General
 *  Public License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *  The full text of the license can be found in the file LICENSE. 
**/

#!/bin/bash

# Arguments: $1 = generated dope.log, $2 = reference log

generated_log="$1"
reference_log="$2"

# Check if both files exist
if [ ! -f "$generated_log" ]; then
    echo "Generated log file '$generated_log' does not exist."
    exit 1
fi

if [ ! -f "$reference_log" ]; then
    echo "Reference log file '$reference_log' does not exist."
    exit 1
fi

# Compare the contents of the two files, ignoring minor differences (up to 5% allowed difference)
similarity_threshold=0.05

# Count the number of lines in each file
num_lines_generated=$(wc -l < "$generated_log")
num_lines_expected=$(wc -l < "$reference_log")

# Calculate the allowed difference
allowed_difference=$(echo "$num_lines_expected * $similarity_threshold" | bc)
difference=$(echo "$num_lines_expected - $num_lines_generated" | bc)

# Compare the generated log file with the expected log file
if diff -q generated_log reference_log; then
    echo "Output matches the expected results."
else
    echo "Output does not match the expected results."
    exit 1
fi

#if (( $(echo "$difference > $allowed_difference" | bc -l) || $(echo "$difference < -$allowed_difference" | bc -l) )); then
#    echo "Output differs from the expected results by more than 5%."
#    exit 1
#else
#    echo "Output is within the 5% tolerance range compared to the expected results."
#fi

echo "Comparison completed successfully."
