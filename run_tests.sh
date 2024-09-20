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

# Run the main program
./JCP_5_1_1

# Check if the Results folder and dope.log file exist
if [ -d "Results" ] && [ -f "Results/dope.log" ]; then
    echo "Output file 'Results/dope.log' exists."
else
    echo "Output file 'Results/dope.log' is missing."
    exit 1
fi

# Compare the generated log file with the expected log file
if diff -q Results/dope.log ../dope_Aug_12_2024.log; then
    echo "Output matches the expected results."
else
    echo "Output does not match the expected results."
    exit 1
fi

echo "All tests passed!"
