#!/bin/bash
# first attempt of the shell that should automize the different pragma loops...let's hope for the best

# configuration
MATRIX="1138_bus/1138_bus.mtx"
RUNS=12
THREADS=(1 2 4 8 16) # change according to your system
EXECUTABLES=("spmv_static" "spmv_dynamic" "spmv_collapse" "spmv_guided" "spmv_runtime" "spmv_auto" "spmv_chunked") # add others if necessary
CHUNK_SIZE=(10 100 1000)

# output file
OUTPUT_FILE="results.csv"
echo "Executable,Threads,Run,Time" > $OUTPUT_FILE 

# loop for each executable
for EXEC in ${EXECUTABLES[@]}; do
    # loop for each thread count
    for T in "${THREADS[@]}"; do
        echo "Running $EXEC with $T threads..."
        export OMP_NUM_THREADS=$T
        # loop for number of runs
        for ((i=1; i<=RUNS; i++)); do
            output=$(./$EXEC $MATRIX | grep "Elapsed time" | awk '{print $3}')
            echo "$EXEC,$T,$i,$output" >> $OUTPUT_FILE
            if [ "$EXEC" == "spmv_chunked" ]; then
                for CHUNK in "${CHUNK_SIZE[@]}"; do
                    TIME=$("./$EXEC" "$INPUT_FILE" "$CHUNK")
                    echo "$EXEC,$T,$i,$TIME" >> $OUTPUT_FILE
                done
            else
                TIME=$("./$EXEC" "$INPUT_FILE")
                echo "$EXEC,$T,$i,$TIME" >> $OUTPUT_FILE
            fi
        done
    done
done    