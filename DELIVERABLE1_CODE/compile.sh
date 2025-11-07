#!/bin/bash
# let's see where this goes...

# Configuration
MATRIX="1138_bus/1138_bus.mtx"
RUNS=12
THREADS=(1 2 4 8 16)  # adjust for your system
EXECUTABLES=("spmv_static" "spmv_dynamic" "spmv_guided")
CHUNK_SIZE=(10 100 1000)

# others to add: "spmv_collapse" "spmv_runtime" "spmv_auto" "spmv_chunked"!!!!!!!!!!!

# Output file
OUTPUT_FILE="results.csv"
echo "Executable,Threads,Run,ChunkSize,Time" > "$OUTPUT_FILE"

# Loop for each executable
for EXEC in "${EXECUTABLES[@]}"; do
    for T in "${THREADS[@]}"; do
        echo "Running $EXEC with $T threads..."
        export OMP_NUM_THREADS=$T

        for ((i=1; i<=RUNS; i++)); do
            if [ "$EXEC" == "spmv_chunked" ]; then
                for CHUNK in "${CHUNK_SIZE[@]}"; do
                    output=$(./"$EXEC" "$MATRIX" "$CHUNK" | grep "Elapsed time" | awk '{print $3}')
                    echo "$EXEC,$T,$i,$CHUNK,$output" >> "$OUTPUT_FILE"
                done
            else
                output=$(./"$EXEC" "$MATRIX" | grep "Elapsed time" | awk '{print $3}')
                echo "$EXEC,$T,$i,NA,$output" >> "$OUTPUT_FILE"
            fi
        done
    done
done
