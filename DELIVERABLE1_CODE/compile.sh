#!/bin/bash
# let's see where this goes...

# Configuration
MATRICES=("1138_bus/1138_bus.mtx" "utm5940/utm5940.mtx")
RUNS=12
THREADS=(1 2 4 8 16 32 64)  # adjust for your system
EXECUTABLES=("spmv_sequential" "spmv_static" "spmv_manual")
CHUNK_SIZES=(10 100 1000)

# Directory for results
RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

# Loop for each matrix
for MATRIX_PATH in "${MATRICES[@]}"; do
    MATRIX_NAME=$(basename "$MATRIX_PATH" .mtx)
    OUTPUT_FILE="${RESULTS_DIR}/${MATRIX_NAME}.csv"

    echo "============================================"
    echo "Testing matrix: $MATRIX_NAME"
    echo "============================================"

    # Initialize output file for this matrix
    echo "Executable,Threads,Run,ChunkSize,Time" > "$OUTPUT_FILE"

    # Loop for each executable
    for EXEC in "${EXECUTABLES[@]}"; do
        echo "Running executable: $EXEC"

        # Loop for each thread count
        for T in "${THREADS[@]}"; do
            export OMP_NUM_THREADS=$T
            echo "  Threads: $T"

            # Run multiple times
            for ((i=1; i<=RUNS; i++)); do
                if [ "$EXEC" == "spmv_chunked" ]; then
                    for CHUNK in "${CHUNK_SIZES[@]}"; do
                        echo "    Run $i | Chunk $CHUNK"
                        output=$(./"$EXEC" "$MATRIX_PATH" "$CHUNK" | grep "Elapsed time" | awk '{print $3}')
                        echo "$EXEC,$T,$i,$CHUNK,$output" >> "$OUTPUT_FILE"
                    done
                else
                    echo "    Run $i"
                    output=$(./"$EXEC" "$MATRIX_PATH" | grep "Elapsed time" | awk '{print $3}')
                    echo "$EXEC,$T,$i,NA,$output" >> "$OUTPUT_FILE"
                fi
            done
        done
    done

    echo "Results for $MATRIX_NAME saved to $OUTPUT_FILE"
done

echo "All tests completed. Results saved in '$RESULTS_DIR/'"
