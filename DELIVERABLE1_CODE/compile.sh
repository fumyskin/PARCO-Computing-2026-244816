#!/bin/bash
# HPC benchmarking script for SpMV variants

# ===========================
# Configuration
# ===========================
MATRICES=("1138_bus/1138_bus.mtx" "utm5940/utm5940.mtx")
RUNS=12
THREADS=(1 2 4 8 16 32 64 96)
EXECUTABLES=("spmv_sequential" "spmv_static" "spmv_manual" "spmv_dynamic" "spmv_guided" \
             "spmv_runtime_static" "spmv_runtime_dynamic" "spmv_runtime_guided")
CHUNK_SIZES=(10 100 1000 10000)

RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

# ===========================
# Main benchmarking loop
# ===========================
for MATRIX_PATH in "${MATRICES[@]}"; do
    MATRIX_NAME=$(basename "$MATRIX_PATH" .mtx)
    echo "============================================"
    echo "Testing matrix: $MATRIX_NAME"
    echo "============================================"

    # Prepare single chunked results file
    CHUNKED_OUTPUT="${RESULTS_DIR}/${MATRIX_NAME}_chunked.csv"
    echo "Executable,Threads,Run,ChunkSize,Time" > "$CHUNKED_OUTPUT"

    # Loop for each executable
    for EXEC in "${EXECUTABLES[@]}"; do
        echo "Running executable: $EXEC"

        if [[ "$EXEC" == spmv_runtime_static || "$EXEC" == spmv_runtime_dynamic || "$EXEC" == spmv_runtime_guided ]]; then
            
            # Determine schedule type
            if [[ "$EXEC" == spmv_runtime_static ]]; then
                SCHED_TYPE="static"
            elif [[ "$EXEC" == spmv_runtime_dynamic ]]; then
                SCHED_TYPE="dynamic"
            else
                SCHED_TYPE="guided"
            fi

            # Loop over chunk sizes
            for CHUNK in "${CHUNK_SIZES[@]}"; do
                export OMP_SCHEDULE="${SCHED_TYPE},${CHUNK}"
                echo "  Schedule: ${SCHED_TYPE}, Chunk: ${CHUNK}"

                # Loop over thread counts
                for T in "${THREADS[@]}"; do
                    export OMP_NUM_THREADS=$T
                    echo "    Threads: $T"

                    for ((i=1; i<=RUNS; i++)); do
                        echo "      Run $i"
                        output=$(./"$EXEC" "$MATRIX_PATH" | grep "Elapsed time" | awk '{print $3}')
                        echo "$EXEC,$T,$i,$CHUNK,$output" >> "$CHUNKED_OUTPUT"
                    done
                done
            done

        else
            # Default case for other executables (non-runtime)
            OUTPUT_FILE="${RESULTS_DIR}/${MATRIX_NAME}_${EXEC}.csv"
            echo "Executable,Threads,Run,ChunkSize,Time" > "$OUTPUT_FILE"

            for T in "${THREADS[@]}"; do
                export OMP_NUM_THREADS=$T
                echo "  Threads: $T"

                for ((i=1; i<=RUNS; i++)); do
                    echo "      Run $i"
                    output=$(./"$EXEC" "$MATRIX_PATH" | grep "Elapsed time" | awk '{print $3}')
                    echo "$EXEC,$T,$i,NA,$output" >> "$OUTPUT_FILE"
                done
            done
        fi
    done

    echo "Finished all executables for matrix $MATRIX_NAME"
done

echo "All tests completed. Results saved in '$RESULTS_DIR/'"
