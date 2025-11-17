#!/bin/bash
# let's see where this goes...

# Configuration
# define exactly, for each NUMA block, the cores used 
# based on specifications found in experiment1.c
# export OMP_NUM_THREADS=96   # set manually number of threads
export OMP_PROC_BIND=spread # evenly spread threads on the NUMA cores, still PINNING them and assuring thread affinity 
export OMP_PLACES="\
{0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92},\
{1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93},\
{2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94},\
{3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79,83,87,91,95}"


MATRICES=("flowmeter0/flowmeter0.mtx" "1138_bus/1138_bus.mtx")
RUNS=10
THREADS=(1 2 4 8 16 32 64 96)
EXECUTABLES=("spmv_exp1")
CHUNK_SIZES=(10 100 1000 10000)

RESULTS_DIR="DELIVERABLE1_RES/ARCH5"
rm -rf "$RESULTS_DIR"        # Remove old results completely
mkdir -p "$RESULTS_DIR"      # Create fresh folder

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

    # Default case for other executables (non-runtime)
    OUTPUT_FILE="${RESULTS_DIR}/${MATRIX_NAME}.csv"
    echo "Executable,Threads,Run,ChunkSize,Time" > "$OUTPUT_FILE"


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
