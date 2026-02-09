#!/bin/bash

# Configuration
SRC_DIR="../src"
INCLUDE_DIR="../include"
EXEC_PARALLEL="../results/parallel_io_reading/spmv_parallel_io.out"
MATRICES_DIR="../matrices"
RESULTS_DIR="../results/parallel_io_reading"
SUMMARY_DIR="../results/parallel_io_reading/summary"

# MPI Configuration
MAX_PROCESSES=128
PROCESSES=(1 2 4 8 16 32 64 128)

# I/O types to test
IO_TYPES=("sequential" "parallel")

# Test parameters
REPEATS=10
ROWS_PER_PROC=10000
NNZ_PER_ROW=50

# Compiler
MPICC="mpicc"
CFLAGS="-ffp-contract=fast -O3 -Wall -I$INCLUDE_DIR"

# Output files
IO_COMPARISON_CSV="$RESULTS_DIR/io_comparison.csv"
STRONG_SCALING_SEQ_CSV="$RESULTS_DIR/strong_scaling_sequential.csv"
STRONG_SCALING_PAR_CSV="$RESULTS_DIR/strong_scaling_parallel.csv"
WEAK_SCALING_SEQ_CSV="$RESULTS_DIR/weak_scaling_sequential.csv"
WEAK_SCALING_PAR_CSV="$RESULTS_DIR/weak_scaling_parallel.csv"

# Find all .mtx files
find_matrices() {
    MATRICES=($(ls "$MATRICES_DIR"/*.mtx 2>/dev/null | xargs -n 1 basename))
    if [ ${#MATRICES[@]} -eq 0 ]; then
        echo "ERROR: No .mtx files found in $MATRICES_DIR"
        exit 1
    fi
    echo "Found ${#MATRICES[@]} matrices: ${MATRICES[@]}"
}

# Compile the parallel I/O comparison executable
compile_code() {
    rm -f "$EXEC_PARALLEL"
    mkdir -p "$RESULTS_DIR"

    echo "Compiling parallel I/O comparison version..."
    
    # Include main_parallel.c and bonus_utils.c
    $MPICC $CFLAGS \
        "$SRC_DIR/main_parallel.c" \
        "$SRC_DIR/parallel_reading.c" \
        "$SRC_DIR/distributions.c" \
        "$SRC_DIR/communication.c" \
        "$SRC_DIR/mmio.c" \
        "$SRC_DIR/specifications.c" \
        -o "$EXEC_PARALLEL" -lm
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Compilation failed!"
        exit 1
    fi
    echo "Compilation successful"
}

# Extract summary paragraph from log
extract_summary_paragraph() {
    local log_file="$1"
    local summary_file="$2"
    local io_type="$3"
    local num_procs="$4"
    
    echo "========== [$io_type I/O] $num_procs processes ==========" >> "$summary_file"
    
    awk '
        /^========================================$/ {
            if (getline > 0 && $0 ~ /COMPREHENSIVE METRICS SUMMARY/) {
                print "========================================"
                print $0
                getline
                print $0
                while (getline > 0 && !/^========================================$/) {
                    print $0
                }
                if (/^========================================$/) print $0
            }
        }
    ' "$log_file" >> "$summary_file"
    
    echo "" >> "$summary_file"
}

# Initialize matrix summary file
initialize_matrix_summary() {
    local matrix_name="$1"
    local summary_file="$SUMMARY_DIR/${matrix_name%.mtx}_io_summary.txt"
    
    > "$summary_file"
}

# Run strong scaling test for a specific I/O type
run_strong_scaling() {
    local matrix="$1"
    local num_procs="$2"
    local io_type="$3"
    local csv_file="$4"

    local log_file="$RESULTS_DIR/logs/strong_${io_type}_${matrix%.mtx}_np${num_procs}.log"
    local summary_file="$SUMMARY_DIR/${matrix%.mtx}_io_summary.txt"
    mkdir -p "$RESULTS_DIR/logs"

    echo ""
    echo "STRONG SCALING [$io_type I/O]: $matrix | $num_procs processes"
    echo ""
    
    mpirun -np "$num_procs" "$EXEC_PARALLEL" "$MATRICES_DIR/$matrix" "$REPEATS" "$io_type" 2>&1 | \
        tee "$log_file" | \
        grep "^\[RESULT\]" | sed 's/\[RESULT\] //' | \
        awk -v matrix="${matrix%.mtx}" '{print $0","matrix}' >> "$csv_file"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: MPI execution failed"
        return 1
    fi

    extract_summary_paragraph "$log_file" "$summary_file" "$io_type" "$num_procs"

    echo " -> Data appended to: $csv_file"
    echo " -> Summary appended to: $summary_file"
    echo " -> Full log saved to: $log_file"
}

# Run weak scaling test for a specific I/O type
run_weak_scaling() {
    local num_procs="$1"
    local io_type="$2"
    local csv_file="$3"
    local total_rows=$((ROWS_PER_PROC * num_procs))

    local log_file="$RESULTS_DIR/logs/weak_${io_type}_np${num_procs}.log"
    mkdir -p "$RESULTS_DIR/logs"

    echo ""
    echo "WEAK SCALING [$io_type I/O]: $num_procs processes | ${total_rows} total rows"
    echo ""
    
    mpirun -np "$num_procs" "$EXEC_PARALLEL" dummy "$REPEATS" "$io_type" "$ROWS_PER_PROC" "$NNZ_PER_ROW" 2>&1 | \
        tee "$log_file" | \
        grep "^\[RESULT\]" | sed 's/\[RESULT\] //' >> "$csv_file"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: MPI execution failed"
        return 1
    fi

    echo " -> Data appended to: $csv_file"
    echo " -> Full log saved to: $log_file"
}

# Extract I/O time from log and append to comparison CSV
extract_io_time() {
    local log_file="$1"
    local matrix="$2"
    local num_procs="$3"
    local io_type="$4"
    
    # Extract I/O time from log
    local io_time=$(grep "I/O TIME:" "$log_file" | awk '{print $3}')
    
    if [ -n "$io_time" ]; then
        echo "${matrix%.mtx},$num_procs,$io_type,$io_time" >> "$IO_COMPARISON_CSV"
    fi
}

# Initialize CSV files
initialize_csv_files() {
    echo "rank,num_procs,io_type,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,io_time,matrix" > "$STRONG_SCALING_SEQ_CSV"
    echo "rank,num_procs,io_type,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,io_time,matrix" > "$STRONG_SCALING_PAR_CSV"
    echo "rank,num_procs,io_type,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,io_time" > "$WEAK_SCALING_SEQ_CSV"
    echo "rank,num_procs,io_type,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,io_time" > "$WEAK_SCALING_PAR_CSV"
    echo "matrix,num_procs,io_type,io_time_seconds" > "$IO_COMPARISON_CSV"
    echo "CSV files initialized"
}

# Main execution

echo ""
echo "========================================"
echo "PARALLEL I/O BENCHMARK SUITE"
echo "========================================"
echo ""

# Cleanup and setup
rm -rf "$RESULTS_DIR"/logs
rm -rf "$SUMMARY_DIR"
rm -f "$IO_COMPARISON_CSV" "$STRONG_SCALING_SEQ_CSV" "$STRONG_SCALING_PAR_CSV"
rm -f "$WEAK_SCALING_SEQ_CSV" "$WEAK_SCALING_PAR_CSV"
mkdir -p "$RESULTS_DIR"
mkdir -p "$SUMMARY_DIR"

find_matrices
compile_code
initialize_csv_files

# Initialize summary files for each matrix
for matrix in "${MATRICES[@]}"; do
    initialize_matrix_summary "$matrix"
done

# ============================================================
# STRONG SCALING TESTS - SEQUENTIAL I/O
# ============================================================
echo ""
echo "========================================"
echo "Strong Scaling Tests - SEQUENTIAL I/O"
echo "========================================"
echo ""

for matrix in "${MATRICES[@]}"; do
    echo ""
    echo "Matrix: $matrix"
    for num_procs in "${PROCESSES[@]}"; do
        if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
            run_strong_scaling "$matrix" "$num_procs" "sequential" "$STRONG_SCALING_SEQ_CSV"
            
            # Extract I/O time for comparison
            log_file="$RESULTS_DIR/logs/strong_sequential_${matrix%.mtx}_np${num_procs}.log"
            extract_io_time "$log_file" "$matrix" "$num_procs" "sequential"
            
            sleep 1
        fi
    done
done

# ============================================================
# STRONG SCALING TESTS - PARALLEL I/O
# ============================================================
echo ""
echo "========================================"
echo "Strong Scaling Tests - PARALLEL I/O"
echo "========================================"
echo ""

for matrix in "${MATRICES[@]}"; do
    echo ""
    echo "Matrix: $matrix"
    for num_procs in "${PROCESSES[@]}"; do
        if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
            run_strong_scaling "$matrix" "$num_procs" "parallel" "$STRONG_SCALING_PAR_CSV"
            
            # Extract I/O time for comparison
            log_file="$RESULTS_DIR/logs/strong_parallel_${matrix%.mtx}_np${num_procs}.log"
            extract_io_time "$log_file" "$matrix" "$num_procs" "parallel"
            
            sleep 1
        fi
    done
done

# ============================================================
# WEAK SCALING TESTS - SEQUENTIAL I/O
# ============================================================
echo ""
echo "========================================"
echo "Weak Scaling Tests - SEQUENTIAL I/O"
echo "========================================"
echo ""

for num_procs in "${PROCESSES[@]}"; do
    if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
        run_weak_scaling "$num_procs" "sequential" "$WEAK_SCALING_SEQ_CSV"
        sleep 1
    fi
done

# ============================================================
# WEAK SCALING TESTS - PARALLEL I/O
# ============================================================
echo ""
echo "========================================"
echo "Weak Scaling Tests - PARALLEL I/O"
echo "========================================"
echo ""

for num_procs in "${PROCESSES[@]}"; do
    if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
        run_weak_scaling "$num_procs" "parallel" "$WEAK_SCALING_PAR_CSV"
        sleep 1
    fi
done

# ============================================================
# GENERATE I/O COMPARISON SUMMARY
# ============================================================
echo ""
echo "========================================"
echo "I/O TIME COMPARISON ANALYSIS"
echo "========================================"
echo ""

# Create a summary report
COMPARISON_REPORT="$RESULTS_DIR/io_speedup_report.txt"
> "$COMPARISON_REPORT"

echo "I/O SPEEDUP ANALYSIS" >> "$COMPARISON_REPORT"
echo "====================" >> "$COMPARISON_REPORT"
echo "" >> "$COMPARISON_REPORT"

for matrix in "${MATRICES[@]}"; do
    echo "Matrix: ${matrix%.mtx}" >> "$COMPARISON_REPORT"
    echo "----------------------------------------" >> "$COMPARISON_REPORT"
    printf "%-12s %-15s %-15s %-10s\n" "Processes" "Sequential(s)" "Parallel(s)" "Speedup" >> "$COMPARISON_REPORT"
    
    for num_procs in "${PROCESSES[@]}"; do
        if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
            # Get times from CSV
            seq_time=$(grep "^${matrix%.mtx},$num_procs,sequential," "$IO_COMPARISON_CSV" | cut -d',' -f4)
            par_time=$(grep "^${matrix%.mtx},$num_procs,parallel," "$IO_COMPARISON_CSV" | cut -d',' -f4)
            
            if [ -n "$seq_time" ] && [ -n "$par_time" ]; then
                speedup=$(echo "scale=2; $seq_time / $par_time" | bc)
                printf "%-12s %-15s %-15s %-10s\n" "$num_procs" "$seq_time" "$par_time" "${speedup}x" >> "$COMPARISON_REPORT"
            fi
        fi
    done
    echo "" >> "$COMPARISON_REPORT"
done

cat "$COMPARISON_REPORT"

# ============================================================
# FINAL SUMMARY
# ============================================================
echo ""
echo "========================================"
echo "BENCHMARK SUMMARY"
echo "========================================"
echo ""
echo "Results Files:"
echo " Sequential Strong Scaling: $STRONG_SCALING_SEQ_CSV"
echo " Parallel Strong Scaling:   $STRONG_SCALING_PAR_CSV"
echo " Sequential Weak Scaling:   $WEAK_SCALING_SEQ_CSV"
echo " Parallel Weak Scaling:     $WEAK_SCALING_PAR_CSV"
echo " I/O Time Comparison:       $IO_COMPARISON_CSV"
echo " I/O Speedup Report:        $COMPARISON_REPORT"
echo " Full Logs:                 $RESULTS_DIR/logs/"
echo " Matrix Summaries:          $SUMMARY_DIR/"
echo ""
echo "Data points collected:"
echo " Sequential Strong: $(tail -n +2 "$STRONG_SCALING_SEQ_CSV" | wc -l) measurements"
echo " Parallel Strong:   $(tail -n +2 "$STRONG_SCALING_PAR_CSV" | wc -l) measurements"
echo " Sequential Weak:   $(tail -n +2 "$WEAK_SCALING_SEQ_CSV" | wc -l) measurements"
echo " Parallel Weak:     $(tail -n +2 "$WEAK_SCALING_PAR_CSV" | wc -l) measurements"
echo " I/O Comparisons:   $(tail -n +2 "$IO_COMPARISON_CSV" | wc -l) timing records"
echo ""
echo "Matrix summary files created:"
for matrix in "${MATRICES[@]}"; do
    summary_file="$SUMMARY_DIR/${matrix%.mtx}_io_summary.txt"
    if [ -f "$summary_file" ]; then
        echo " - ${matrix%.mtx}_io_summary.txt"
    fi
done
echo ""
echo "I/O Types tested: ${IO_TYPES[@]}"
echo ""