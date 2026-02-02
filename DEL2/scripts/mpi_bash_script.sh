#!/bin/bash

# configuration
SRC_DIR="../src"
INCLUDE_DIR="../include"
EXEC_MPI="../results/spmv_mpi.out"
EXEC_HYBRID="../results/spmv_hybrid.out"
MATRICES_DIR="../matrices"
RESULTS_DIR="../results"

# MPI Configuration
PROCESSES_PER_NODE=24
MAX_PROCESSES=128
PROCESSES=(1 2 4 8 16 32 64 128)

# distribution strategies
DISTRIBUTIONS=("1D" "2D" "2D_cyclic")

# OpenMP 
OMP_SCHEDULE="dynamic,64"
OMP_PROC_BIND="close"

# test parameters
REPEATS=10
ROWS_PER_PROC=10000
NNZ_PER_ROW=50

# compiler
MPICC="mpicc"
CFLAGS_MPI="-ffp-contract=fast -O3 -Wall -I$INCLUDE_DIR"
CFLAGS_HYBRID="-fopenmp -DHYBRID -ffp-contract=fast -O3 -Wall -I$INCLUDE_DIR" 

# output files
STRONG_SCALING_CSV="$RESULTS_DIR/strong_scaling_all.csv"
WEAK_SCALING_CSV="$RESULTS_DIR/weak_scaling_all.csv"
STRONG_SCALING_HYBRID_CSV="$RESULTS_DIR/strong_scaling_hybrid.csv"
WEAK_SCALING_HYBRID_CSV="$RESULTS_DIR/weak_scaling_hybrid.csv"

find_matrices() {
    MATRICES=($(ls "$MATRICES_DIR"/*.mtx 2>/dev/null | xargs -n 1 basename))
    if [ ${#MATRICES[@]} -eq 0 ]; then
        echo "ERROR: No .mtx files found in $MATRICES_DIR"
        exit 1
    fi
    echo "Found ${#MATRICES[@]} matrices: ${MATRICES[@]}"
}

compile_code() {
    rm -f "$EXEC_MPI" "$EXEC_HYBRID"
    mkdir -p "$RESULTS_DIR"

    # Compile MPI version
    echo "Compiling MPI-only version..."
    $MPICC $CFLAGS_MPI "$SRC_DIR"/*.c -o "$EXEC_MPI" -lm
    if [ $? -ne 0 ]; then
        echo "ERROR: MPI compilation failed!"
        exit 1
    fi
    echo "MPI compilation successful"

    #Compile MPI+OpenMP version TO DO LATER
    echo "Compiling MPI+OpenMP hybrid version..."
    $MPICC $CFLAGS_HYBRID "$SRC_DIR"/*.c -o "$EXEC_HYBRID" -lm
    if [ $? -ne 0 ]; then
        echo "ERROR: Hybrid compilation failed!"
        exit 1
    fi
    echo "Hybrid compilation successful"
}

run_strong_scaling() {
    local matrix="$1"
    local num_procs="$2"
    local exec="$3"
    local csv_file="$4"
    local mode="$5"
    local dist_type="$6"

    local log_file="$RESULTS_DIR/logs/strong_${mode}_${dist_type}_${matrix%.mtx}_np${num_procs}.log"
    mkdir -p "$RESULTS_DIR/logs"

    if [ "$mode" = "HYBRID" ]; then
        local threads=$((PROCESSES_PER_NODE / num_procs))
        if [ $threads -lt 1 ]; then
            threads=1
        fi
        
        echo ""
        echo "STRONG SCALING [$mode/$dist_type]: $matrix | $num_procs processes | $threads threads/process"
        echo ""
        
        OMP_NUM_THREADS=$threads \
        OMP_SCHEDULE="$OMP_SCHEDULE" \
        OMP_PROC_BIND="$OMP_PROC_BIND" \
        mpirun -np "$num_procs" \
        "$exec" "$MATRICES_DIR/$matrix" "$REPEATS" "$dist_type" 2>&1 | tee "$log_file" | \
            grep "^\[RESULT\]" | sed 's/\[RESULT\] //' | \
            awk -v matrix="${matrix%.mtx}" '{print $0","matrix}' >> "$csv_file"
    else
        echo ""
        echo "STRONG SCALING [$mode/$dist_type]: $matrix | $num_procs processes"
        echo ""
        
        mpirun -np "$num_procs" "$exec" "$MATRICES_DIR/$matrix" "$REPEATS" "$dist_type" 2>&1 | tee "$log_file" | \
            grep "^\[RESULT\]" | sed 's/\[RESULT\] //' | \
            awk -v matrix="${matrix%.mtx}" '{print $0","matrix}' >> "$csv_file"
    fi

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: MPI execution failed"
        return 1
    fi

    echo " -> Data appended to: $csv_file"
    echo " -> Full log saved to: $log_file"
}

run_weak_scaling() {
    local num_procs="$1"
    local exec="$2"
    local csv_file="$3"
    local mode="$4"
    local dist_type="$5"
    local total_rows=$((ROWS_PER_PROC * num_procs))

    local log_file="$RESULTS_DIR/logs/weak_${mode}_${dist_type}_np${num_procs}.log"
    mkdir -p "$RESULTS_DIR/logs"

    if [ "$mode" = "HYBRID" ]; then
        local threads=$((PROCESSES_PER_NODE / num_procs))
        if [ $threads -lt 1 ]; then
            threads=1
        fi
        
        echo ""
        echo "WEAK SCALING [$mode/$dist_type]: $num_procs processes | ${total_rows} total rows | $threads threads/process"
        echo ""
        
        OMP_NUM_THREADS=$threads \
        OMP_SCHEDULE="$OMP_SCHEDULE" \
        OMP_PROC_BIND="$OMP_PROC_BIND" \
        mpirun -np "$num_procs" \
               "$exec" dummy "$REPEATS" "$dist_type" "$ROWS_PER_PROC" "$NNZ_PER_ROW" 2>&1 | \
            tee "$log_file" | grep "^\[RESULT\]" | sed 's/\[RESULT\] //' >> "$csv_file"
    else
        echo ""
        echo "WEAK SCALING [$mode]: $num_procs processes | ${total_rows} total rows"
        echo ""
        
        mpirun -np "$num_procs" "$exec" dummy "$REPEATS" "$dist_type" "$ROWS_PER_PROC" "$NNZ_PER_ROW" 2>&1 | \
            tee "$log_file" | grep "^\[RESULT\]" | sed 's/\[RESULT\] //' >> "$csv_file"
    fi

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: MPI execution failed"
        return 1
    fi

    echo " -> Data appended to: $csv_file"
    echo " -> Full log saved to: $log_file"
}

initialize_csv_files() {
    # Headers for all CSV files
    echo "rank,num_procs,distribution,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,matrix" > "$STRONG_SCALING_CSV"
    echo "rank,num_procs,distribution,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops" > "$WEAK_SCALING_CSV"
    echo "rank,num_procs,distribution,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops,matrix" > "$STRONG_SCALING_HYBRID_CSV"
    echo "rank,num_procs,distribution,run,elapsed_time,comm_time,local_nz,ghost_entries,local_flops" > "$WEAK_SCALING_HYBRID_CSV"
    echo "CSV files initialized"
}

# main 

rm -rf "$RESULTS_DIR"/logs
rm -f "$STRONG_SCALING_CSV" "$WEAK_SCALING_CSV" "$STRONG_SCALING_HYBRID_CSV" "$WEAK_SCALING_HYBRID_CSV"
mkdir -p "$RESULTS_DIR"

find_matrices
compile_code
initialize_csv_files

# Strong Scaling - MPI (test all distribution types)
echo ""
echo "========================================"
echo "Strong Scaling Tests (MPI-only)"
echo "========================================"
echo ""

for dist_type in "${DISTRIBUTIONS[@]}"; do
    echo ""
    echo "Testing distribution: $dist_type"
    echo "----------------------------------------"
    
    for matrix in "${MATRICES[@]}"; do
        echo ""
        echo "Matrix: $matrix"
        for num_procs in "${PROCESSES[@]}"; do
            if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
                # Skip 2D distributions for 1 process (not meaningful)
                if [ "$num_procs" -eq 1 ] && [[ "$dist_type" == "2D"* ]]; then
                    echo "Skipping $dist_type for 1 process (using 1D instead)"
                    continue
                fi
                
                run_strong_scaling "$matrix" "$num_procs" "$EXEC_MPI" "$STRONG_SCALING_CSV" "MPI" "$dist_type"
                sleep 1
            fi
        done
    done
done

# Weak Scaling - MPI (test all distribution types)
echo ""
echo "========================================"
echo "Weak Scaling Tests (MPI-only)"
echo "========================================"
echo ""

for dist_type in "${DISTRIBUTIONS[@]}"; do
    echo ""
    echo "Testing distribution: $dist_type"
    echo "----------------------------------------"
    
    for num_procs in "${PROCESSES[@]}"; do
        if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
            # Skip 2D distributions for 1 process
            if [ "$num_procs" -eq 1 ] && [[ "$dist_type" == "2D"* ]]; then
                echo "Skipping $dist_type for 1 process"
                continue
            fi
            
            run_weak_scaling "$num_procs" "$EXEC_MPI" "$WEAK_SCALING_CSV" "MPI" "$dist_type"
            sleep 1
        fi
    done
done


echo ""
echo "========================================"
echo "Strong Scaling Tests (MPI+OpenMP)"
echo "========================================"
echo ""

for dist_type in "${DISTRIBUTIONS[@]}"; do
    echo ""
    echo "Testing distribution: $dist_type"
    echo "----------------------------------------"
    
    for matrix in "${MATRICES[@]}"; do
        echo ""
        echo "Matrix: $matrix"
        for num_procs in "${PROCESSES[@]}"; do
            if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
                # Skip 2D distributions for 1 process
                if [ "$num_procs" -eq 1 ] && [[ "$dist_type" == "2D"* ]]; then
                    echo "Skipping $dist_type for 1 process (using 1D instead)"
                    continue
                fi
                
                run_strong_scaling "$matrix" "$num_procs" "$EXEC_HYBRID" "$STRONG_SCALING_HYBRID_CSV" "HYBRID" "$dist_type"
                sleep 1
            fi
        done
    done
done

# Weak Scaling - MPI+OpenMP (test all distribution types)
echo ""
echo "========================================"
echo "Weak Scaling Tests (MPI+OpenMP)"
echo "========================================"
echo ""

for dist_type in "${DISTRIBUTIONS[@]}"; do
    echo ""
    echo "Testing distribution: $dist_type"
    echo "----------------------------------------"
    
    for num_procs in "${PROCESSES[@]}"; do
        if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
            # Skip 2D distributions for 1 process
            if [ "$num_procs" -eq 1 ] && [[ "$dist_type" == "2D"* ]]; then
                echo "Skipping $dist_type for 1 process"
                continue
            fi
            
            run_weak_scaling "$num_procs" "$EXEC_HYBRID" "$WEAK_SCALING_HYBRID_CSV" "HYBRID" "$dist_type"
            sleep 1
        fi
    done
done

echo ""
echo "========================================"
echo "BENCHMARK SUMMARY"
echo "========================================"
echo ""
echo "Results Files:"
echo " MPI Strong Scaling: $STRONG_SCALING_CSV"
echo " MPI Weak Scaling: $WEAK_SCALING_CSV"
echo " MPI+OpenMP Strong Scaling: $STRONG_SCALING_HYBRID_CSV"
echo " MPI+OpenMP Weak Scaling: $WEAK_SCALING_HYBRID_CSV"
echo " Full Logs: $RESULTS_DIR/logs/"
echo ""
echo "Data points collected:"
echo " MPI Strong: $(tail -n +2 "$STRONG_SCALING_CSV" | wc -l) measurements"
echo " MPI Weak: $(tail -n +2 "$WEAK_SCALING_CSV" | wc -l) measurements"
echo " MPI+OpenMP Strong: $(tail -n +2 "$STRONG_SCALING_HYBRID_CSV" | wc -l) measurements"
echo " MPI+OpenMP Weak: $(tail -n +2 "$WEAK_SCALING_HYBRID_CSV" | wc -l) measurements"
echo ""
echo "Distributions tested: ${DISTRIBUTIONS[@]}"
echo ""




# Strong Scaling - MPI
# echo ""
# echo "Strong Scaling Tests (MPI-only):"
# echo ""

# for matrix in "${MATRICES[@]}"; do
#     echo ""
#     echo "Testing matrix: $matrix"
#     echo "-----------------------------------------"
#     for num_procs in "${PROCESSES[@]}"; do
#         if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
#             run_strong_scaling "$matrix" "$num_procs" "$EXEC_MPI" "$STRONG_SCALING_CSV" "MPI"
#             sleep 1
#         fi
#     done
# done

# # Weak Scaling - MPI
# echo ""
# echo "Weak Scaling Tests (MPI-only):"
# echo ""

# for num_procs in "${PROCESSES[@]}"; do
#     if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
#         run_weak_scaling "$num_procs" "$EXEC_MPI" "$WEAK_SCALING_CSV" "MPI"
#         sleep 1
#     fi
# done

# # Strong Scaling - MPI+OpenMP
# echo ""
# echo "Strong Scaling Tests (MPI+OpenMP):"
# echo ""

# for matrix in "${MATRICES[@]}"; do
#     echo ""
#     echo "Testing matrix: $matrix"
#     echo "-----------------------------------------"
#     for num_procs in "${PROCESSES[@]}"; do
#         if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
#             run_strong_scaling "$matrix" "$num_procs" "$EXEC_HYBRID" "$STRONG_SCALING_HYBRID_CSV" "HYBRID"
#             sleep 1
#         fi
#     done
# done

# # Weak Scaling - MPI+OpenMP
# echo ""
# echo "Weak Scaling Tests (MPI+OpenMP):"
# echo ""

# for num_procs in "${PROCESSES[@]}"; do
#     if [ "$num_procs" -le "$MAX_PROCESSES" ]; then
#         run_weak_scaling "$num_procs" "$EXEC_HYBRID" "$WEAK_SCALING_HYBRID_CSV" "HYBRID"
#         sleep 1
#     fi
# done

# echo ""
# echo "BENCHMARK SUMMARY"
# echo ""
# echo "Results Files:"
# echo " MPI Strong Scaling: $STRONG_SCALING_CSV"
# echo " MPI Weak Scaling: $WEAK_SCALING_CSV"
# echo " MPI+OpenMP Strong Scaling: $STRONG_SCALING_HYBRID_CSV"
# echo " MPI+OpenMP Weak Scaling: $WEAK_SCALING_HYBRID_CSV"
# echo " Full Logs: $RESULTS_DIR/logs/"
# echo ""
# echo "Data points collected:"
# echo " MPI Strong: $(tail -n +2 "$STRONG_SCALING_CSV" | wc -l) measurements"
# echo " MPI Weak: $(tail -n +2 "$WEAK_SCALING_CSV" | wc -l) measurements"
# echo " MPI+OpenMP Strong: $(tail -n +2 "$STRONG_SCALING_HYBRID_CSV" | wc -l) measurements"
# echo " MPI+OpenMP Weak: $(tail -n +2 "$WEAK_SCALING_HYBRID_CSV" | wc -l) measurements"
# echo ""