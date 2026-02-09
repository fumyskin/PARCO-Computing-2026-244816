# DELIVERABLE 2 run instructions

## Overview
This folder contains the DELIVERABLE 2 source code.
The following code is a playground to test different MPI approaches to the classical SpMV kernel.
Included in this file are the following implementation:

- **BASELINE**: 1D cyclic distribution with halo exchange management
- parallel I/O reading approach to Matrix Market files 
- 2D block and 2D cyclic distributions
- Hybrid MPI approach

However, the DELIVERABLE2 content will focus ONLY on exploring, apart from the baseline, **parallel I/O reading** approach and **2D block partitioning**, studying them with **bonus metrics**.  


## FOLDER ORGANIZATION
The repository is organized as follows:
``` bash
.
├── Deliverable2-MPI.pdf
├── Parco_ADY_2526 _InstrutionReport.pdf
├── Parco_ADY_2526_ReportExample.pdf
├── README.md
├── include
│   ├── communication.h
│   ├── distributions.h
│   ├── mmio.h
│   ├── parallel_reading.h
│   ├── specifications.h
│   └── tests.h
├── matrices                            #input matrices
│   └── 1138_bus.mtx
├── plots                               #plots
│   ├── batch1                          #plots of the 3 smallest matrices
│   ├── batch2                          #plots of the 3 biggest matrices
│   └── parallel_io_reading             #plots for sequential vs parallel reading
├── results
│   ├── batch1                          #results of the 3 smallest matrices
│   ├── batch2                          #results of the 3 biggest matrices
│   └── parallel_io_reading             #results for sequential vs parallel reading
├── scripts
│   ├── mpi_bash_baseline.sh            #bash script for main_baseline.c
│   ├── mpi_bash_parallel.sh            #bash script for main_parallel.c
│   ├── mpi_job_script_baseline.pbs     #pbs script for mpi_bash_baseline.sh
│   ├── mpi_job_script_parallel.pbs     #pbs script for mpo_bash_parallel.sh
│   ├── plot_parallel.py                #python script for parallel_io_reading results
│   └── plot_script.py                  #python script for 1D_and_2D results
└── src
    ├── communication.c                 #utils for finding/exchanging ghost entries
    ├── distributions.c                 #utils for 1D, 2D (block and cyclic) distributions
    ├── main_baseline.c                 #main for 1D, 2D (block and cyclic), HYBRID MPI runs
    ├── main_parallel.c                 #main for sequential and parallel I/O reading comparison
    ├── mmio.c                          #utils for reading Matrix Market files
    ├── parallel_reading.c              #utils for parallel reading I/O
    ├── specifications.c                #utils for SpMV, COO->CSR operations
    └── tests.c                         #tests for 1D, 2D block and 2D cyclic
```

## HOW TO REPRODUCE
### MATRIX EXTRACTION
> [!WARNING]
> The full matrices used in the experiments are not included in this repository due to size constraints.

To reproduce the results:
- Download the following matrices from <https://sparse.tamu.edu/>:
    - Andrews.mtx
    - bcsstk27.mtx
    - msdoor.mtx
    - Hook_1498.mtx
    - tmt_sym.mtx
- Place the `.tar.gz` files in the `/matrices` folder
- extract them by using the following command:
``` bash
tar -xzf <name_of_matrix>.tar.gz
```
- Extract the .mtx file from their respective folders. Make sure that in the `/matrices` folder are just the `.mtx` files

## JOB SUBMISSION
To submit the .pbs that executes 1D and 2D block distributions on SpMV go into the `/scripts` folder and type:
```bash
qsub mpi_job_script_baseline.pbs
```

To submit the .pbs that executes sequential vs parallel I/O reading of Matrix Market file go into the `/scripts` folder and type:
```bash
qsub mpi_job_script_parallel.pbs
```

The results of the simulations are going to be collected into the `/results` folder :
- `/1D_and_2D` folder contains the results of the 1D and 2D distribution simulations (as well as 2D_cyclic and HYBRID, but we're not concerned with those): 
    - `1D_and_2D/logs` -> contains single logs
    - `1D_and_2D/summary` -> **SUMMARY, FOR EACH MATRIX ANALYIZED, OF ALL THE METRICS MEASURED FOR EACH ITERATION** 
    - `strong_scaling_all.csv` -> contains strong scaling information for all matrices measured
    - `weak_scaling_all.csv` -> contains weak scaling information for all matrices measured
    - `weak_scaling_hybrid.csv`, `weak_scaling_hybrid.csv` -> strong and weak scaling information of HYRBID MPI approach, but for the deliverable we are NOT concerned with these
- `/parallel_io_reading` folder contains information regarding sequential and parallel I/O reading
    - `parallel_io_reading/logs` -> contains single logs
    - `parallel_io_reading/summary` -> **SUMMARY, FOR EACH MATRIX ANALYIZED, OF ALL THE METRICS MEASURED FOR EACH ITERATION** 
    - `strong/weak_scaling_sequential.csv` 
    - `strong/weak_scaling_parallel.csv`

> [!WARNING]
> The mpi_job_script_baseline.pbs will run also the 2D_cyclic and HYBRID MPI experiments as well. Since the cluster HPC2 didn't work anymore I preferred to keep the original experimental setup as it is in order to be sure the source code is reproducible. 2D_cyclic and HYBRID MPI data are NOT to be considered

> [!IMPORTANT]
> The `../summary` folders contain the **TOTAL BENCHMARK INFORMATION** that is human readable. **PLEASE REFER TO THESE to check metrics** for a given matrix.


## PLOTTING
To plot the results for 1D and 2D distribution SpMV go into the `/scripts` folder and type:
``` bash
module load python-3.10.14_gcc91  #command to make sure python3 is installed
python3 plot_script.py
```
The contents plotted are going to be:
- `comm_volume_strong` -> heatmap + respective `.csv` of ghost entries per rank number
- `comp_vs_comm_strong` -> computation vs communication histograms for strong scaling
- `comparison` -> comparison for MPI strong scaling for different matrices
- `data_reduction` -> `.csv` files associated for each iteration for each matrix
- `speedup_strong` -> speedup for each matrix for strong scaling
- `weak_scaling` -> weak scaling info

**NB: the results labelled with _Hybrid are the hyrbid MPI results so should not be considered**

To compute plots for sequential and parallel I/O reading go into the `/scripts` folder and type:
``` bash
module load python-3.10.14_gcc91
python3 plot_parallel.py
```

The contents plotted are going to be speedup, io time comparison and computation vs communication breakdown

## GUIDELINES: HOW TO NAVIGATE RESULTS AND PLOTS
Regarding results, the MOST TELLING file I recommend checking for benchmark analysis is the `/results/summary` folder: for each process count are summarized:
- `Time statistics`: computation vs communcation time, total execution time
- `Load balance metrics`: min/avg/max nnz entries
- `Communication metrics`: ghost entries number
- `GFLOPS`

Regarding plots, I recommend checking `speedup_strong`, `comp_vs_comm_strong`, `weak_scaling` the most: they illustrate clearly the most important conclusions of each simulation.

> [!IMPORTANT]
> In the final repository, due to long exeuction times required to process large matrices, the 1D_and_2D folder has been divided in two: 
> - `/results/batch1` contains the results of the smaller matrices: `1138_bus.mtx`, `Andrews.mtx` and `bcsstk27.mtx`
> - `/results/batch2` contains the results of the biggest matrices: `Hook_1498.mtx`, `msdoor.mtx`, `Andrews.mtx`
> In case you want to re run the testbenches, the instructions above for obtaining results and plots remain unchanged and you will obtain one single 1D_and_2D folder containing all information 
