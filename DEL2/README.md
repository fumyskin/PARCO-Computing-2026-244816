# DELIVERABLE 2 run instructions

## Overviews
This folder contains the DELIVERABLE 2 source code.
The following code is a playground to test different MPI approaches to the classical SpMV kernel.
Included in this file are the following implementation:

- **BASELINE**: 1D cyclic distribution with halo exchange management
- parallel I/O reading approach to Matrix Market files 
- 2D block and 2D cyclic distributions
- Hybrid MPI approach

However, the DELIVERABLE content will focus only on exploring, apart from the baseline, **parallel I/O reading** approach and **2D block partitioning**, studying them with **bonus metrics**.


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
│   ├── batch1
│   ├── batch2
│   └── parallel_io_reading             #plots for sequential vs parallel reading
├── results
│   ├── batch1
│   ├── batch2
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
To submit the .pbs that executes 1D and 2D block distributions on SpMV :
```bash
qsub mpi_job_script_baseline.pbs
```

To submit the .pbs that executes sequential vs parallel I/O reading of Matrix Market file :
```bash
qsub mpi_job_script_parallel.pbs
```

> [!WARNING]
> The mpi_job_script_baseline.pbs will run also the 2D_cyclic and HYRBID MPI experiments as well. Since the cluster HP2 didn't work anymore I preferred to kept the original experimental setup as it is in order to be sure the source code is reproducible