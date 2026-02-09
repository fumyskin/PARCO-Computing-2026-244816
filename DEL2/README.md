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
