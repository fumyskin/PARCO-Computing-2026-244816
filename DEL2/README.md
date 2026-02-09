# DELIVERABLE 2 run instructions

This folder contains the DELIVERABLE 2 source code.
The following code is a playground to test different MPI approaches to the classical SpMV kernel.
Included in this file are the following implementation:

- **BASELINE**: 1D cyclic distribution with halo exchange management
- parallel I/O reading approach to Matrix Market files 
- 2D block and 2D cyclic distributions
- Hybrid MPI approach

However, the DELIVERABLE content will focus only on exploring, apart from the baseline, **parallel I/O reading** approach and **2D block partitioning**, studying them with **bonus metrics**.


## FOLDER ORGANIZATION
inserisci diagramma ad albero della repo
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
├── matrices
│   └── 1138_bus.mtx
├── plots
│   ├── batch1
│   ├── batch2
│   └── parallel_io_reading
├── results
│   ├── batch1
│   ├── batch2
│   └── parallel_io_reading
├── scripts
│   ├── mpi_bash_baseline.sh            #bash script for main_baseline.c
│   ├── mpi_bash_parallel.sh            #bash script for main_parallel.c
│   ├── mpi_job_script_baseline.pbs     #pbs script for mpi_bash_baseline.sh
│   ├── mpi_job_script_parallel.pbs     #pbs script for mpo_bash_parallel.sh
│   ├── plot_parallel.py                #python script for parallel_io_reading results
│   └── plot_script.py                  #python script for 1D_and_2D results
└── src
    ├── communication.c
    ├── distributions.c
    ├── main_baseline.c
    ├── main_parallel.c
    ├── mmio.c
    ├── parallel_reading.c
    ├── specifications.c
    └── tests.c
```

## HOW TO REPRODUCE
