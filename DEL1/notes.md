
# NOTES ON PARALLELIZATION
## CACHE

- efficient parallel code we need to assure that SEQUENTIAL PARTS ARE PERFORMANT
- TRY TO COORDINATE SEQUENTIAL PARTS in such a way that exploit cache as much as possible (eg row major, ...)
- **WARNING**:
    - **FALSE SHARING MAY BE A PROBLEM**
    -> if you use cache, pay attention to cache lines: if data is invalidated, the entire content of cache line is
    -> FORCE VARIABLES WHICH ARE ACCESSED BY DIFFERENT THREADS TO BE ON DIFFERENT CACHE LINES
    
        struct alignTo64ByteCacheLine {
            int _onCacheLine1 __attribute__((aligned(64)))
            int _onCacheLine2 __attribute__((aligned(64)))
        }

    - **BRANCH PREDICTION**
    -> Random data leads to unpredictable branches, slowing execution
    -> SORTING data can improve branch prediction and speed up execution


## OPENMP
- Thread team size control -> either define a FIXED NUMBER OF THREADS by using:
• $ echo ${OMP_NUM_THREADS} # to query the value
• $ export OMP_NUM_THREADS=4 # to set it in BASH

or define the number of threads at program level (which may be better ?)

- verify that the code you want to parallelize with openmp CAN be parallelized:
    • The variable index must have integer or pointer type (e.g., it can’t be a float).
    • The expressions start, end, and incr must have a compatible type. For example, if index is a pointer, then incr must have integer type.
    • The expressions start, end, and incr must not change during execution of the loop.
    • During execution of the loop, the variable index can only be modified by the “increment expression” in the for statement.

- look at : #pragma omp parallel for collaps(n); #pragma omp parallel for (); #pragma 
- make sure work is assigned almost equally
- to run pragma omp parallel , take a look at the cyclic schedule with static keyword (page 42 pp 9)
- maybe omp parallel with chunk division is better and more general purpose than simple pragma omp parallel for -> **dynamic or guided is better for managing chunks that change in size** -> IN THE DELIVERABLE CONSIDER ALSO DIFFERENT CHUNK SIZES -> EVALUATE WHICH CHUNK DIVISION is better 
- REALLY IMPORTANT: BARRIER DIRECTIVE




- ** NOTE: RACE CONDITIONS ARE STILL POSSIBLE EVEN IN CASE OF pragma omp parallel schedule() **
- POSSIBLE PROBLEM: A REGISTER FROM A CORE MAY WANT TO ACCESS THE MEMORY OF ANOTHER CORE (**NUMA EFFECT**): mem_barrier; cpu_barrier MAY BE SOLUTIONS TO THIS
- lscpu !


try these:
- #pragma omp parallel for reduction(+:sum) schedule(static, chunk_sizes[c])
- #pragma omp parallel for reduction(+:temp_sum)
- #pragma omp parallel for collapse(2) reduction(+:sum)
- see schedule.c 


- CONSIDER THE AVERAGE OF THE MOST RECURRENT RESULTS (NO OUTLIERS)

- consider tying the first code run directly with the cache somehow....






# REFERENCES
https://stackoverflow.com/questions/28482833/understanding-the-collapse-clause-in-openmp -> for understanding better the COLLAPSE keyword for pragmas

https://ieeexplore.ieee.org/document/10444348 -> paper to Efficient COO to CSR Conversion for Accelerating Sparse Matrix Processing on FPGA

https://arxiv.org/abs/2510.13412 -> paper to Formal Verification of COO to CSR Sparse Matrix Conversion (Invited Paper) -> https://www.cs.princeton.edu/~appel/papers/coo-csr.pdf

https://stackoverflow.com/questions/23583975/convert-coo-to-csr-format-in-c -> Convert COO to CSR format in c++ stackoverflow

https://www.intel.com/content/www/us/en/developer/articles/technical/intel64-and-ia32-architectures-optimization.html -> intel manual for caching optimization under intel processor (which I have, i guess it's still useful to refer to it even if cores are not intel)

https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp -> reference for openmp static/dynamic scheduling

https://arxiv.org/html/2502.19284v1#alg1 -> regarding the optimization of SpMV algorithm

https://michalpitr.substack.com/p/optimizing-matrix-multiplication -> inspo for perf stat functionaity

https://arxiv.org/html/2404.06047v1#S8 -> THIS ON IS A GIGANTIC SURVEY THAT SUMMARIZES SIMULATIONS