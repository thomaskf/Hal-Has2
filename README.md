# Introduction

Molecular phylogenetic studies often assume that the evolutionary process behind the divergence of homologous genes was globally stationary, reversible, and homogeneous (SRH) and that a model of evolution comprising one or several site-specific, time-reversible rate matrices (e.g., the GTR rate matrix) is sufficient to accurately model the evolution of the data over the whole tree. However, an increasing body of data suggests that evolution under globally SRH conditions is the exception, rather than the norm. Therefore, we introduce a family of non-stationary and non-homogeneous mixture models that approximate the rate Heterogeneity Across Lineages (<b>HAL</b>) and the rate Heterogeneity Across Sites (<b>HAS</b>) without the assumption of an underlying predefined statistical distribution. We also develop an algorithm for searching model space and identifying a model that is less likely to over- or under-parameterize the data.

## Installation of the software

The software was written in C++, and it has been tested under linux and MacOS platform. You need to have <b>C++</b> compiler and <b>CMake</b> installed in the machine in order to compile the source codes. The compilation steps are shown as follows:

### To compile HAL2 and HAS2

```
$ tar -zxvf Hal-Has2-x.x.x.tar.gz
$ cd Hal-Has2-x.x.x
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```

Then two executable files will appear:
* HAL2	: The HAL2 program based on normal BIC formula
* HAS2	: The HAS2 program based on normal BIC formula

### To compile HAL2-P and HAS2-P

```
$ tar -zxvf Hal-Has2-x.x.x.tar.gz
$ cd Hal-Has2-x.x.x
$ mkdir build_p
$ cd build_p
$ cmake -DFLAGS=C_REG ..
$ make -j
```

Then two executable files will appear:
* HAL2-P	: The HAL2 program based on BIC formula with penalty for convergence
* HAS2-P	: The HAS2 program based on BIC formula with penalty for convergence

### To compile static files (for Linux only)

In order to compile a static file, use the option `-DFLAGS=static`

For example, to compile static files HAL2 and HAS2:

```
$ tar -zxvf Hal-Has2-x.x.x.tar.gz
$ cd Hal-Has2-x.x.x
$ mkdir build
$ cd build
$ cmake -DFLAGS=static ..
$ make -j
```
For example, to compile static files HAL2-P and HAS2-P

```
$ tar -zxvf Hal-Has2-x.x.x.tar.gz
$ cd Hal-Has2-x.x.x
$ mkdir build_p
$ cd build_p
$ cmake -DFLAGS=C_REG,static ..
$ make -j
```

## HAL2 / HAL2-P

Every unique partition of edges represents a HAL model. Since the total number of <b>HAL</B> models is a Bell number, an exhaustive search of all models is infeasible. Therefore we use an algorithm that searches a subset of the <b>HAL</B> models to identify the optimal and/or near optimal models. HAL includes two algorithms: a bottom-up algorithm, followed by a top-down algorithm.

### Usage

```
Syntax: ./HAL2   <alignment file> <topology file> <other options>
        ./HAL2-P <alignment file> <topology file> <other options>

other options:

  <alignment file>     : Multiple alignment file

  <topology file>      : Topology file in Newick format

  -f <format of file>  : The format of the multiple alignment file
                         1 - FASTA format
                         2 - Sequential PHYLIP format (default)

  -u <# of CPUs>       : Number of CPU threads used (default: 4)

  -o <output prefix>   : Prefix for output files
                         (default: <alignment file> w/o .ext)

  -r <checkpoint file> : Resume from the last execution;
                         <checkpoint file>, which was outputted from the last
                         execution of the program, stores all the immediate
                         results (with file extension: '.chkpt.txt')

  -i <# of iterations> : Number of iterations performed for each-time
                         parameter tuning
                         (default: -1 [no limit, run until converge])

  -y <matrix file>     : Optimize the parameters according to the matrix file
                         (no search algorithm will be performed)

  -b <info criteria>   : Information Criteria to be used (default: 3)
                         1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC

  -g <gap handling>    : 1 - ignore all columns with gaps;
                         2 - treat gaps as missing data (default)

  -precise             : More precise optimisation, but needs more time.
                         (default: disable)

  -l <start matrix>    : The rate matrix to start from. It is not a file,
                         but a string representing the rate matrix.
                         For example: -l 1,1,1,1,1,2,1,1,2,2,2,2,2,2

  -q                   : 0 - Not run top-down algorithm after bottom-up algorithm
                         1 - Running top-down algorithm after bottom-up algorithm
                         (default: 1)

  -h                   : This help page

```

### Output files

* \<output prefix\>.NEW-BU.chkpt.txt, which stores all the intermediate results for resuming the program when necessary.

* \<output prefix\>.NEW-BU.\<info criteria\>.result.txt, which is the resulting optimal rate matrices after the bottom-up algorithm.

* \<output prefix\>.NEW-BU.HAL.\<info criteria\>.result.txt, which is the resulting optimal rate matrices after both the bottom-up and the top-down algorithm.

### Example of running HAL2 / HAL2-P program

Example alignment file: data.phy
Example topology file: tree.txt
To execute HAL2 program by using default setting:
```
$ ./HAL2 data.phy tree.txt
```
To execute HAL2-P program by using default setting:
```
$ ./HAL2-P data.phy tree.txt
```
The result file: data.NEW-BU.HAL.BIC.result.txt 

### Difference between HAL2 and HAL2-P

Unlike the original version of HAL program, HAL2/HAL2-P considers the occurrence of convergence (i.e. existence of the same rate matrix along different lineage paths) when computing the optimal arrangement of the rate matrices along the tree.

By default, the information criteria used in HAL2 is BIC. The formula is as follows:

$BIC = -2 ln⁡(L) + k ln(n)$

where $L$ is log-likelihood value, $k$ is degree of freedoms, and $n$ is the number of sites in the alignment.

However, we found that sometimes the convergences happened more than it should be. Therefore, in HAL2-P, it gives an additional penalty for each occurrence of convergence. The formula of the information criteria being used is:

$BIC_p = -2 ln⁡(L) + (k + mc) ln(n)$

where $L$ is log-likelihood value, $k$ is degree of freedoms, $n$ is the number of sites in the alignment, $m$ is the number of occurrences of convergence, and $c$ is a constant to control the magnitude of the penalty. According to the experiments, we found $c=1$ is appropriate. Thus we set $c=1$ in the program.


## HAS2 / HAS2-P

Once the optimal HAL model has been identified, HAS2 / HAS2-P program can be used to group the variable sites into multiple categories ($K > 1$) and apply the same HAL model to all the rate categories but with different sets of parameters. The value of $K$ is increased, by default, from 1 to 4. HAS2 (or HAS2-P) program uses the BIC (or BICp) value to identify the optimal combination of $K$ and HAS model that best fits the data. This model is referred as the optimal HAL-HAS model as it takes into account rate-heterogeneity across lineages and across sites.

### Usage

```
Syntax: ./HAS2   <alignment file> <topology file> <other options>
        ./HAS2-P <alignment file> <topology file> <other options>

  <alignment file>     : Multiple alignment file

  <topology file>      : Topology file in Newick format

other options:

  -t <HAL result file> : The resulting file from HAL program
                         If no HAL result file is provided, then all edges
                         are assumed to be the same rate group

  -f <format of file>  : The format of the multiple alignment file
                         1 - FASTA format
                         2 - Sequential PHYLIP format (default)

  -u <# of CPUs>       : Number of CPU threads used (default: 4)

  -o <output prefix>   : Prefix for output files
                         (default: <alignment file> w/o .ext)

  -r <checkpoint file> : Resume from the last execution;
                         <checkpoint file>, which was outputted from the last
                         execution of the program, stores all the immediate
                         results (with file extension: '.chkpt.txt')

  -i <# of iterations> : Number of iterations performed for each-time
                         parameter tuning
                         (default: -1 [no limit, run until converge])

  -m <min category #>  : Minimum number of site categories allowed
                         (default: 1)

  -n <max category #>  : Maximum number of site categories allowed
                         (default: 4)

  -w <step category #> : Step of change of site categories
                         (default: 1)

  -precise             : More precise optimzation, but needs more time.
                         (default: disable)

  -b <info criteria>   : Information Criteria to be used (default: 3)
                         1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 – CAIC

  -h                   : This help page
```

### Output files

* \<output prefix\>.HAS.chkpt.txt, which stores all the intermediate results for resuming the program when necessary.

* \<output prefix\>.HAS.\<info criteria\>.result.txt, which includes the optimal number of site categories, the corresponding proportion and the details of parameters for each site category.

### Example of running HAS2 / HAS2-P program

Example alignment file: data.phy
Example topology file: tree.txt
HAL result file: data.NEW-BU.HAL.BIC.result.txt
To execute HAS2 program:
```
$ ./HAS2 data.phy tree.txt –t data.NEW-BU.HAL.BIC.result.txt
```
To execute HAS2-P program:
```
$ ./HAS2-P data.phy tree.txt –t data.NEW-BU.HAL.BIC.result.txt
```
The result file: data.HAS.BIC.result.txt 

## Contact person

Dr Lars Jermiin

Email: Lars.Jermiin@universityofgalway.ie

Dr Thomas Wong

Email: Thomas.Wong@anu.edu.au

