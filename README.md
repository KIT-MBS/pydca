# About `pydca`
`pydca` is Python implementation of direct coupling analysis (DCA) of residue coevolution for protein and RNA sequence families using the **__mean-field__** and **__pseudolikelihood maximization__** algorithms. Given multiple sequence alignment (MSA) files in FASTA format, `pydca` computes the coevolutionary scores of pairs of sites in the alignment. In addition, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA. The software provides command line utilities or it can be used as a library. 

# Prerequisites
`pydca` is implemented mainly in Python with the pseudolikelihood maximization parameter inference part implemented using C++ backend for optimization. To install pydca and successfully carry out DCA computations, the following are required. 
* Python 3, version 3.5 or later.
* C++ compiler that supports C++11 (e.g. the GNU compiler collection).
* Optionally, OpenMP for multithreading support.


# Installing
To install the current version of `pydca` from PyPI, run on the command line
```bash
$ pip install pydca
```
or you can use the `install.sh` bash script as 
```bash 
$ source install.sh
```

# Using `pydca` as a Python Library
After installation, pydca can be imported into other Python source codes and used. 
[Here is IPython Notebook example](https://github.com/KIT-MBS/pydca/blob/master/examples/pydca_demo.ipynb). 
If you encounter a problem opening the Ipython Notebook example, copy and past the URL [here](https://nbviewer.jupyter.org/).

# Running `pydca` From Command Line
When `pydca` is installed, it provides three main command. Namely `pydca`, `plmdca`, and `mfdca`. 
The command `pydca` is used for tasks such as trimming alignment data before DCA computation, and 
visualization of contact maps or true positive rates. The other two command are associated with 
DCA computation with the pseudolikelihood maximization algorithm (plmDCA) or the mean-field algorithm (mfDCA).
Below we show some usage examples of all the three commands.
## Trimming MSA data 
Trim gaps by reference sequence:
```bash
$ pydca trim_by_refseq <biomolecule>  <alignment.fa>  <refseq_file.fa> --remove_all_gaps --verbose
```
Trim by percentage of gaps in MSA columns:
```bash 
$ pydca trim_by_gap_size <alignmnet.fa> --max_gap 0.9 --verbose
```
### DCA Computation
#### Using `pydca`'s Pseudolikelihood Maximization Algorithm
```bash 
$ plmdca compute_fn <biomolecule> <alignment.fa> --max_iterations 500 --num_threads 6 --apc --verbose 
```
We can also the values of regularization parameters 
```bash
$ plmdca compute_fn <biomolecule> <alignment.fa> --apc --lambda_h 1.0 --lambda_J 50.0 --verbose 
```
The command `compute_fn` computes DCA scores obtained from the Frobenius norm of the couplings. `--apc` performs
average product correction (APC). To obtain DCA scores from direct-information (DI) we replace the subcommand 
`compute_fn` by `compute_di`. 
#### Using `pydca`'s Mean-Field Algorithm 
```bash
$ mfdca compute_fn <biomolecule> <alignment.fa> --apc --pseudocount 0.5 --verbose
```
### Contact Map Visualization 
When protein/RNA sequence family has a resolved PDB structure, we can evaluate the 
performance of `pydca` by contact map visualization. Example:
```bash
$ pydca plot_contact_map <biomolecule> <PDB_chain_name> <PDB_id/PDB_file.PDB> <refseq.fa> <DCA_file.txt> --verbose  
```
### Plotting True Positive Rate
In addition to contact map we can evaluate the performance of `pydca` by plotting 
the true positive rate. 
```bash
$ pydca plot_contact_map <biomolecule> <PDB_chain_name> <PDB_id/PDB_file.PDB> <refseq.fa> <DCA_file.txt> --verbose
```
To get help message about a (sub)command  we use, for example, 
```bash
$ pydca --help
```
```bash
$ plmdca compute_fn  --help
```

# References
### If you use pydca for your work please cite the following references
1. Zerihun, MB., Pucci, F, Peter, EK, and Schug, A. <br>
pydca: v1.0: a comprehensive software for direct coupling analysis of RNA and protein sequences <br>
 Bioinformatics, btz892, doi.org/10.1093/bioinformatics/btz892

2. Morcos, F., Pagnani, A., Lunt, B., Bertolino, A., Marks, DS., Sander, C., Zecchina, R., Onuchic, JN., Hwa, T., and Weigt, M. <br>
Direct-coupling analysis of residue coevolution captures native contacts across many protein families <br>
PNAS December 6, 2011 108 (49) E1293-E1301, doi:10.1073/pnas.1111471108

2. Ekeberg, M., Lövkvist, C., Lan, Y., Weigt, M., & Aurell, E. (2013). <br>
Improved contact prediction in proteins: Using pseudolikelihoods to infer Potts models. <br>
Physical Review E, 87(1), 012707, doi:10.1103/PhysRevE.87.012707