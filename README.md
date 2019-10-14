# About `pydca`
`pydca` is Python implementation of direct coupling analysis (DCA) of residue coevolution for protein and RNA sequence families using the **__mean-field__** and **__pseudolikelihood maximization__** algorithms. Given multiple sequence alignment (MSA) files in FASTA format, `pydca` computes the coevolutionary scores of pairs of sites in the alignment. In addition, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA. The software provides command line utilities or it can be used as a library. 

# Prerequisites
`pydca` is implemented mainly in Python with the pseudolikelihood maximization parameter inference part implemented using C++ backend for optimization. To install pydca and successfully carry out DCA computations, the following are required. 
* Python 3, version 3.5 or latter.
* C++ compiler that supports C++11 (we recommend GCC).
* Optionally, OpenMP for multithreading support.


# Installing
To install the current version of `pydca` from PyPI, run on the command line
```bash
$ pip install pydca
```  
# Running DCA Computation
When `pydca` is installed, it provides a main command called `mfdca`. Executing this command displays  help messages.
```bash
(demo) $ mfdca
usage: mfdca [-h]
              {compute_di,compute_fn,compute_couplings,compute_fields,compute_fi ...
....
```
In the above, we have truncated the help message.

To compute DCA scores summarized by direct information score, we can execute
```bash
(demo) $ mfdca compute_di <biomolecule> <msa_file>  --apc --verbose
```
`<biomolecule>` takes either `protein` or `rna` (case insensitive). The `<msa_file>` is
a FASTA formated multiple sequence alignment file. The optional argument `--apc` allows to make average product correction (APC) of DCA scores and  `--verbose` triggers
logging messages to be displayed on the screen as DCA computation progresses. Thus, an exaple of DCA computation for an MSA file `alignment.fa` containing protein sequences would be:
```bash
(demo) $ mfdca compute_di protein alignment.fa --apc --verbose
```
There are a few subcommand in the pydca, e.g., to compute parameters of the global probability function, frequencies of alignment data,  DCA scores summarized by Frobenius norm.
One can lookup the help messages corresponding to a particular sub(command). Example:
```bash
(demo) $ mfdca compute_di --help
usage: mfdca compute_di [-h] [--verbose] [--apc] [--seqid SEQID]
                        [--pseudocount PSEUDOCOUNT]
                        [--refseq_file REFSEQ_FILE] [--output_dir OUTPUT_DIR]
                        [--force_seq_type]
                        biomolecule msa_file
...................................................
```
Finally, the current virtual environment can be deactivated using `exit` command.

# Commonly Used Commands
Information about the commands  and subcommands can be obtained using the `--help` optional argument from the command line.  Below is a summary of commonly used commands.
##  Computing DCA Scores
In `pydca` DCA scores can be computed from the direct information score or from the Frobenius norm of the couplings using the `compute_di` or `compute_fn` subcommands. Example:
```bash
$ mfdca compute_di <biomolecule> <msa_file> --verbose
```
or
```bash
$ mfdca  compute_fn <biomolecule> <msa_file> --verbose
```
To compute the average product corrected DCA score, we can use the `--apc` optional argument. Example:

```bash
$ mfdca compute_di <biomolecule> <msa_file> --apc --verbose
```
We can also set the values of the relative pseudocount and sequence identity. The default values are 0.5 and 0.8, respectively.

```bash
$ mfdca  compute_di <biomolecule> <msa_file> --pseudocount 0.3 --seqid 0.9 --verbose
```
Furthermore, we can supply a FASTA formatted file containing a reference sequence so that DCA scores corresponding to residue pairs of that particular sequence are computed. Example:

```bash
$ mfdca compute_di <biomolecule> <msa_file> --refseq_file <refseq_file> --verbose
```
## Plotting Contact Maps or True Positive Rates

Residue pairs ranked by DCA score can be visualized and compared with an existing PDB contact map. Example:

```bash
$ mfdca plot_contact_map <biomolecule> <pdb_chain_id>  <pdb_file> <refseq_file> <dca_file> --linear_dist 5 --contact_dist 9.5 --num_dca_contacts 100 --verbose
```
In the above the optional argument `--linear_dist` filters out residue pairs that are not at least 5 residues apart in the sequence, `--contact_dist` sets the maximum distance (in Angstrom) between two residue pairs to be considered contacts in the PDB structure, and `--num_dca_contacts` sets the number of top DCA ranked residue pairs to be taken for contact map comparison. Two residues in a PDB structure are considered to be contacts if they have at least a pair of heavy atoms that are less than the distance set by `--contact_dist ` parameter.

A similar command can be used to plot the true positive rates of DCA ranked residues pairs, except that there is no restriction on the number of DCA ranked pairs. Example:

```bash
$ mfdca plot_tp_rate <biomolecule> <pdb_chain_id> <pdb_file> <ref_seq_file> --linear_dist 5 --contact_dist 8.5 --verbose
```

## Computing Couplings and Fields

The couplings and fields can be computed in a single step or separately. To compute these parameters in one run we can use the `compute_params` subcommand. Example:

```bash
$ mfdca compute_params <biomolecule> <msa_file>  --pseudocount 0.4 --seqid 0.9 --verbose
```

or to compute the couplings alone
```bash
$ mfdca compute_couplings <biomolecule> <msa_file> --pseudocount 0.4 --seqid 0.9 --verbose
```
and the fields

```bash
$ mfdca compute_fields <biomolecule> <msa_file> --pseudocount 0.3 --seqid  0.75 --verbose
```

## Computing Single- and Pair-site Frequencies

By default `pydca` uses a relative pseudocount of 0.5 thus the frequencies are regularized. To compute non-regularized frequencies we need to set the pseudocount to zero. Example:

 ```bash
 $ mfdca compute_fi <biomolecule> <msa_file> --pseudocount 0 --verbose
 ```
The `compute_fi` subcommand computes single-site frequencies. To compute pair-site frequencies, we use the `compute_fij` subcommand.

```bash
$ mfdca compute_fij <biomolecule> <msa_file> --pseudocount 0 --verbrose
```
