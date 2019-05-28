# About `pydca`
`pydca` is a Python implementation of direct coupling analysis (DCA) of residue coevolution for protein and RNA sequence families using the mean-field approximation algorithm. Given multiple sequence alignment (MSA) files in FASTA format, `pydca` computes the coevolutionary scores of pairs of sites in the alignment. In addition, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA. Furthermore, `pydca`  provides commands to compute the parameters of the energy function, compare and visualize contact map or true positive rates of DCA predicted residue pairs.

# Prerequisites
* Python3, version 3.5 or latter.
* python3-tk


# Installing
Here we show installation of `pydca` in a Linux (Ubuntu) machine. A similar procedure is used for other Linux/Unix machines as well as Mac or Windows platforms.
### Installing python3-tk
```bash
$ sudo apt install python3-tk
```
### Create a virtual environment for Python3
`pydca` depends on python based libraries, e.g., `numpy`, `scipy`, `matplotlib`, or `biopython`. To avoid conflicts that can arise between different versions of libraries, its recommended that we install `pydca`  separate from other global installations. Here, we use a separate virtual environment using [pipenv](https://docs.pipenv.org/en/latest/). To install `pipenv` use use `pip`. Thus, first we need to install `pip`  for `Python3`.

```bash
$ sudo apt install python3-pip
```  
Now we can install pipenv as follows
```bash
$ pip3 install pipenv
```
Once `pipenv` is installed we can install `pydca` using `pipenv`. To demonstrate this, lets create a
directory called `demo`, you can call it another name. We will create this directory in current user's `home` directory. Then, `cd` to this directory and install `pydca` using `pipenv`.
```bash
$ mkdir demo && cd demo
```
We can see our current directory's path using `pwd` command.
```bash
$ pwd
$ /home/username/demo
```
Now we install pydca
```bash
$ pipenv install pydca
```
At this time, `pipenv` creates a virtual environment for us and installs `pydca` in this newly created virtual environment. After installation is completed `pipenv` has created two files in
the directory `/home/username/demo`. We can see them as
```bash
$ ls
Pipfile  Pipfile.lock
```
These files contain metadata about the virtual environment. Now we can activate the virtual environment
```bash
$ pipenv shell
Launching subshell in virtual environment…
$  . /home/username/.local/share/virtualenvs/demo-HO4EpenA/bin/activate
(demo) $
```
The virtual environment is now activated. Lets check that we are executing Python from the virtual environment by starting Python shell, importing `sys` module and executing `sys.executable`.
```bash
(demo) $ python
Python 3.5.2 (default, Nov 12 2018, 13:43:14)
[GCC 5.4.0 20160609] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> sys.executable
'/home/username/.local/share/virtualenvs/demo-HO4EpenA/bin/python'
>>>
```
Lets exit the Python shell. We do this using the `exit()` command.
```bash
>>> exit()
(demo) $
```
Now we are back in the virtual environment's shell.
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
