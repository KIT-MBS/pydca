# About pydca
pydca is a Python implementation of direct coupling analysis of residue coevolution for proteins and RNAs using the mean-field approximation algorithm. Given multiple sequence alignment (MSA) files in FASTA format, pydca computes the coevolutionary scores of pairs of sites in the alignment. Furthermore, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA.

# Prerequisites
* Python3, version 3.5 or latter.
* python3-tk


# Installing
Before we install pydca, we need to make sure that we have Python3 interpreter and python3-tk installed. Here we show installation in an ubuntu/debian machine. A similar procedure is used for other Linux/Unix machines as well as Mac or Windows platforms.
### Installing python3-tk
```bash
$ sudo apt install python3-tk
```
### Create a virtual environment for Python3
pydca depends python based libraries, e.g., *numpy*, *scipy*, *matplotlib*, or *biopython*. To avoid conflicts that can arise between different versions of these libraries, its recommended to install pydca separate from other global installations. Here, we use a separate virtual environment using [pipenv](https://docs.pipenv.org/en/latest/). To install *pipenv* use use *pip*. Thus, we need to install *pip*  for Python3.
```bash
$ sudo apt install python3-pip
```  
Now we can install pipenv as follows
### Activate the virtual environment
```bash
$ pip3 install pipenv
```
Once `pipenv` is installed we can install pydca using pipenv. To demonstrate this, lets create a
directory called `demo` in, you can call it another name. We will create this directory in current user's `home` directory. Then, `cd` to this directory and install `pydca` using `pipenv`.
```bash
$ mkdir demo && cd demo
```
We can see our current directory's path using `pwd` command.
```bash
$ pwd
$ /home/mehari/demo
```
Now we install pydca
```bash
$ pipenv install pydca
```
At this time, `pipenv` creates a virtual environment for us and installs pydca in this newly created virtual environment. After installation is completed `pipenv` has created two files in
the directory `/home/mehari/demo`. We can see them as
```bash
$ ls
Pipfile  Pipfile.lock
```
These files contain metadata about virtual environment. Now we can activate the virtual environment
```bash
$ pipenv shell
Launching subshell in virtual environment…
$  . /home/mehari/.local/share/virtualenvs/demo-HO4EpenA/bin/activate
(demo) $
```
The virtual environment is now activated. Lets check that we are executing Python from the virtual environment by starting python shell, importing `sys` module and executing `sys.executable`.
```bash
(demo) $ python
Python 3.5.2 (default, Nov 12 2018, 13:43:14)
[GCC 5.4.0 20160609] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> sys.executable
'/home/mehari/.local/share/virtualenvs/demo-HO4EpenA/bin/python'
>>>
```
Lets exit the Python shell. We do this using the `exit()` command.
```bash
>>> exit()
(demo) $
```
Now we are back in the virtual environment's shell.
# Running DCA computation
When `pydca` is installed, it provides main command called `mfdca`. Executing this command displays messages.
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
a FASTA formated multiple sequence alignment file. The optional argument `--apc` allows to make average product correction of DCA scores and  `--verbose` triggers
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
