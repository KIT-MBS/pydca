# About pymfDCA
pymfDCA is a Python implementation of direct coupling analysis of residue coevolution for proteins and RNAs using the mean-field approximation algorithm. Given multiple sequence alignment (MSA) files in FASTA format, pymDCA computes the coevolutionary scores of pairs of sites in the alignment. Furthermore, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA.

# Prerequisites
pymfDCA requires Python3, version 3.5 or latter.

# Installing
In order to avoid conflicts that might arise from different version of dependencies its recommended to install pymfDCA as a separate installation, e.g., in a separate virtual environment. Here is an example of installing pymfDCA in a separate virtual environment.  

1. Create a virtual environment for Python3. I will call the name of the virtual environment "demo".
```bash
vitrualenv -p python3 demo
```  
2. Activate the virtual environment
```bash
source demo/bin/activate
```
3. Now lets check if we are using Python from the newly created virtual environment.
