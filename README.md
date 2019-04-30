# About pydca
pydca is a Python implementation of direct coupling analysis of residue coevolution for proteins and RNAs using the mean-field approximation algorithm. Given multiple sequence alignment (MSA) files in FASTA format, pydca computes the coevolutionary scores of pairs of sites in the alignment. Furthermore, when an optional file containing a reference sequence is supplied, scores corresponding to pairs of sites of this reference sequence are computed by mapping the reference sequence to the MSA.

# Prerequisites
* Python3, version 3.5 or latter.
* python3-tk


# Installing
Here, we demonstrate how to install pydca in a Linux/Unix machine. A similar procedure follows for Mac or Windows.

In order to avoid conflicts that might arise from different version of dependencies its recommended to install pydca as a separate installation, e.g., in a separate virtual environment. Here is an example of installing pydca in a separate virtual environment.  

### Create a virtual environment for Python3.
I will call the name of the virtual environment *python3venv*.
```bash
$ virtualenv -p python3 python3venv
```  
### Activate the virtual environment
```bash
$ source python3venv/bin/activate
(python3venv) $
```
The virual environment, *python3venv* is activated.
Now lets check if we are really using Python from the activated virtual environment. For this, first I start Python shell using `python` command, `ìmport sys` module and run `sys.executable` as follows.
```bash
(python3venv) $ python
Python 3.5.2 (default, Nov 12 2018, 13:43:14)
[GCC 5.4.0 20160609] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> sys.executable
'/home/mehari/python3venv/bin/python'
>>>
```
The output of `sys.executable`, i.e., *'/home/mehari/python3venv/bin/python'* shows that I am using Python from the activate virtual environment. Now, we can exit the Python shell using `exit()` command. Lets check installed softwares in the current virtual environment using `pip freeze` as:
```bash
(python3venv) $ pip freeze
(python3venv) $
```
As we can see the virtual environment is clean, i.e., `pip freeze` outputs nothing.

### Install pydca

Now we can install pydca in the current virtual environment from PyPi as:
```bash
(python3venv) $ pip install pydca
```
You can also clone/download the [current version of pydca](https://github.com/KIT-MBS/pydca) and manually install it.

# Running DCA computation
===============
