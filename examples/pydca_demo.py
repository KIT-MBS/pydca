#import pydca  or selected modules
import pydca 
from pydca.plmdca import plmdca 
from pydca.meanfield_dca import meanfield_dca 
from pydca.sequence_backmapper import sequence_backmapper
from pydca.msa_trimmer import msa_trimmer 

"""Demonstrates usage of pydca within Python scrips
"""

protein_msa_file = ''
protein_refseq_file = ''
rna_msa_file = ''
rna_refseq_file = ''


# Creating PlmDCA instance and using it 

plmdca = plmdca.PlmDCA('protein', msa_file, 
    seqid=0.7,
    lambda_h=0.5, 
    lambda_J = 50.0,
    num_threads=6, 
)

# In the above we created plmdca instance for a protein sequence whose MSA data
# is going to be read from the parameter msa_file. In addition we have set the 
# values of optional paramters. seqid sets the threshold value for sequence 
# similarity, lambda_h sets the values of L2 norm regularization parameter for fields
# and lambda_J for that of the the couplings. Finally we set the number of threads 
# for parallel execution to be 6. Note that we can only use multiple threads 
# only when pydca installation was done with OpenMP support. If OpenMP is not 
# supported we can still carry out plmDCA computation althought it can only be 
# done using a single threads which might be slow for large alignment data.


# compute sorted Frobenius norm of the couplings, average product corrected 
fn_scores_apc = plmdca.compute_sorted_FN_APC()

# compute sorted Frobenius norm of couplings, without average product correction
fn_scores_raw = plmdca.compute_sorted_FN() 

# compute DCA scores summarized by direct information (DI), average product corrected. 