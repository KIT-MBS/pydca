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

plmdca_protein = plmdca.PlmDCA('protein', msa_file, 
    seqid=0.7,
    lambda_h=0.5, 
    lambda_J = 50.0,
    num_threads=6, 
)

# compute sorted Frobenius norm of the couplings, average product corrected 
fn_scores_apc = plmdca_protein.compute_sorted_FN_APC()

# compute sorted Frobenius norm of couplings, without average product correction
fn_scores_raw = plmdca_protein.compute_sorted_FN() 

# compute DCA scores summarized by direct information (DI), average product corrected. 