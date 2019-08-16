from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.fasta_reader.fasta_reader import get_alignment_from_fasta_file
import ctypes 
import logging
import os
import glob 

"""Implements Python wrapper for the plmDCA implemented in c++. In addition,
it defines command line interface to carry out DCA computation using the 
pseudolikelihood maximization algorithm.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)


class PlmDCAException(Exception):
    """Implements exceptions related to PlmDCA computation
    """

class PlmDCA:
    """Wraps  the c++ implementation of plmDCA.

    Attributes
    ----------
        plmdca_so : str 
            Path to the plmDCA shared object created from the C++ source code
    """
    plmdca_so_paths  = glob.glob(
            os.path.join(os.path.dirname(__file__), '_plmdcaBackend*')
    )
    logger.info('\n\tplmdca backend path: {}'.format(plmdca_so_paths))
    try:
        plmdca_lib_path = plmdca_so_paths[0]
    except IndexError:
        logger.error('\n\tUnable to find plmdca dynamic library path.' 
            '\nAre you running  pydca as a module before installation?'
            '\nIn this case you need to build the plmdca shared object.'
        )
        raise 


    def __init__(self, biomolecule, msa_file, seqid=None, lambda_h=None, 
            lambda_J=None, num_iterations=None):
        """Initilizes plmdca instances
        """

        self.__biomolecule = biomolecule.strip().upper()
        if self.__biomolecule not in ('PROTEIN', 'RNA'):
            logger.error('\n\tInvalid biomolecule type {}'.format(self.__biomolecule))
            raise PlmDCAException
        self.__msa_file = msa_file
        self.__biomolecule_int = 1 if self.__biomolecule == 'PROTEIN' else 2
        self.__num_site_states = 21 if self.__biomolecule== 'PROTEIN' else 5
        self.__num_seqs, self.__seqs_len = self._get_num_and_len_of_seqs()
        self.__seqid = 0.8 if seqid is None else seqid 
        if self.__seqid <= 0 or self.__seqid > 1.0: 
            logger.error('\n\t{} is an invalid value of sequences identity (seqid) parameter'.format(self.__seqid))
            raise PlmDCAException 
        self.__lambda_h= 1 if lambda_h is None else lambda_h
        if self.__lambda_h < 0 :
            logger.error('\n\tlambda_h must be a positive number. You passed lambda_h={}'.format(self.__lambda_h))
            raise PlmDCAException  
        self.__lambda_J= 0.2*(self.__seqs_len - 1) if lambda_J is None else lambda_J
        if self.__lambda_J < 0: 
            logger.error('\n\tlambda_J must be a positive number. You passed lambda_J={}'.format(self.__lambda_J))
            raise PlmDCAException
        self.__num_iterations = num_iterations if num_iterations is not None else  (self.__seqs_len**2)*(self.__num_site_states **2) 
        # plmdcaBackend interface
        # extern "C" double* plmdcaBackend(const unsigned int biomolecule, unsigned int num_site_states, const char* msa_file, 
        # const unsigned int seqs_len, const double seqid, const double lambda_h, const double lambda_j, int num_iter_steps)
        self.__plmdca = ctypes.CDLL(self.plmdca_lib_path)
        self.__plmdcaBackend = self.__plmdca.plmdcaBackend 
        self.__plmdcaBackend.argtypes = (ctypes.c_uint, ctypes.c_uint, ctypes.c_char_p, ctypes.c_uint, 
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_uint)
        data_size = (self.__seqs_len * (self.__seqs_len - 1) * (self.__num_site_states ** 2))/2 + self.__seqs_len * self.__num_site_states
        self.__data_size = int(data_size) 
        self.__plmdcaBackend.restype = ctypes.POINTER(ctypes.c_double * self.__data_size)
        log_message="""Created plmDCA instance with:
            biomolecule: {}
            MSA sequence length: {}
            Total number of (unweighted) sequences: {}
            sequence identity: {}
            lambda_h: {}
            lambda_J: {}
        """.format(self.__biomolecule, self.__seqs_len, self.__num_seqs, self.__seqid, self.__lambda_h, self.__lambda_J)
        logger.info(log_message)
        return None

    
    def __del__(self):
        """
        """
        # TODO 
        pass 
    

    def _get_num_and_len_of_seqs(self):
        """Obtains the length of sequences in MSA data.

        Parameters
        ----------
            self : PlmDCA 
                Instance of PlmDCA class
        
        Returns
        -------
            num_seqs, seqs_len: tuple  
                A tuple of the number and length of sequences as read from 
                the MSA file
        """
        msa_data = get_alignment_from_fasta_file(self.__msa_file)
        num_seqs = len(msa_data)
        seqs_len = len(msa_data[0])
        return num_seqs, seqs_len
    
    
    def compute_params(self):
        """
        """
        couplings_ptr = self.__plmdcaBackend(self.__biomolecule_int, self.__num_site_states, self.__msa_file.encode('utf-8'), 
            self.__seqs_len,  self.__seqid, self.__lambda_h, self.__lambda_J, self.__num_iterations)
        #couplings = [val for val in ]
        return None #couplings 
    

if __name__ == '__main__':
    """
    """
    pass