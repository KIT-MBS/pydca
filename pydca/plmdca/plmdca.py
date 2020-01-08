#from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.fasta_reader.fasta_reader import get_alignment_from_fasta_file
from pydca.fasta_reader.fasta_reader import get_alignment_int_form
from . import msa_numerics
import ctypes 
import logging
import os
import glob 
import numpy as np


"""Python wrapper for psuedolikelihood maximization direct coupling analysis (plmDCA).
The gradient decent algorithm is implemented using c++ backend.

Authors: Mehari B. Zerihun, Fabrizio Pucci
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


    def __init__(self, msa_file, biomolecule, seqid = None, lambda_h=None, 
            lambda_J = None, max_iterations = None, num_threads = None, 
            verbose = False):
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
        self.__lambda_h= 0.2*(self.__seqs_len - 1) if lambda_h is None else lambda_h
        if self.__lambda_h < 0 :
            logger.error('\n\tlambda_h must be a positive number. You passed lambda_h={}'.format(self.__lambda_h))
            raise PlmDCAException  
        self.__lambda_J=  0.2*(self.__seqs_len - 1) if lambda_J is None else lambda_J
        if self.__lambda_J < 0: 
            logger.error('\n\tlambda_J must be a positive number. You passed lambda_J={}'.format(self.__lambda_J))
            raise PlmDCAException
        self.__max_iterations = max_iterations if max_iterations is not None else 100
        self.__num_threads = 1 if num_threads is None else num_threads
        self.__verbose = True if verbose else False  
        # plmdcaBackend interface
        # extern "C" float* plmdcaBackend(unsigned short const biomolecule, unsigned short const num_site_states, 
        # const char* msa_file, unsigned int const seqs_len, float const seqid, float const lambda_h, 
        # float const lambda_J, unsigned int const max_iterations)
        self.__plmdca = ctypes.CDLL(self.plmdca_lib_path)
        self.__plmdcaBackend = self.__plmdca.plmdcaBackend 
        self.__plmdcaBackend.argtypes = (ctypes.c_ushort, ctypes.c_ushort, ctypes.c_char_p, ctypes.c_uint, 
            ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_uint, ctypes.c_uint, ctypes.c_bool)
        data_size = (self.__seqs_len * (self.__seqs_len - 1) * (self.__num_site_states ** 2))/2 + self.__seqs_len * self.__num_site_states
        self.__data_size = int(data_size) 
        self.__plmdcaBackend.restype = ctypes.POINTER(ctypes.c_float * self.__data_size)
        #extern "C" void freeFieldsAndCouplings(float*& fields_and_couplings)
        self.freeFieldsAndCouplings = self.__plmdca.freeFieldsAndCouplings
        #self.freeFieldsAndCouplings.argtypes = (ctypes.POINTER(ctypes.c_float),) 
        self.freeFieldsAndCouplings.restype = None 
        log_message="""Created plmDCA instance with:
            biomolecule: {}
            MSA sequence length: {}
            total number of (unweighted) sequences: {}
            sequence identity: {}
            lambda_h: {}
            lambda_J: {}
            gradient decent iterations: {}
            number of threads: {}
        """.format(self.__biomolecule, self.__seqs_len, self.__num_seqs, 
            self.__seqid, self.__lambda_h, self.__lambda_J, self.__max_iterations,
            self.__num_threads,
        )
        logger.info(log_message)
        return None


    @property
    def biomolecule(self):
        """PlmDCA biomolecule attribute 
        """
        return self.__biomolecule


    @property
    def sequence_identity(self):
        """PlmDCA sequence_identity attribute
        """
        return self.__seqid


    @property
    def lambda_h(self):
        """PlmDCA lambda_h attribute
        """
        return self.__lambda_h
    

    @property
    def lambda_J(self):
        """PlmDCA lambda_J attribute
        """
        return self.__lambda_J


    @property
    def max_iterations(self):
        """PlmDCA max_iterations attribute
        """
        return self.__max_iterations


    @property
    def sequences_len(self):
        """PlmDCA sequences_len attribute
        """
        return self.__seqs_len


    @property
    def num_sequences(self):
        """PlmDCA num_sequences attribute
        """
        return self.__num_seqs
    
    
    @property
    def effective_num_sequences(self):
        """PlmDCA effective_num_sequences attribute
        """
        raise NotImplementedError
        

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


    def map_index_couplings(self, i, j, a, b):
        """Couplings index mapper.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class.
            i : int 
                Site in site-pair (i, j) such that j > i. 
            j : int 
                Site in site-pair (i, j) such that j > i.
        """
        q = self.__num_site_states
        L = self.__seqs_len
        site = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q)
        k =  L * q + site + b + a * q
        return k
    
    
    def get_fields_and_couplings_from_backend(self):
        """Compute fields and couplings using the C++ backend of plmDCA

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class

        Returns
        -------
            fields_and_couplings : np.array
                A one-dimensional array of the fields and couplings
        """
        logger.info('\n\tComputing fields and couplings using gradient decent')
        h_J_ptr = self.__plmdcaBackend(
            self.__biomolecule_int, self.__num_site_states, self.__msa_file.encode('utf-8'), 
            self.__seqs_len,  self.__seqid, self.__lambda_h, self.__lambda_J, self.__max_iterations,
            self.__num_threads, self.__verbose
        )
        
        fields_and_couplings = np.zeros((self.__data_size,), dtype=np.float32)
        #fields_and_couplings = [c for c in h_J_ptr.contents]
        counter = 0
        for i, h_or_J in enumerate(h_J_ptr.contents):
            fields_and_couplings[i] = h_or_J
            counter += 1
        #Free fields and couplings data from PlmDCABackend
        h_J_ptr_casted = ctypes.cast(h_J_ptr, ctypes.POINTER(ctypes.c_void_p))
        self.freeFieldsAndCouplings(h_J_ptr_casted)

        #fields_and_couplings = np.array(fields_and_couplings, dtype=np.float32)
        
        try:
            assert fields_and_couplings.size == counter 
        except AssertionError:
            logger.error('\n\tData size mismatch from the plmDCA backend')
            raise
        else:
            logger.info('\n\tData size expected from plmDCA backend: {}'.format(self.__data_size))
            logger.info('\n\tData size obtained from plmDCA backend to Python: {}'.format(fields_and_couplings.size))

        return fields_and_couplings 


    def get_couplings_no_gap_state(self, fields_and_couplings_all):
        """Extract the couplings from fields and couplings numpy array

        Parameters
        ----------
            self : PlmDCA 
                An instanc of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        Returns
        -------
            couplings : np.array
                A one-dimensional array of the couplings
        """
        couplings = list()
        for i in range(self.__seqs_len - 1):
            for j in range(i + 1, self.__seqs_len):
                for a in range(self.__num_site_states - 1):
                    for b in range(self.__num_site_states -1):
                        indx =  self.map_index_couplings(i, j , a, b)
                        couplings.append(fields_and_couplings_all[indx])
        return np.array(couplings)

    
    def get_fields_no_gap_state(self, fields_and_couplings_all):
        """Extracts the fields from fields and couplings numpy array.

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        Returns
        --------
            fields_no_gap_state : list 
                A  list of fields excluding fields corresponding to gap states.
        """
        
        fields_all = fields_and_couplings_all[:self.__seqs_len * self.__num_site_states]
        fields_no_gap_state = list()
        for i in range(self.__seqs_len):
            for a in range(self.__num_site_states - 1): # iterate over q - 1 states to exclude gaps
                fields_no_gap_state.append(fields_all[a + i * self.__num_site_states])
        return fields_no_gap_state

    
    def get_fields_and_couplings_no_gap_state(self, fields_and_couplings_all):
        """Computes fields and couplings excluding gap state fields and gap state
        couplings.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        
        Returns
        --------
            fields_no_gap_state, couplings_no_gap_state : tuple
                A tuple of list of fields and couplings
        """

        
        logger.info('\n\tObtaining fields and couplings excluding gap state.')
        fields_no_gap_state = self.get_fields_no_gap_state(fields_and_couplings_all)
        couplings_no_gap_state = self.get_couplings_no_gap_state(fields_and_couplings_all)
        
        return fields_no_gap_state, couplings_no_gap_state

    def shift_couplings(self, couplings_ij):
        """Shifts the couplings value.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class
            couplings_ij : np.array
                1d array of couplings for site pair (i, j)
        Returns
        -------
            shifted_couplings_ij : np.array
                A 2d array of the couplings for site pair (i, j)
        """
        qm1 = self.__num_site_states - 1
        couplings_ij = np.reshape(couplings_ij, (qm1,qm1))
        avx = np.mean(couplings_ij, axis=1)
        avx = np.reshape(avx, (qm1, 1))
        avy = np.mean(couplings_ij, axis=0)
        avy = np.reshape(avy, (1, qm1))
        av = np.mean(couplings_ij)
        couplings_ij = couplings_ij -  avx - avy + av
        return couplings_ij 


    def compute_params(self, seqbackmapper=None, ranked_by = None, linear_dist = None, num_site_pairs =None):
        """Computes fields and couplings with the couplings ranked by DCA score.

        Parameters
        ----------
            self : PlmDCA
                An instanc of PlmDCA class
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class
            ranked_by : str
                DCA score type usef to rank the couplings by their site pairs.
                By default they are ranked by the Frobenius Norm of couplings with
                average product correction.
            linear_dist : int
                Minimum separation beteween site pairs (i, j).
            num_site_pairs : int 
                Number of site pairs whose couplings are to be otained. 
        
        Returns
        -------
            fields, couplings : tuple 
                A tuple of lists of fields and couplings. 
        """
        if ranked_by is None: ranked_by = 'fn_apc'
        if linear_dist is None: linear_dist = 4

        RANKING_METHODS = ('FN', 'FN_APC', 'DI', 'DI_APC')
        ranked_by = ranked_by.strip().upper()
        if ranked_by not in RANKING_METHODS:
            logger.error('\n\tInvalid ranking criterion {}.\nChoose from {}'.format(ranked_by, RANKING_METHODS))
            raise PlmDCAException
        if ranked_by == 'FN': dca_scores = self.compute_sorted_FN(seqbackmapper=seqbackmapper)
        if ranked_by == 'FN_APC': dca_scores = self.compute_sorted_FN_APC(seqbackmapper=seqbackmapper)
        if ranked_by == 'DI': dca_scores = self.compute_sorted_DI(seqbackmapper=seqbackmapper)
        if ranked_by == 'DI_APC': dca_scores = self.compute_sorted_DI_APC(seqbackmapper=seqbackmapper)

        fields = self.get_fields_no_gap_state(self.__fields_and_couplings_all)
        couplings = self.get_couplings_no_gap_state(self.__fields_and_couplings_all)
        
        qm1 = self.__num_site_states - 1 
        L = self.__seqs_len
        if seqbackmapper is not None:
            # mapping_dict has keys from MSA sites and values from refseq sites
            # we need to reverse this mapping as the fields and couplings are from MSA sites
            mapping_dict = {
                value : key for key, value in self.__refseq_mapping_dict.items()
            }
        else:
            mapping_dict = {
                i : i for i in range(self.__seqs_len)
            }
        # set default number of site pairs whose couplings are to be extracted
        if num_site_pairs is None :
            num_site_pairs = len(seqbackmapper.ref_sequence) if seqbackmapper is not None else len(mapping_dict.keys()) 
        # we need only the fields corresponding to mapped sites 
        fields_mapped = list()
        logger.info('\n\tExtracting fields')
        for i in mapping_dict.keys():
            site_in_msa = mapping_dict[i]
            fields_im = fields[qm1 * site_in_msa : qm1 * site_in_msa + qm1]
            site_fields = i, fields_im
            fields_mapped.append(site_fields)
        # extract couplings
        logger.info('\n\tExtracting couplings for top {} site pairs (i, j) with |i - j| > {} and ranked by {}'.format(
            num_site_pairs, linear_dist, ranked_by)
        )
        couplings_ranked_by_dca_score = list()
        count_pairs = 0
        for pair, score in dca_scores:
            site_1_in_refseq, site_2_in_refseq = pair[0], pair[1]
            if abs(site_1_in_refseq - site_2_in_refseq) > linear_dist:
                count_pairs += 1
                if count_pairs > num_site_pairs: break 
                i, j = mapping_dict[site_1_in_refseq], mapping_dict[site_2_in_refseq]
                if(i > j): 
                    logger.error('\n\tInvalid site pair. Site pair (i, j) should be ordered in i < j')
                    raise PlmDCAException
                start_indx = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * qm1 * qm1)
                end_indx = start_indx + qm1 * qm1
                couplings_ij = couplings[start_indx:end_indx]
                couplings_ij = self.shift_couplings(couplings_ij)
                couplings_ij = np.reshape(couplings_ij, (qm1*qm1,))
                pair_couplings_ij = pair, couplings_ij   
                couplings_ranked_by_dca_score.append(pair_couplings_ij)
        if count_pairs < num_site_pairs:
            logger.warning('\n\tObtained couplings for only {} ranked site pairs.' 
                '\n\tThis is the maximum number of site paris we can obtain under ' 
                'the given criteria'.format(count_pairs)
            )
        return tuple(fields_mapped), tuple(couplings_ranked_by_dca_score) 


    def compute_sorted_FN(self, seqbackmapper=None):
        """Computes DCA scores using the Frobenius norm of the couplings. Note 
        that this score is not average product corrected.

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            
        Returns
        -------
            dca_scores_not_apc : List of tuples of site pairs and DCA scores.
                The contents are sorted according to DCA score, in descending 
                order.
        """
        dca_scores_not_apc = list()
        L = self.__seqs_len
        q = self.__num_site_states
        fields_and_couplings = self.get_fields_and_couplings_from_backend()
        # add attribute fields and couplings
        self.__fields_and_couplings_all = fields_and_couplings
        couplings = self.get_couplings_no_gap_state(fields_and_couplings)
        qm1 = q - 1
        logger.info('\n\tComputing non-APC sorted DCA score') 
        for i in range(self.__seqs_len -1):
            for j in range(i + 1, self.__seqs_len):
                start_indx = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * qm1 * qm1)
                end_indx = start_indx + qm1 * qm1
                couplings_ij = couplings[start_indx:end_indx]
                
                couplings_ij = np.reshape(couplings_ij, (qm1,qm1))
                avx = np.mean(couplings_ij, axis=1)
                avx = np.reshape(avx, (qm1, 1))
                avy = np.mean(couplings_ij, axis=0)
                avy = np.reshape(avy, (1, qm1))
                av = np.mean(couplings_ij)
                couplings_ij = couplings_ij -  avx - avy + av
                dca_score = np.sum(couplings_ij * couplings_ij)
                dca_score = np.sqrt(dca_score)
                data = ((i, j), dca_score)
                dca_scores_not_apc.append(data)
        dca_scores_not_apc = sorted(dca_scores_not_apc, key=lambda k : k[1], reverse=True)
        if seqbackmapper is not None:
            dca_scores_not_apc = self.get_mapped_site_pairs_dca_scores(dca_scores_not_apc, seqbackmapper)
        return dca_scores_not_apc


    def compute_sorted_FN_APC(self, seqbackmapper=None):
        """Performs average product correction (APC) of the Frobenius norm of the
        couplings for the plmDCA. If A SequenceBackmapper object is provided, 
        dca scores corresponding to mapped site-pairs are computed.

        Parameters
        ---------
            self : PlmDCA
                An instance of PlmDCA class
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns 
        -------
            sorted_FN_APC : list of tuples
                A list of tuples of  site of pairs and APC DCA scores
        """ 
    
        # compute the average score of each site
        av_score_sites = list()
        N = self.__seqs_len
        scores_plmdca = self.compute_sorted_FN()
        logger.info('\n\tPerforming average product correction (APC) of FN  of DCA scores')
        for i in range(N):
            i_scores = [score for pair, score in scores_plmdca if i in pair]
            assert len(i_scores) == N - 1
            i_scores_sum = sum(i_scores)
            i_scores_ave = i_scores_sum/float(N - 1)
            av_score_sites.append(i_scores_ave)
        # compute average product corrected DI
        av_all_scores = sum(av_score_sites)/float(N)
        sorted_FN_APC = list()
        for pair, score in scores_plmdca:
            i, j = pair
            score_apc = score - av_score_sites[i] * (av_score_sites[j]/av_all_scores)
            sorted_FN_APC.append((pair, score_apc))
        # sort the scores as doing APC may have disrupted the ordering
        sorted_FN_APC = sorted(sorted_FN_APC, key = lambda k : k[1], reverse=True)
        if seqbackmapper is not None:
            sorted_FN_APC = self.get_mapped_site_pairs_dca_scores(sorted_FN_APC, seqbackmapper)
        return sorted_FN_APC

    
    def  get_mapped_site_pairs_dca_scores(self, sorted_dca_scores, seqbackmapper):
        """Filters mapped site pairs with a reference sequence. 

        Parameters
        -----------
            self : PlmDCA
                An instance of PlmDCA class
            sorted_dca_scores : tuple of tuples
                A tuple of tuples of site-pair and DCA score sorted by DCA scores 
                in reverse order.
            seqbackmapper : SequenceBackmapper 
                An instance of SequenceBackmapper class
        
        Returns
        -------
            sorted_scores_mapped : tuple
                A tuple of tuples of site pairs and dca score
        """
        mapping_dict = seqbackmapper.map_to_reference_sequence()
        # Add attribute __reseq_mapping_dict
        self.__refseq_mapping_dict = mapping_dict 
        sorted_scores_mapped = list()
        num_mapped_pairs = 0
        for pair, score in sorted_dca_scores:
            try:
                mapped_pair = mapping_dict[pair[0]], mapping_dict[pair[1]]
            except  KeyError:
                pass 
            else:
                current_pair_score = mapped_pair, score 
                sorted_scores_mapped.append(current_pair_score)
                num_mapped_pairs += 1
        # sort mapped pairs in case they were not
        sorted_scores_mapped = sorted(sorted_scores_mapped, key = lambda k : k[1], reverse=True)
        logger.info('\n\tTotal number of mapped sites: {}'.format(num_mapped_pairs))
        return tuple(sorted_scores_mapped)

    
    def compute_seqs_weight(self):
        """Computes sequences weight

        Parameters
        ----------
            self: PlmDCA
                An instance of PlmDCA class
        
        Returns 
        -------
            seqs_weight : np.array
                A 1d numpy array containing sequences weight.
        """
        logger.info('\n\tComputing sequences weight with sequence identity {}'.format(self.__seqid))
        alignment_data = np.array(
            get_alignment_int_form(self.__msa_file, 
            biomolecule = self.__biomolecule)
        )
        seqs_weight = msa_numerics.compute_sequences_weight(
            alignment_data=alignment_data, 
            sequence_identity = self.__seqid
        )
        Meff = np.sum(seqs_weight)
        logger.info('\n\tEffective number of sequences: {}'.format(Meff))
        self.__seqs_weight = seqs_weight 
        self.__eff_num_seqs = Meff
        return seqs_weight 


    def get_single_site_freqs(self):
        """Computes single site frequencies from MSA data

        Parameters 
        ----------
            self : PlmDCA
                An instance of PlmDCA class
        
        Returns
        -------
            single_site_freqs :
                A 2d numpy array of type float64. The shape of this array is
                (seqs_len, num_site_states) where seqs_len is the length of sequences
                in the alignment data.
        """
        alignment_data = np.array(get_alignment_int_form(self.__msa_file, biomolecule = self.__biomolecule))
        seqs_weight = msa_numerics.compute_sequences_weight(
            alignment_data=alignment_data, 
            sequence_identity = self.__seqid
        )
        logger.info('\n\tComputing single site frequencies')
        
        single_site_freqs = msa_numerics.compute_single_site_freqs(
            alignment_data = alignment_data, 
            num_site_states = self.__num_site_states, 
            seqs_weight = seqs_weight 
        )
        return single_site_freqs
      

    def get_reg_single_site_freqs(self):
        """Regularizes single-site frequencies
        
        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
        
        Returns
        -------
            fi_reg : np.array 
                A 2d numpy array of shape (seqs_len, num_site_states) of single site
                frequencies after they are regularized.
        """
        # TODO  add attribute pseudocount at the __init__ method for DI computation.
        # Currently DI scores for plmDCA are computed using pseudocount = 0.5 
        
        fi = self.get_single_site_freqs()
        logger.info('\n\tComputing regularized single-site frequencies')
        fi_reg = msa_numerics.get_reg_single_site_freqs(
            single_site_freqs = fi, 
            seqs_len = self.__seqs_len, 
            num_site_states = self.__num_site_states, 
            pseudocount = 0.5
        )
        return fi_reg 

    
    def compute_two_site_model_fields(self, couplings):
        """Computes two-site model fields. These fields are used for the 
        computation of the direct information.

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            couplings : np.array 
                A 1 D array of couplings
            
        Returns
        -------
            two_site_model_fields : np.array 
                A numpy array of shape (P, 2, num_site_states), where P is the number
                of unique site pairs excluding self pairings.
                P = seqs_len * (seqs_len - 1)/2.
        """
        
        reg_fi = self.get_reg_single_site_freqs()

        logger.info('\n\tComputing two-site model fields') 
        two_site_model_fields = msa_numerics.compute_two_site_model_fields(
            couplings = couplings, 
            reg_fi = reg_fi, 
            seqs_len = self.__seqs_len, 
            num_site_states = self.__num_site_states
            )
        return two_site_model_fields


    def compute_direct_info_unsorted_DI(self):
        """Compute plmDCA direct information score

        Parameters
        ---------- 
            self : PlmDCA
                An instance of PlmDCA class
        
        Returns
        -------
            di_scores : np.array 
                 A 1d numpy array of shape (P, ) containing the values of
                direct informations (DI).  P is the total number of unique site pairs.
                Example, index P = 0 contains DI of pair (0, 1),index P = 1 that
                of (0, 2) and so on. The last pair is (L-2, L-1).  Note that the
                direct information is computed from couplings and fields that involve
                residues, although the direct probability is computed for all couplings
                and new fields. The couplings involving a gap are set to 0. The fields
                of gap states are not necessarily zero, they are  the new fields as
                computed by two site model. If Pdir is the direct probabiliy of shape
                (q, q), we use Pdir[:q-1, :q-1] when computing the direct information.
        """

        fields_and_couplings = self.get_fields_and_couplings_from_backend()
        #add attribute to capture the fields and couplings obtained from backend
        self.__fields_and_couplings_all = fields_and_couplings
        couplings = self.get_couplings_no_gap_state(fields_and_couplings)
        reg_fi = self.get_reg_single_site_freqs()
        two_site_model_fields = self.compute_two_site_model_fields(couplings)
        logger.info('\n\tComputing direct information')
        di_scores = msa_numerics.compute_direct_info(
            couplings = couplings, 
            fields_ij = two_site_model_fields,
            reg_fi = reg_fi, 
            seqs_len = self.__seqs_len, 
            num_site_states = self.__num_site_states
        )
        return di_scores 
        
    
    def compute_sorted_DI(self, seqbackmapper=None):
        """Sorted the DI score

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            seqbackmapepr : SequenceBackmapepr
                An instance of SequenceBackmapper class
        
        Returns
        -------
            sorted_di : tuple of site pairs and DCA scores of the form
                (((0, 1), score1), ((0,2), score2 ...)
        """
        di_scores= self.compute_direct_info_unsorted_DI()
        di_scores_dict = dict()
        pair_counter = 0
        for i in range(self.__seqs_len - 1):
            for j in range(i + 1, self.__seqs_len):
                site_pair = (i, j)
                di_scores_dict[site_pair] = di_scores[pair_counter]
                pair_counter += 1
        sorted_di = sorted(di_scores_dict.items(), key =lambda k : k[1], reverse=True)
        if seqbackmapper is not None:
            sorted_di = self.get_mapped_site_pairs_dca_scores(sorted_di, seqbackmapper)
        return sorted_di

    
    def compute_sorted_DI_APC(self, seqbackmapper=None):
        """Performs average product correction of DI scores and sorts them in reverse order.

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class
        
        Returns
        -------
            sorted_DI_apc : list
                A list of tuples the tuples containing (pair, score). The list is
                sorted by score in descending order.
        """
        sorted_DI = self.compute_sorted_DI()
        logger.info('\n\tPerforming average product correction (APC) of DI scores')
        # compute the average score of each site
        av_score_sites = list()
        N = self.__seqs_len 
        for i in range(N):
            i_scores = [score for pair, score in sorted_DI if i in pair]
            assert len(i_scores) == N - 1
            i_scores_sum = sum(i_scores)
            i_scores_ave = i_scores_sum/float(N - 1)
            av_score_sites.append(i_scores_ave)
        # compute average product corrected DI
        av_all_scores = sum(av_score_sites)/float(N)
        sorted_DI_apc = list()
        for pair, score in sorted_DI:
            i, j = pair
            score_apc = score - av_score_sites[i] * (av_score_sites[j]/av_all_scores)
            sorted_DI_apc.append((pair, score_apc))
        # sort the scores as doing APC may have disrupted the ordering
        sorted_DI_apc = sorted(sorted_DI_apc, key = lambda k : k[1], reverse=True)
        if seqbackmapper is not None:
            sorted_DI_apc = self.get_mapped_site_pairs_dca_scores(sorted_DI_apc, seqbackmapper)
        return sorted_DI_apc 
        


if __name__ == '__main__':
    """
    """
    