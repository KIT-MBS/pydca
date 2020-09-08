from __future__ import absolute_import, division
from . import msa_numerics
from pydca.fasta_reader import fasta_reader
import logging
import numpy as np

from Bio import Align

"""This module implements Direc Coupling Analysis (DCA) of residue coevolution
for protein and RNA sequences using the mean-field algorithm. The final
coevolution score is computed from the direct probability. The general steps
carried out are outlined as follows

For a detailed information about Direct Coupling Analysis, one can refer to the
following articles:

    a)  Identification of direct residue contacts in protein-protein interaction
        by message-passing
        Martin Weigt, Robert A White, Hendrik Szurmant, James A Hoch, Terence Hwa
        Journal: Proceedings of the National Academy of Sciences
        Volume: 106
        Issue: 1
        Pages: 67-72
    b)  Direct-coupling analysis of residue coevolution captures native contacts
        across many protein families
        Faruck Morcos, Andrea Pagnani, Bryan Lunt, Arianna Bertolino,
        Debora S Marks, Chris Sander, Riccardo Zecchina, Jose N Onuchic,
        Terence Hwa, Martin Weigt
        Journal: Proceedings of the National Academy of Sciences
        Volume: 108
        Issue: 49
        Pages: E1293-E1301

Author(s)  Mehari B. Zerihun, Alexander Schug
"""

logger = logging.getLogger(__name__)

class MeanFieldDCAException(Exception):
    """
    """

class MeanFieldDCA:
    """MeanFieldDCA class. Instances of this class are used to carry out Direct
    Coupling Analysis (DCA) of residue coevolution using the mean-field DCA
    algorithm.
    """
    def __init__(self, msa, biomolecule, pseudocount=None, seqid=None):
        """MeanFieldDCA object class initializer
        Parameters
        ----------
            msa:
                either
                name of a FASTA formatted alignment file
                or
                Bio.Align.MultipleSeqAlignment object

            biomolecule : str
                Type of biomolecule (must be protein or RNA, lower or
                upper case)
            pseudocount : float
                Parameter for regularizing data before DCA analysis.
                Default value is 0.5
            seqid : float
                This parameter's value measure the maximum
                similarity two or more sequences can have so that they can be
                considered distinct, or lumped together otherwise.
        Returns
        -------
            None : None
        """

        self.__pseudocount = pseudocount  if pseudocount is not None else 0.5
        self.__seqid = seqid if seqid is not None else 0.8
        #Validate the value of pseudo count incase user provide an invalid one
        if self.__pseudocount >= 1.0 or self.__pseudocount < 0:
            logger.error('\n\tValue of relative pseudo-count must be'
            ' between 0 and 1.0. Typical value is 0.5')
            raise ValueError
        #Validate the value of sequence identity
        if self.__seqid > 1.0 or self.__seqid <= 0.0:
            logger.error('\n\tValue of sequence-identity must'
            ' not exceed 1 nor less than 0. Typical values are 0.7, 0.8., 0.9')
            raise ValueError
        biomolecule = biomolecule.strip().upper()
        self.__msa = msa
        if biomolecule=='RNA':
            self.__num_site_states = 5
        elif biomolecule=='PROTEIN':
            self.__num_site_states = 21
        else:
            logger.error(
                '\n\tUnknown biomolecule ... must be protein (PROTEIN) or rna (RNA)',
            )
            raise ValueError

        if type(self.__msa) == str:
            self.__sequences = fasta_reader.get_alignment_int_form(
                self.__msa,
                biomolecule=biomolecule,
            )
        elif type(self.__msa) == Align.MultipleSeqAlignment:
            sequences = [record.seq.strip().upper() for record in self.__msa if record.seq]
            self.__sequences = fasta_reader.alignment_letter2int(sequences, biomolecule)
        else:
            raise ValueError("Alignment input parameter is invalid")

        self.__num_sequences = len(self.__sequences)
        self.__sequences_len = len(self.__sequences[0])
        self.__biomolecule = biomolecule
        if self.__seqid < 1.0:
            self.__sequences_weight = self.compute_sequences_weight()
        else :
            # assign each sequence a weight of one
            self.__sequences_weight = np.ones((self.__num_sequences,), dtype = np.float64)
        self.__effective_num_sequences = np.sum(self.__sequences_weight)
        #sometimes users might enter the wrong biomolecule type
        #verify biomolecule type

        mf_dca_info = """\n\tCreated a MeanFieldDCA object with the following attributes
        \tbiomolecule: {}
        \ttotal states at sites: {}
        \tpseudocount: {}
        \tsequence identity: {}
        \talignment length: {}
        \ttotal number of unique sequences (excluding redundant sequences with 100 percent similarity): {}
        \teffective number of sequences (with sequence identity {}): {}
        """.format(
            biomolecule,
            self.__num_site_states,
            self.__pseudocount,
            self.__seqid,
            self.__sequences_len,
            self.__num_sequences,
            self.__seqid,
            self.__effective_num_sequences,
        )
        logger.info(mf_dca_info)
        return None


    def __str__(self):
        """Describes the MeanFieldDCA object.

        Parameters
        ----------
            self: MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            description : str
                A representation about objects created from
                the MeanFieldDCA class.
        """
        description = '<instance of MeanFieldDCA>'
        return description


    def __call__(self, pseudocount = 0.5 , seqid = 0.8):
        """Resets the value of pseudo count and sequence identity through
        the instance.

        Parameters
        ----------
            self : MeanFieldDCA
                MeanFieldDCA instance.
            pseudocount : float
                The value of the raltive pseudo count. It must be between
                0 and 1. Default value is 0.5.
            seqid : float
                Threshold sequence similarity for computing sequences weight.
                This parameter must be between 0 and 1. Typical values are
                0.7, 0.8, 0.9 or something in between these numbers.

        Returns
        -------
                None : None
        """

        #warn the user that paramertes are being reset
        self.__pseudocount = pseudocount
        self.__seqid = seqid
        logger.warning('\n\tYou have changed one of the parameters (pseudo count or sequence identity)'
        '\n\tfrom their default values'
        '\n\tpseudocount: {} \n\tsequence_identity: {}'.format(
            self.__pseudocount, self.__seqid,
            )
        )
        return None


    @property
    def alignment(self):
        """Alignment data getter.
        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        --------
            self.__sequences : list
                A 2d list of alignment sequences in integer representation.
        """

        return self.__sequences

    @property
    def biomolecule(self):
        """Sequence type getter

        Parameters
        ----------
            Self : MeanFieldDCA
                Instance of MeanFieldDCA class
        Returns
        -------
            self.__biomolecule : str
                Biomolecule type (protein or RNA)
        """
        return self.__biomolecule
    @property
    def sequences_len(self):
        """Sequences length getter.

        Parameters
        ---------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__sequences_len : int
                Sequences length in alignment data
        """

        return self.__sequences_len


    @property
    def num_site_states(self):
        """Get number of states for an MSA (eg. 5 for RNAs and 21 for proteins)

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__num_site_states : int
                Maximum number of states in a sequence site
        """

        return self.__num_site_states

    @property
    def num_sequences(self):
        """Getter for the number of sequences read from alignment file

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__num_sequences : int
                The total number of sequences in alignment data
        """

        return self.__num_sequences


    @property
    def sequence_identity(self):
        """Getter for the value of sequence indentity.

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__seqid : float
                Cut-off value for sequences similarity above which sequences are
                considered identical
        """

        return self.__seqid


    @property
    def pseudocount(self):
        """Getter for value of pseudo count

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__pseudocount : float
                Value of pseudo count usef for regularization
        """

        return self.__pseudocount


    @property
    def sequences_weight(self):
        """Getter for the weight of each sequences in alignment data.

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            self.__sequences_weight : np.array(dtype=np.float64)
                A 1d numpy array containing the weight of each sequences in the
                alignment.
        """

        return self.__sequences_weight


    @property
    def effective_num_sequences(self):
        """Getter for the effective number of sequences.

        Parameters
        ----------
            self : MeanFieldDCA
                Instance of MeanFieldDCA class

        Returns
        -------
            np.sum(self.__sequences_weight) : float
                The sum of each sequence's weight.
        """

        return np.sum(self.__sequences_weight)


    def compute_sequences_weight(self):
        """Computes the weight of each sequences in the alignment. If the
        sequences identity is one, each sequences has equal weight and this is
        the maximum weight a sequence in the alignment data can have. Whenever
        the sequence identity is set a value less than one, sequences that have
        similarity beyond the sequence identity are lumped together. If there are
        m similar sequences, their corresponding weight is the reciprocal.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance

        Returns
        -------
            weights : np.array
                A 1d numpy array of size self.__num_sequences containing the
                weight of each sequence.
        """

        logger.info('\n\tComputing sequences weights')
        weights = msa_numerics.compute_sequences_weight(
            alignment_data= np.array(self.__sequences, dtype=np.int32),
            seqid = self.__seqid,
        )
        return weights


    def get_single_site_freqs(self):
        """Computes single site frequency counts.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.

        Returns
        -------
            single_site_freqs : np.array
                A 2d numpy array of shape (L, q) containing the frequency
                count of residues at sequence sites. L is the length of
                sequences in the alignment, and q is the maximum possible
                states a site can accommodate. The last state (q) of each
                site represents a gap.
        """

        logger.info('\n\tComputing single site frequencies')

        single_site_freqs = msa_numerics.compute_single_site_freqs(
            alignment_data = np.array(self.__sequences),
            num_site_states = self.__num_site_states,
            seqs_weight = self.__sequences_weight,
            )
        return single_site_freqs


    def get_reg_single_site_freqs(self):
        """Regularizes single site frequencies.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance

        Returns
        -------
            reg_single_site_freqs : np.array
                A 2d numpy array of shape (L, q) containing regularized single
                site frequencies. L and q are the sequences length and maximum
                number of site-states respectively.
        """

        single_site_freqs = self.get_single_site_freqs()

        logger.info('\n\tRegularizing single site frequencies')

        reg_single_site_freqs = msa_numerics.get_reg_single_site_freqs(
            single_site_freqs = single_site_freqs,
            seqs_len = self.__sequences_len,
            num_site_states = self.__num_site_states,
            pseudocount = self.__pseudocount,
        )
        return reg_single_site_freqs


    def get_pair_site_freqs(self):
        """Computes pair site frequencies

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.

        Returns
        -------
            pair_site_freqs : np.array
                A 2d numpy array of pair site frequncies. It has a shape of
                (N, q-1, q-1) where N is the number of unique site pairs and q
                is the maximum number of states a site can accommodate. Note
                site pairig is performed in the following order: (0, 0), (0, 1),
                ..., (0, L-1), ...(L-1, L) where L is the sequences length. This
                ordering is critical that any computation involding pair site
                frequencies must be implemented in the righ order of pairs.
        """

        logger.info('\n\tComputing pair site frequencies')
        pair_site_freqs = msa_numerics.compute_pair_site_freqs(
        alignment_data = np.array(self.__sequences),
        num_site_states = self.__num_site_states,
        seqs_weight = self.__sequences_weight,
        )
        return pair_site_freqs


    def get_reg_pair_site_freqs(self):
        """Regularizes pair site frequencies

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.

        Returns
        -------
            reg_pair_site_freqs : np.array
                A 3d numpy array of shape (N, q-1, q-1) containing regularized
                pair site frequencies. N is the number of unique site pairs and
                q is the maximum number of states in a sequence site. The
                ordering of pairs follows numbering like (unregularized) pair
                site frequencies.
        """

        pair_site_freqs = self.get_pair_site_freqs()
        logger.info('\n\tRegularizing pair site frequencies')
        reg_pair_site_freqs = msa_numerics.get_reg_pair_site_freqs(
            pair_site_freqs = pair_site_freqs,
            seqs_len = self.__sequences_len,
            num_site_states = self.__num_site_states,
            pseudocount = self.__pseudocount,
        )
        return reg_pair_site_freqs


    def construct_corr_mat(self, reg_fi, reg_fij):
        """Constructs the correlation matrix from regularized frequencies.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.
            reg_fi : np.array
                Regularized single site frequencies.
            reg_fij : np.array
                Regularized pair site frequncies.

        Returns
        -------
            corr_mat : np.array
                A 2d numpy array of (N, N) where N = L*(q-1) where L and q are
                the length of sequences and number of states in a site
                respectively.
        """

        logger.info('\n\tConstructing the correlation matrix')
        corr_mat = msa_numerics.construct_corr_mat(
            reg_fi = reg_fi,
            reg_fij = reg_fij,
            seqs_len = self.__sequences_len,
            num_site_states = self.__num_site_states,
        )
        return corr_mat


    def compute_couplings(self, corr_mat):
        """Computing couplings by inverting the matrix of correlations. Note that
        the couplings are the negative of the inverse of the correlation matrix.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.
            corr_mat : np.array
                The correlation matrix formed from regularized  pair site and
                single site frequencies.

        Returns
        -------
            couplings : np.array
                A 2d numpy array of the same shape as the correlation matrix.
        """

        logger.info('\n\tComputing couplings')
        try:
            couplings = msa_numerics.compute_couplings(corr_mat = corr_mat)
        except Exception as e:
            logger.error('\n\tCorrelation {}\n\tYou set the pseudocount {}.'
                ' You might need to increase it.'.format(e, self.__pseudocount)
            )
            raise
        # capture couplings to avoid recomputing
        self.__couplings = couplings
        logger.info('\n\tMaximum and minimum couplings: {}, {}'.format(
            np.max(couplings), np.min(couplings)))
        return couplings


    def compute_two_site_model_fields(self, couplings, reg_fi):
        """Computes two site model fields by fitting the marginal probabilities
        of the direct probability with the empirical data obtained from the
        alignment

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.
            couplings : np.array
                A 2d numpy array of couplings computed from the correlation matrix.
            reg_fi : np.array
                A 3d numpy array of regularized single site frequencies.
        Returns
        -------
            two_site_model_fields : np.array
                A 3d numpy array of shape (N, q, q) where N is the total number
                of unique site pairs and q is the maximum number of states a site
                can accommodate. The ordering of site pairs is the same as those
                in pair site frequencies.
        """

        logger.info('\n\tComputing two site model fields')
        two_site_model_fields = msa_numerics.compute_two_site_model_fields(
            couplings = couplings,
            reg_fi = reg_fi,
            seqs_len = self.__sequences_len,
            num_site_states = self.__num_site_states,
        )
        return two_site_model_fields


    def compute_fields(self, couplings=None):
        """Computes the local fields of the global probability of sequence space.

        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class

            couplings : np.array
                A 2d numpy array of the couplings. If not give, will be computed.

        Returns
        -------
            fields : dict
                A dictionary of fields whose keys are sites in MSA and whose values
                are arrays of fields per site.
        """

        if couplings is None:
            reg_fi = self.get_reg_single_site_freqs()
            reg_fij = self.get_reg_pair_site_freqs()
            corr_mat = self.construct_corr_mat(reg_fi, reg_fij)
            couplings = self.compute_couplings(corr_mat)
        else:
            reg_fi = self.get_reg_single_site_freqs()
        q = self.__num_site_states
        fields = dict()
        logger.info('\n\tComputing local fields of the global probability function')
        for i in range(self.__sequences_len):
            pi = reg_fi[i]
            piq = pi[-1]
            sum = np.zeros((q-1, 1))
            row_start = i * (q - 1)
            row_end = row_start + (q - 1)
            for j in range(self.__sequences_len):
                if j != i:
                    pj = reg_fi[j]
                    col_start = j * (q - 1)
                    col_end = col_start + (q - 1)
                    couplings_ij = couplings[row_start:row_end, col_start:col_end]
                    pj_col_vec = np.reshape(pj[:-1], (q-1, 1))
                    sum += np.dot(couplings_ij, pj_col_vec)

            fields_i = np.log(pi[:-1]/piq) - np.reshape(sum, (q-1, ))
            fields[i] = fields_i
        return fields


    def shift_couplings(self, couplings_ij):
        """Shifts the couplings value.

        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class
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


    def compute_params(self, seqbackmapper=None, ranked_by=None, linear_dist=None, num_site_pairs=None):
        """Computes fields and couplings with the couplings ranked by DCA score.

        Parameters
        ----------
            self : MeanFieldDCA
                An instanc of MeanFieldDCA class
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
            raise MeanFieldDCAException
        if ranked_by == 'FN': dca_scores = self.compute_sorted_FN(seqbackmapper=seqbackmapper)
        if ranked_by == 'FN_APC': dca_scores = self.compute_sorted_FN_APC(seqbackmapper=seqbackmapper)
        if ranked_by == 'DI': dca_scores = self.compute_sorted_DI(seqbackmapper=seqbackmapper)
        if ranked_by == 'DI_APC': dca_scores = self.compute_sorted_DI_APC(seqbackmapper=seqbackmapper)

        fields = self.compute_fields(couplings=self.__couplings)

        qm1 = self.__num_site_states - 1

        if seqbackmapper is not None:
            # mapping_dict has keys from MSA sites and values from refseq sites
            # we need to reverse this mapping as the fields and couplings are from MSA sites
            mapping_dict = {
                value : key for key, value in self.__refseq_mapping_dict.items()
            }
        else:
            mapping_dict = {
                i : i for i in range(self.__sequences_len)
            }
        # set default number of site pairs whose couplings are to be extracted
        if num_site_pairs is None :
            num_site_pairs = len(seqbackmapper.ref_sequence) if seqbackmapper is not None else len(mapping_dict.keys())
        # we need only the fields corresponding to mapped sites
        fields_mapped = list()
        logger.info('\n\tExtracting fields')
        for i in mapping_dict.keys():
            site_in_msa = mapping_dict[i]
            fields_im = fields[site_in_msa]
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
                    raise MeanFieldDCAException
                row_start = i * qm1
                row_end = row_start + qm1
                column_start = j * qm1
                column_end = column_start + qm1
                couplings_ij = self.__couplings[row_start:row_end, column_start:column_end]
                couplings_ij = self.shift_couplings(couplings_ij) # now couplings_ij is a 2d numpy array
                couplings_ij = np.reshape(couplings_ij, (qm1*qm1,))
                pair_couplings_ij = pair, couplings_ij
                couplings_ranked_by_dca_score.append(pair_couplings_ij)
        if count_pairs < num_site_pairs:
            logger.warning('\n\tObtained couplings for only {} ranked site pairs.'
                '\n\tThis is the maximum number of site paris we can obtain under '
                'the given criteria'.format(count_pairs)
            )

        return tuple(fields_mapped), tuple(couplings_ranked_by_dca_score)


    def  get_mapped_site_pairs_dca_scores(self, sorted_dca_scores, seqbackmapper):
        """Filters mapped site pairs with a reference sequence.

        Parameters
        -----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class
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


    def get_site_pair_di_score(self):
        """Obtains computed direct information (DI) scores from backend and
        puts them a list of tuples of in (site-pair, score) form.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.

        Returns
        -------
            site_pair_di_score : list
                A list of tuples containing site pairs and DCA score, i.e., the
                list [((i, j), score), ...] for all unique ite pairs (i, j)
                such that j > i.
        """
        reg_fi = self.get_reg_single_site_freqs()
        reg_fij = self.get_reg_pair_site_freqs()
        corr_mat = self.construct_corr_mat(reg_fi, reg_fij)
        couplings = self.compute_couplings(corr_mat)
        fields_ij = self.compute_two_site_model_fields(couplings, reg_fi)
        logger.info('\n\tComputing direct information')
        unsorted_DI = msa_numerics.compute_direct_info(
            couplings = couplings,
            fields_ij = fields_ij,
            reg_fi = reg_fi,
            seqs_len = self.__sequences_len,
            num_site_states = self.__num_site_states,
        )

        site_pair_di_score= dict()
        pair_counter = 0
        for i in range(self.__sequences_len - 1):
            for j in range(i + 1, self.__sequences_len):
                site_pair = (i , j)
                site_pair_di_score[site_pair] = unsorted_DI[pair_counter]
                pair_counter += 1
        return site_pair_di_score

    def compute_sorted_DI(self, seqbackmapper=None):
        """Computes direct informations for each pair of sites and sorts them in
        descending order of DCA score.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class.

        Returns
        -------
            sorted_DI : list
                A list of tuples containing site pairs and DCA score, i.e., the
                contents of sorted_DI are [((i, j), score), ...] for all unique
                site pairs (i, j) such that j > i.
        """
        unsorted_DI = self.get_site_pair_di_score()
        sorted_DI = sorted(unsorted_DI.items(), key = lambda k : k[1], reverse=True)
        if seqbackmapper is not None:
            sorted_DI = self.get_mapped_site_pairs_dca_scores(sorted_DI, seqbackmapper)
        return sorted_DI


    def compute_sorted_DI_APC(self, seqbackmapper=None):
        """Computes the average DI score for every site.

        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class.
        Returns
        -------
            sorted_DI_APC : list
                A list of tuples containing site pairs and DCA score, i.e., the
                contents of sorted_DI are [((i, j), score), ...] for all unique
                site pairs (i, j) such that j > i. These DI scores are average
                product corrected.
        """

        sorted_DI = self.compute_sorted_DI() # we must not supply seqbackmapper at this point.
        # the backmapping is done at the end of APC step
        logger.info('\n\tPerforming average product correction (APC) of DI scores')
        # compute the average score of each site
        av_score_sites = list()
        N = self.__sequences_len
        for i in range(N):
            i_scores = [score for pair, score in sorted_DI if i in pair]
            assert len(i_scores) == N - 1
            i_scores_sum = sum(i_scores)
            i_scores_ave = i_scores_sum/float(N - 1)
            av_score_sites.append(i_scores_ave)
        # compute average product corrected DI
        av_all_scores = sum(av_score_sites)/float(N)
        sorted_DI_APC = list()
        for pair, score in sorted_DI:
            i, j = pair
            score_apc = score - av_score_sites[i] * (av_score_sites[j]/av_all_scores)
            sorted_DI_APC.append((pair, score_apc))
        # sort the scores as doing APC may have disrupted the ordering
        sorted_DI_APC = sorted(sorted_DI_APC, key = lambda k : k[1], reverse=True)
        # Now we must do backmapping if seqbackmapper is provided.
        if seqbackmapper is not None:
            sorted_DI_APC = self.get_mapped_site_pairs_dca_scores(sorted_DI_APC, seqbackmapper)
        return sorted_DI_APC


    def compute_sorted_FN(self, seqbackmapper=None):
        """Computes the Frobenius norm of couplings.
        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class.
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class.

        Returns
        -------
            fn_sorted   : list
                A list of tuples containing site pairs and DCA score, i.e., the
                list [((i, j), score), ...] for all unique
                site pairs (i, j) such that j > i.
        """
        reg_fi = self.get_reg_single_site_freqs()
        reg_fij = self.get_reg_pair_site_freqs()
        corr_mat = self.construct_corr_mat(reg_fi, reg_fij)
        couplings = self.compute_couplings(corr_mat)
        logger.info('\n\tComputing Frobenius norm of couplings')
        num_sites = self.__sequences_len
        q = self.__num_site_states
        frobenius_norm = list()
        for i in range(num_sites):
            row_start = i * (q - 1)
            row_end = row_start + (q - 1)
            for j in range(i + 1, num_sites):
                site_pair = (i, j)
                col_start = j * (q - 1)
                col_end = col_start + (q - 1)
                cij = couplings[row_start:row_end, col_start:col_end]
                cij_mean_1 = np.reshape(np.mean(cij, axis=0), (1, q-1))
                cij_mean_2 = np.reshape(np.mean(cij, axis=1), (q-1, 1))
                cij_mean =  np.mean(cij)
                cij_new = cij - cij_mean_1 - cij_mean_2 + cij_mean
                fn_ij = np.sqrt(np.sum(cij_new * cij_new))
                frobenius_norm.append((site_pair, fn_ij))
        fn_sorted = sorted(frobenius_norm, key = lambda x : x[1], reverse=True)
        if seqbackmapper is not None:
            fn_sorted = self.get_mapped_site_pairs_dca_scores(fn_sorted, seqbackmapper)
        return fn_sorted


    def compute_sorted_FN_APC(self, seqbackmapper = None):
        """Performs average product correction (APC) on DCA scores

        Parameters
        ----------
            self    : MeanFieldDCA
                An instance of MeanFieldDCA class.
            seqbackmapper : SequenceBackmapper
                An instance of SequenceBackmapper class.

        Returns
        -------
            sorted_FN_APC : list
                A list of tuples containing site pairs and DCA score, i.e., the
                list [((i, j), score), ...] for all unique site pairs (i, j)
                such that j > i. The DCA scores are average product corrected.
        """
        raw_FN = self.compute_sorted_FN() # Must not supply seqbackmapper at this stage.
        logger.info('\n\tPerforming average product correction (APC) to Frobenius'
            ' norm of couplings.'
        )

        # compute the average score of each site
        av_score_sites = list()
        N = self.__sequences_len
        for i in range(N):
            i_scores = [score for pair, score in raw_FN if i in pair]
            assert len(i_scores) == N - 1
            i_scores_sum = sum(i_scores)
            i_scores_ave = i_scores_sum/float(N - 1)
            av_score_sites.append(i_scores_ave)
        # compute average product corrected DI
        av_all_scores = sum(av_score_sites)/float(N)
        sorted_FN_APC = list()
        for pair, score in raw_FN:
            i, j = pair
            score_apc = score - av_score_sites[i] * (av_score_sites[j]/av_all_scores)
            sorted_FN_APC.append((pair, score_apc))
        sorted_FN_APC = sorted(sorted_FN_APC, key=lambda x : x[1], reverse=True)
        # Must do backmapping is sebackmapper is not None
        if seqbackmapper is not None:
            sorted_FN_APC = self.get_mapped_site_pairs_dca_scores(sorted_FN_APC, seqbackmapper)
        return sorted_FN_APC


if __name__ == '__main__':
    """
    """
