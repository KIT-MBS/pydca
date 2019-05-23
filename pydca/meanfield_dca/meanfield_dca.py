from __future__ import absolute_import, division
from . import msa_numerics
from pydca.fasta_reader import fasta_reader
import logging
import numpy as np

"""This module implements Direc Coupling Analysis (DCA) of residue coevolution
for protein and RNA sequences using the mean-field algorithm. The final
coevolution score is computed from the direct probability. The general steps
carried out are outlined as follows

1. Creation of a MeanFieldDCA object
    In this step a MeanFieldDCA object is instantiated with a few attributes.
    The sequences weight is used for computing single and pair site frequncies
    which, in turn are used for the construction of correlation matrix. In addi-
    tion, single site frequencies are used for computing two-site model fields.

2. Construction of correlation matrix
    In this step, the matrix of correlation is computed from regularized pair-site
    and single site frequency counts.

3. Computing couplings
    Once the matrix of correlations is computed, couplings are computed from the
    inverse of the correlation matrix. Note that the couplings are the negative
    of the inverse of the correlation matrix.

4. Computing direct information
    The direct information is computed from the direct probabiliy. To compute
    the direct probablity, new local fields are introduced such that the
    emperical (regularized) single site frequency counts match to the
    corresponding marginal probabilites of the direct probability.

5. The final step is sorting pair of sites according to coevolutionary score
    in descending order. This is dome by  MeanFieldDCA.compute_sorted_DI method.

6. Optionally, one can compute the average product corrected DCA score.  This has
    shown to improve contact prediction in DCA and related methods.


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
    def __init__(
            self, msa_file_name, biomolecule, pseudocount=None,
            sequence_identity=None, force_seq_type=False):
        """MeanFieldDCA object class initializer
        Parameters
        ----------
            msa_file : str
                Name of the FASTA formatted file containing alignmnet
            biomolecule : str
                Type of biomolecule (must be protein or RNA, lower or
                upper case)
            pseudocount : float
                Parameter for regularizing data before DCA analysis.
                Default value is 0.5
            sequence_identity : float
                This parameter's value measure the maximum
                similarity two or more sequences can have so that they can be
                considered distinct, or lumped together otherwise.
            force_seq_type : bool
                Sequences types are typically distinguished by the residue
                representations in the MSA file. Typically, the anticipated
                number of residues plus gap state is 21 for proteins and 5 for
                RNAs. However, it is possible that MSA data can contain
                contain non-standred residues. To avoid confusion about the
                biomolecule type, the default behaviour is to raise error when
                the number of residues plus gap deviates significantly. If the
                user is sure about the biomolecule type is set correctly, they
                can enforce force_seq_type = True option to by pass the error and
                proceed the computation.
        Returns
        -------
            None : None
        """

        self.__pseudocount = pseudocount  if pseudocount is not None else 0.5
        self.__sequence_identity = sequence_identity if sequence_identity is not None else 0.8
        #Validate the value of pseudo count incase user provide an invalid one
        if self.__pseudocount >= 1.0 or self.__pseudocount < 0:
            logger.error('\n\tValue of relative pseudo-count must be'
            ' between 0 and 1.0. Typical value is 0.5')
            raise ValueError
        #Validate the value of sequence identity
        if self.__sequence_identity > 1.0 or self.__sequence_identity <= 0.0:
            logger.error('\n\tValue of sequence-identity must'
            ' not exceed 1 nor less than 0. Typical values are 0.7, 0.8., 0.9')
            raise ValueError
        biomolecule = biomolecule.strip().upper()
        self.__msa_file_name = msa_file_name
        if biomolecule=='RNA'.strip():
            self.__num_site_states = 5
        elif biomolecule=='PROTEIN'.strip():
            self.__num_site_states = 21
        else:
            logger.error(
                '\n\tUnknown biomolecule ... must be protein (PROTEIN) or rna (RNA)',
            )
            raise ValueError
        #verify biomolecule type
        raw_data = fasta_reader.get_alignment_from_fasta_file(self.__msa_file_name)
        num_unique_res = len(set([res for seq in raw_data for res in seq]))
        logger.info('\n\tTotal number of unique site states'
            ' found in alignment data: {}'.format(num_unique_res))
        if num_unique_res < self.__num_site_states or num_unique_res > 2*self.__num_site_states:
            if not force_seq_type:
                logger.error('\n\tThe total number of unique residues plus gap\n\t'
                'found in the alignment file is {}. For {} the anticipated\n\t'
                'number is {}. This error can happen, for example, when\n\t'
                'there are too many non-standared residues in alignment\n\t'
                'data or if protein is entered instead of RNA. To disregared\n\t'
                'the error and continue, Use --force_seq_type command or\n\t'
                'enter correct biomolecule type'.format(num_unique_res,
                    biomolecule, self.__num_site_states))
                raise ValueError
        self.__sequences = fasta_reader.get_alignment_int_form(self.__msa_file_name,
            biomolecule=biomolecule,
            )

        self.__num_sequences = len(self.__sequences)
        self.__sequences_len = len(self.__sequences[0])
        self.__biomolecule = biomolecule
        if self.__sequence_identity < 1.0:
            self.__sequences_weight = self.compute_sequences_weight()
        else :
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
        \ttotal number of (non-effective) sequences: {}
        \teffective number of sequences: {}
        """.format(
            biomolecule,
            self.__num_site_states,
            self.__pseudocount,
            self.__sequence_identity,
            self.__sequences_len,
            self.__num_sequences,
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


    def __call__(self, pseudocount = 0.5 , sequence_identity = 0.8):
        """Resets the value of pseudo count and sequence identity through
        the instance.

        Parameters
        ----------
            self : MeanFieldDCA
                MeanFieldDCA instance.
            pseudocount : float
                The value of the raltive pseudo count. It must be between
                0 and 1. Default value is 0.5.
            sequence_identity : float
                Threshold sequence similarity for computing sequences weight.
                This parameter must be between 0 and 1. Typical values are
                0.7, 0.8, 0.9 or something in between these numbers.

        Returns
        -------
                None : None
        """

        #warn the user that paramertes are being reset
        self.__pseudocount = pseudocount
        self.__sequence_identity = sequence_identity
        logger.warning('\n\tYou have changed one of the parameters (pseudo count or sequence identity)'
        '\n\tfrom their default values'
        '\n\tpseudocount: {} \n\tsequence_identity: {}'.format(
            self.__pseudocount, self.__sequence_identity,
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
            self.__sequence_identity : float
                Cut-off value for sequences similarity above which sequences are
                considered identical
        """

        return self.__sequence_identity


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
            sequence_identity = self.__sequence_identity,
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
            fields :
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


    def compute_hamiltonian(self):
        """Computes the couplings and fields

        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class

        Returns
        -------
            couplings : np.array
                A 2d numpy array of the couplings
            fields : dict
                A dictionary of fields whose keys are sites and values are a numpy
                array of fields.
        """
        reg_fi = self.get_reg_single_site_freqs()
        reg_fij = self.get_reg_pair_site_freqs()
        corr_mat = self.construct_corr_mat(reg_fi, reg_fij)
        logger.info('\n\tComputing the Hamiltonian')
        couplings = self.compute_couplings(corr_mat)
        fields = self.compute_fields(couplings=couplings)
        return fields, couplings


    def compute_sorted_DI(self):
        """Computes direct informations for each pair of sites and sorts them in
        descending order of DCA score.

        Parameters
        ----------
            self : MeanFieldDCA
                The instance.

        Returns
        -------
            sorted_DI : list
                A list of tuples containing site pairs and DCA score, i.e., the
                contents of sorted_DI are [((i, j), score), ...] for all unique
                site pairs (i, j) with corresponding DCA score as computed from
                mean-field DCA. Note that i and j start from 0.

             couplings : np.array
                A 2d numpy array of the couplings.
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

        sorted_DI = dict()
        pair_counter = 0
        for i in range(self.__sequences_len - 1):
            for j in range(i + 1, self.__sequences_len):
                site_pair = (i , j)
                sorted_DI[site_pair] = unsorted_DI[pair_counter]
                pair_counter += 1
        sorted_DI = sorted(sorted_DI.items(), key = lambda k : k[1], reverse=True)
        return sorted_DI, couplings


    def compute_sorted_DI_APC(self):
        """Computes the average DI score for every site.

        Parameters
        ----------
            self : MeanFieldDCA
                An instance of MeanFieldDCA class
        Returns
        -------
            sorted_DI_apc : list
                A list of tuples the tuples containing (pair, score). The list is
                sorted by score in descending order.

            couplings : np.array
                A 2d numpy array of the couplings.

        """

        sorted_DI, couplings = self.compute_sorted_DI()
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
        sorted_DI_apc = list()
        for pair, score in sorted_DI:
            i, j = pair
            score_apc = score - av_score_sites[i] * (av_score_sites[j]/av_all_scores)
            sorted_DI_apc.append((pair, score_apc))
        # sort the scores as doing APC may have disrupted the ordering
        sorted_DI_apc = sorted(sorted_DI_apc, key = lambda k : k[1], reverse=True)
        return sorted_DI_apc, couplings


    def compute_sorted_FN(self):
        """Computes the Frobenius norm of couplings.
        Parameters
        ----------

        Returns
        -------
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
                cij_mean_1 = np.reshape(np.mean(cij, axis=0), (q-1, 1))
                cij_mean_2 = np.reshape(np.mean(cij, axis=1), (1, q-1))
                cij_mean =  np.mean(cij)
                cij_new = cij - cij_mean_1 - cij_mean_2 + cij_mean
                fn_ij = np.sqrt(np.sum(cij_new * cij_new))
                frobenius_norm.append((site_pair, fn_ij))
        frobenius_norm_sorted = sorted(frobenius_norm, key = lambda x : x[1], reverse=True)
        return frobenius_norm_sorted


    def compute_sorted_FN_APC(self):
        """
        """
        raw_fn = self.compute_sorted_FN()
        logger.info('\n\tPerforming average product correction (APC) to Frobenius'
            ' norm of couplings.'
        )
        fn_site_averages = list()
        for i in range(self.__sequences_len):
            scores_i = [score for pair, score in raw_fn if i in pair]
            scores_i_av = np.mean(scores_i)
            fn_site_averages.append(scores_i_av)
        frobenius_norm_apc = list()
        fn_av = np.mean(fn_site_averages)
        for pair, fn_ij in raw_fn:
            fn_i = fn_site_averages[pair[0]]
            fn_j = fn_site_averages[pair[1]]
            fn_ij_apc = fn_ij - (fn_i * fn_j)/fn_av
            frobenius_norm_apc.append(tuple([pair, fn_ij_apc]))
        frobenius_norm_apc = sorted(frobenius_norm_apc, key=lambda x : x[1], reverse=True)
        return frobenius_norm_apc


if __name__ == '__main__':
    """
    import logging
    from argparse import ArgumentParser
    logging.basicConfig(level=logging.DEBUG)
    parser = ArgumentParser(description='reads MSA file name and biomolecule from command line')
    parser.add_argument('msa_file', help='name of alignment FASTA file')
    parser.add_argument('biomolecule', help='biomolecule type', choices=['protein', 'PROTEIN', 'rna', 'RNA'])
    args=parser.parse_args()
    mf_dca = MeanFieldDCA(
        args.msa_file,
        args.biomolecule,
        sequence_identity = 0.8,
        pseudocount = 0.5,
    )

    logger.info('--------------------------HEADER-------------------------------------')
    logger.info('ALIGNMENT FILE: {}'.format(args.msa_file))
    logger.info('SEQUENCE IDENTITY: {}'.format(mf_dca.sequence_identity))
    logger.info('PSEUDO COUNT: {}'.format(mf_dca.pseudocount))
    logger.info('TOTAL NUMBER OF SEQUENCES: {}'.format(mf_dca.num_sequences))
    logger.info('EFFECTIVE NUMBER OF SEQUENCES: {}'.format(mf_dca.effective_num_sequences))
    logger.info('---------------------------------------------------------------------')

    logger.info('\n\tComputing two site model fields')
    sorted_di, couplings = mf_dca.compute_sorted_DI_APC()
    """
