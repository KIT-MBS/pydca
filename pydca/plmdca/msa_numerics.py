import numpy as np 
from numba import jit 
from numba import prange as parallel_range

"""Computes the direct information (DI) for pseudolikelihood maximization direct
coupling analysis.

Author: Mehari B. Zerihun
"""


@jit(nopython=True, parallel=True)
def compute_sequences_weight(alignment_data=None, sequence_identity=None):
    """Computes weight of sequences. The weights are calculated by lumping
    together sequences whose identity is greater that a particular threshold.
    For example, if there are m similar sequences, each of them will be assigned
    a weight of 1/m. Note that the effective number of sequences is the sum of
    these weights.

    Parameters
    ----------
        alignmnet_data : np.array()
            Numpy 2d array of the alignment data, after the alignment is put in
            integer representation
        sequence_identity : float
            Value at which beyond this sequences are considered similar. Typical
            values could be 0.7, 0.8, 0.9 and so on

    Returns
    -------
        seqs_weight : np.array()
            A 1d numpy array containing computed weights. This array has a size
            of the number of sequences in the alignment data.
    """
    alignment_shape = alignment_data.shape
    num_seqs = alignment_shape[0]
    seqs_len = alignment_shape[1]
    seqs_weight = np.zeros((num_seqs,), dtype=np.float64)
    #count similar sequences
    for i in parallel_range(num_seqs):
        seq_i = alignment_data[i]
        for j in range(num_seqs):
            seq_j = alignment_data[j]
            iid = np.sum(seq_i==seq_j)
            if np.float64(iid)/np.float64(seqs_len) > sequence_identity:
                seqs_weight[i] += 1
    #compute the weight of each sequence in the alignment
    for i in range(num_seqs): seqs_weight[i] = 1.0/float(seqs_weight[i])
    return seqs_weight


@jit(nopython=True)
def compute_single_site_freqs(alignment_data=None,
        num_site_states=None, seqs_weight=None):
    """Computes single site frequency counts for a particular aligmnet data.

    Parameters
    ----------
        alignment_data : np.array()
            A 2d numpy array of alignment data represented in integer form.

        num_site_states : int
            An integer value fo the number of states a sequence site can have
            including a gap state. Typical value is 5 for RNAs and 21 for
            proteins.

        seqs_weight : np.array()
            A 1d numpy array of sequences weight

    Returns
    -------
        single_site_freqs : np.array()
            A 2d numpy array of of data type float64. The shape of this array is
            (seqs_len, num_site_states) where seqs_len is the length of sequences
            in the alignment data.
    """
    alignment_shape = alignment_data.shape
    #num_seqs = alignment_shape[0]
    seqs_len = alignment_shape[1]
    m_eff = np.sum(seqs_weight)
    single_site_freqs = np.zeros(shape = (seqs_len, num_site_states),
        dtype = np.float64)
    for i in range(seqs_len):
        for a in range(1, num_site_states + 1):#we need gap states single site freqs too
            column_i = alignment_data[:,i]
            freq_ia = np.sum((column_i==a)*seqs_weight)
            single_site_freqs[i, a-1] = freq_ia/m_eff
    return single_site_freqs


@jit(nopython=True)
def get_reg_single_site_freqs(single_site_freqs = None, seqs_len = None,
        num_site_states = None, pseudocount = None):
    """Regularizes single site frequencies.

    Parameters
    ----------
        single_site_freqs : np.array()
            A 2d numpy array of single site frequencies of shape
            (seqs_len, num_site_states). Note that gap state frequencies are
            included in this data.
        seqs_len : int
            The length of sequences in the alignment data
        num_site_states : int
            Total number of states that a site in a sequence can accommodate. It
            includes gap states.
        pseudocount : float
            This is the value of the relative pseudo count of type float.
            theta = lambda/(meff + lambda), where meff is the effective number of
            sequences and lambda is the real pseudo count.

    Returns
    -------
        reg_single_site_freqs : np.array()
            A 2d numpy array of shape (seqs_len, num_site_states) of single site
            frequencies after they are regularized.
    """
    reg_single_site_freqs = single_site_freqs
    theta_by_q = np.float64(pseudocount)/np.float64(num_site_states)
    for i in range(seqs_len):
        for a in range(num_site_states):
            reg_single_site_freqs[i, a] = theta_by_q + \
                (1.0 - pseudocount)*reg_single_site_freqs[i, a]
    return reg_single_site_freqs


@jit(nopython=True)
def slice_couplings(couplings = None, site_pair=None, num_site_states=None, seqs_len = None):
    """Constructs couplings array suitable for computing two-site-model fields as well 
    as DI scores. 

    Parameters
    ----------
        couplings : np.array
            A 1 array of couplings excluding gap state couplings
        site_pair : tuple
            Site pair (i, j) such that j > i with o <= i < seqs_len
        num_site_states : int 
            Number of site states for sequence 
        seqs_len : int 
            Length of sequences in MSA data 
    """
    i, j = site_pair[0], site_pair[1]
    q = num_site_states
    qm1 = q - 1
    pair_loc = int((seqs_len * (seqs_len - 1)/2) - (seqs_len - i) * ((seqs_len - i) - 1)/2  + j  - i - 1)
    start_indx = pair_loc * qm1 * qm1  
    end_indx = start_indx + qm1 * qm1
    couplings_ij = np.zeros((q, q), dtype = np.float64)
    couplings_tmp = couplings[start_indx:end_indx]
    couplings_ij[:q-1, :q-1]  = np.reshape(couplings_tmp, shape = (q-1, q-1))
    return couplings_ij 


@jit(nopython=True)
def compute_two_site_model_fields(couplings = None, reg_fi = None,
        seqs_len = None, num_site_states = None):
    """Computes two-site model fields iteratively.

    Parameters
    ----------
        couplings : np.array
            A numpy array of couplings of shape (N, N) where
            N = seqs_len * (num_site_states - 1)

        reg_fi : np.array
            A numpy array of regularized single site frequncies of shape
            (seqs_len, num_site_states)

        seqs_len : int
            Length of sequences in alignment data

        num_site_states : int
            Total number of states a site in a sequence can accommodate,
            including gap state.

    Returns
    -------
        two_site_model_fields : np.array
            A numpy array of shape (P, 2, num_site_states), where P is the number
            of unique site pairs excluding self pairings.
            P = seqs_len * (seqs_len - 1)/2.
    """
    num_unique_pairs = seqs_len * (seqs_len -1)
    num_unique_pairs /= 2
    q = num_site_states
    two_site_model_fields = np.zeros((np.int64(num_unique_pairs), 2, q), dtype=np.float64)
    TOLERANCE = 1.0e-4
    pair_counter = 0
    for i in range(seqs_len - 1):
        freq_i = np.reshape(reg_fi[i], (q, 1))
        for j in range(i + 1, seqs_len):
            site_pair = (i, j)
            freq_j = np.reshape(reg_fi[j], (q, 1))
            couplings_ij = np.exp(slice_couplings(couplings = couplings,
                site_pair = site_pair, num_site_states = q, seqs_len=seqs_len)
            )
            fields_i_old = np.full((q, 1), 1.0/np.float64(q))
            fields_j_old = np.full((q, 1), 1.0/np.float64(q))
            max_fields_change = 10.0
            while max_fields_change > TOLERANCE:
                x_i = np.dot(couplings_ij , fields_j_old)
                x_j = np.dot(np.transpose(couplings_ij), fields_i_old)

                fields_i_new =  freq_i / x_i
                fields_i_new /= np.sum(fields_i_new)
                fields_j_new = freq_j / x_j
                fields_j_new /= np.sum(fields_j_new)

                delta_fields_i = np.max(np.absolute(fields_i_new - fields_i_old))
                delta_fields_j = np.max(np.absolute(fields_j_new - fields_j_old))
                max_fields_change = np.max(np.array([delta_fields_i, delta_fields_j]))

                fields_i_old = fields_i_new
                fields_j_old = fields_j_new
            #capture computed fields after iteration is converged
            two_site_model_fields[pair_counter][0] = fields_i_new.T
            two_site_model_fields[pair_counter][1] = fields_j_new.T
            pair_counter += 1
    return two_site_model_fields


@jit(nopython=True)
def compute_direct_info(couplings = None, fields_ij = None, reg_fi = None,
        seqs_len = None, num_site_states = None):
    """Computes the direct information from direct probabilities.

    Parameters
    ----------
        couplings : np.array
            A 2d numpy array of shape (L(q-1), L(q-1)), where L and q are the
            length of sequences in MSA and number of site-states respectively.
            Note that the couplings are the negative of the inverse of the
            correlation matrix.

        fields_ij : np.array
            A 3d numpy array of two-site model fields. The shape of this array
            is (P, 2, q). Where P is the number of unique site pairs and q is the
            total number of site states. The ordering of site-pairs is very
            important. For example index P=0 refers to site pairs (0, 1), and
            as p increase the pairs are (0, 2), ... ,(0, L-1), (1, 2), ...,
            (1, L-1), ..., (L-2, L-1). the first index of the second dimension
            refers to the first site in site pair. Example, fields_ij[0][0]
            contains the fields of site 0 when its paired with site 1, and
            fields_ij[0][1] contains those of site 1 in the same pair, and so on.

        reg_fi : np.array
            A 2d numpy array of regularized single site frequencies. It has
            a shape of (L, q) where L and q are the length of the sequences
            in alignment data and number of total site states respectively.
            Example, reg_fi[0] contains the frequencies of the first column in
            MSA.

        seqs_len : int
            The length of sequences in MSA.

        num_site_states : int
            The total number of residues plus gap.

    Returns
    -------
        unsorted_DI : np.array
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
    num_unique_pairs = np.int64(seqs_len * (seqs_len - 1)/2)
    unsorted_DI = np.zeros(num_unique_pairs, dtype=np.float64)
    q = num_site_states
    EPSILON = 1.0e-20
    pair_counter = 0
    for i in range(seqs_len - 1):
        fi = reg_fi[i]
        for j in range(i + 1, seqs_len):
            site_pair = (i, j)
            fj = reg_fi[j]
            #h_i = fields_ij[pair_counter][0]
            #h_j = fields_ij[pair_counter][1]
            hij = np.dot(np.reshape(fields_ij[pair_counter][0], (q, 1)),
                np.transpose(np.reshape(fields_ij[pair_counter][1], (q, 1))),
            )

            couplingsij = np.exp(slice_couplings(couplings = couplings,
                site_pair = site_pair, num_site_states = q, seqs_len = seqs_len)
            )
            #Compute direct information
            pdir_ij = couplingsij * hij
            pdir_ij /= np.sum(pdir_ij)
            #Compute product of single site frequencies
            fij = np.dot(np.reshape(fi, (q, 1)),
                np.transpose(np.reshape(fj, (q, 1)))
            )
            #Only take into account residue residue interactions for computing
            #direct information
            fij_residues = fij[:q-1, :q-1] + EPSILON # + operator creats a copy
            pdir_ij_residues = pdir_ij[:q-1, :q-1] + EPSILON
            pdir_by_fij_residues =  pdir_ij_residues/fij_residues
            #Compute direct information
            DI_ij = np.sum(pdir_ij_residues * np.log(pdir_by_fij_residues))
            unsorted_DI[pair_counter] = DI_ij
            #Move to the next site pair
            pair_counter += 1

    return unsorted_DI

