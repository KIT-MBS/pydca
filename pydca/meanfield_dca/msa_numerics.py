
import numpy as np
from numba import jit
from numba import prange as parallel_range


"""This module implements computationally constly routines while performing
Direct Coupling Analysis.

Author : Mehari B. Zerihun
"""

@jit(nopython=True, parallel=True)
def compute_sequences_weight(alignment_data=None, seqid=None):
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
        seqid : float
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
            if np.float64(iid)/np.float64(seqs_len) > seqid:
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


# This function is replaced by the parallelized version below
@jit(nopython=True)
def compute_pair_site_freqs_serial(alignment_data=None, num_site_states=None,
        seqs_weight=None):
    
    """Computes pair site frequencies for an alignmnet data.

    Parameters
    ----------
        alignment_data : np.array()
            A 2d numpy array conatining alignment data. The residues in the
            alignment are in integer representation.
        num_site_states : int
            The number of possible states including gap state that sequence
            sites can accomodate. It must be an integer
        seqs_weight:
            A 1d numpy array of sequences weight

    Returns
    -------
        pair_site_freqs : np.array()
            A 3d numpy array of shape
            (num_pairs, num_site_states, num_site_states) where num_pairs is
            the number of unique pairs we can form from sequence sites. The
            pairs are assumed to in the order (0, 1), (0, 2) (0, 3), ...(0, L-1),
            ... (L-1, L). This ordering is critical and any change must be
            documented.
    """
    alignment_shape = alignment_data.shape
    num_seqs = alignment_shape[0]
    seqs_len = alignment_shape[1]
    num_site_pairs = (seqs_len -1)*seqs_len/2
    num_site_pairs = np.int64(num_site_pairs)
    m_eff = np.sum(seqs_weight)
    pair_site_freqs = np.zeros(
        shape=(num_site_pairs, num_site_states - 1, num_site_states - 1),
        dtype = np.float64)
    pair_counter = 0
    for i in range(seqs_len-1):
        column_i = alignment_data[:, i]
        for j in range(i+1, seqs_len):
            column_j = alignment_data[:, j]
            for a in range(1, num_site_states):
                count_ai = column_i==a
                for b in range(1, num_site_states):
                    count_bj = column_j==b
                    count_ai_bj = count_ai * count_bj
                    freq_ia_jb = np.sum(count_ai_bj*seqs_weight)
                    pair_site_freqs[pair_counter, a-1, b-1] = freq_ia_jb/m_eff
            #move to the next site pair (i, j)
            pair_counter += 1
    return pair_site_freqs


@jit(nopython=True, parallel=True)
def compute_pair_site_freqs(alignment_data=None, num_site_states=None, seqs_weight=None):
    """Computes pair site frequencies for an alignmnet data.

    Parameters
    ----------
        alignment_data : np.array()
            A 2d numpy array conatining alignment data. The residues in the
            alignment are in integer representation.
        num_site_states : int
            The number of possible states including gap state that sequence
            sites can accomodate. It must be an integer
        seqs_weight:
            A 1d numpy array of sequences weight

    Returns
    -------
        pair_site_freqs : np.array()
            A 3d numpy array of shape
            (num_pairs, num_site_states, num_site_states) where num_pairs is
            the number of unique pairs we can form from sequence sites. The
            pairs are assumed to in the order (0, 1), (0, 2) (0, 3), ...(0, L-1),
            ... (L-1, L). This ordering is critical and any change must be
            documented.
    """
    alignment_shape = alignment_data.shape
    num_seqs = alignment_shape[0]
    seqs_len = alignment_shape[1]
    num_site_pairs = (seqs_len -1)*seqs_len/2
    num_site_pairs = np.int64(num_site_pairs)
    m_eff = np.sum(seqs_weight)
    pair_site_freqs = np.zeros(
        shape=(num_site_pairs, num_site_states - 1, num_site_states - 1),
        dtype = np.float64
    )
    for i in parallel_range(seqs_len - 1):
        column_i = alignment_data[:, i]
        for j in range(i+1, seqs_len):
            pair_site = int((seqs_len * (seqs_len - 1)/2) - (seqs_len - i) * ((seqs_len - i) - 1)/2  + j  - i - 1)
            column_j = alignment_data[:, j]
            for a in range(1, num_site_states):
                count_ai = column_i==a
                for b in range(1, num_site_states):
                    count_bj = column_j==b
                    count_ai_bj = count_ai * count_bj
                    freq_ia_jb = np.sum(count_ai_bj*seqs_weight)
                    pair_site_freqs[pair_site, a-1, b-1] += freq_ia_jb/m_eff
    return pair_site_freqs 

@jit(nopython=True)
def get_reg_pair_site_freqs(pair_site_freqs = None, seqs_len = None,
        num_site_states = None, pseudocount = None):
    """Regularizes pair site frequencies

    Parameters
    ----------
        pair_site_freqs : np.array()
            A 3d numpy array of shape (num_unique_site_pairs, num_site_states -1,
            num_site_states -1) containing raw pair site frequency counts where
            num_unique_site_pairs is the total number of unique site pairs
            excluding self pairing. Note that the order in with the pairing is
            done is important. It must be taken in (0, 1), (0,2), ...,
            (0, seqs_len-1), (1, 2)... order. Note that this data does not
            contain pairings with gap states.
        seqs_len : int
            The length of sequences in the alignment.
        num_site_states : int
            The total number of states that a site in the sequences can
            accommodate. This includes gap states.

    Returns
    -------
        reg_pair_site_freqs : np.array()
            A numpy array of shape the same as pair_site_freqs
    """
    reg_pair_site_freqs = pair_site_freqs
    theta_by_qsqrd = pseudocount/float(num_site_states * num_site_states)
    pair_counter = 0
    for i in range(seqs_len - 1):
        for j in range(i + 1, seqs_len):
            for a in range(num_site_states-1):
                for b in range(num_site_states-1):
                    reg_pair_site_freqs[pair_counter, a, b] = theta_by_qsqrd + \
                        (1.0 - pseudocount)*reg_pair_site_freqs[pair_counter, a, b]
            pair_counter += 1
    return reg_pair_site_freqs


@jit(nopython=True)
def construct_corr_mat(reg_fi = None, reg_fij = None, seqs_len = None,
        num_site_states = None):
    """Constructs correlation matrix from regularized frequency counts.

    Parameters
    ----------
        reg_fi : np.array()
            A 2d numpy array of shape (seqs_len, num_site_states) of regularized
            single site frequncies. Note that only fi[:, 0:num_site_states-1] are
            used for construction of the correlation matrix, since values
            corresponding to fi[:, num_site_states]  are the frequncies of gap
            states.
        reg_fij : np.array()
            A 3d numpy array of shape (num_unique_pairs, num_site_states -1,
            num_site_states - 1), where num_unique_pairs is the total number of
            unique site pairs execluding self-pairings.
        seqs_len : int
            The length of sequences in the alignment
        num_site_states : int
            Total number of states a site in a sequence can accommodate.

    Returns
    -------
        corr_mat : np.array()
            A 2d numpy array of shape (N, N)
            where N = seqs_len * num_site_states -1
    """
    corr_mat_len = seqs_len * (num_site_states - 1)
    corr_mat = np.zeros((corr_mat_len, corr_mat_len), dtype=np.float64)
    pair_counter = 0
    for i in range(seqs_len):
        site_i = i * (num_site_states - 1)
        for j in range(i, seqs_len):
            site_j = j * (num_site_states - 1)
            for a in range(num_site_states - 1):
                row = site_i + a
                for b in range(num_site_states -1):
                    col = site_j + b
                    if i==j:
                        fia, fib = reg_fi[i, a], reg_fi[i, b]
                        corr_ij_ab = fia*(1.0 - fia) if a == b else -1.0*fia*fib
                    else:
                        corr_ij_ab = reg_fij[pair_counter, a, b] - reg_fi[i, a] * reg_fi[j, b]
                    corr_mat[row, col] = corr_ij_ab
                    corr_mat[col, row] = corr_ij_ab
            if i != j: pair_counter += 1

    return corr_mat


@jit(nopython=True)
def compute_couplings(corr_mat = None):
    """Computes the couplings by inverting the correlation matrix

    Parameters
    ----------
        corr_mat : np.array()
            A numpy array of shape (N, N) where N = seqs_len *(num_site_states -1)
            where seqs_len  is the length of sequences in the alignment data and
            num_site_states is the total number of states a site in a sequence
            can accommodate, including gapped states.

    Returns
    -------
        couplings : np.array()
            A 2d numpy array of the same shape as the correlation matrix. Note
            that the couplings are the negative of the inverse of the
            correlation matrix.
    """
    couplings = np.linalg.inv(corr_mat)
    couplings = -1.0*couplings
    return couplings


@jit(nopython=True)
def slice_couplings(couplings = None, site_pair = None, num_site_states=None):
    """Returns couplings corresponding to site pair (i, j). Note that the
    the couplings involving gaps are included, but they are set to zero.

    Parameters
    ----------
        couplings : np.array
            A 2d numpy array of couplings. It has a shape of (L(q-1), L(q-1))
            where L and q are the length of sequences in alignment data and total
            number of standard residues plus gap.
        site_pair : tuple
            A tuple of site pairs. Example (0, 1), (0, L-1), ..., (L-2, L-1).
        num_site_states : int
            The value of q.

    Returns
    -------
        couplings_ij : np.array
            A2d numpy array of shape (q, q) containing the couplings. Note that
            couplings_ij[q, :] and couplings[:, q] are set to zero.
    """
    q = num_site_states
    couplings_ij = np.zeros((q, q), dtype = np.float64)
    row_begin = site_pair[0] * (q - 1)
    row_end = row_begin + q - 1
    column_begin = site_pair[1] * (q -1)
    column_end = column_begin + q - 1
    couplings_ij[:q-1, :q-1] = couplings[row_begin:row_end, column_begin:column_end]
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
                site_pair = site_pair, num_site_states = q))
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
                site_pair = site_pair, num_site_states = q),
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
