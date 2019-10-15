from Bio import AlignIO
from ..sequence_backmapper.sequence_backmapper import SequenceBackmapper
import logging
"""Trims MSA data by gap percentage or removing all gaps corresponding to best
matching sequence to a reference sequence.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)

class MSATrimmerException(Exception):
    """Raises exceptions related to MSA trimming
    """

class MSATrimmer:

    def __init__(self, msa_file, biomolecule=None,max_gap=None, refseq_file=None):
        """
        Parameters
        ----------
            self : MSATrimmer
                An instance of MSATrimmer class
            msa_file : str
                Path to the FASTA formatted MSA file
            biomolecule : str
                Type of biomolecule (protein or RNA)
        """
        self.__msa_file = msa_file
        self.__refseq_file = refseq_file
        self.__max_gap = 0.5 if max_gap is None else max_gap
        if self.__max_gap > 1.0 or self.__max_gap < 0.0:
            logger.error('\n\tThe value of max_gap should be between 0 and 1')
            raise MSATrimmerException
        if biomolecule is not None:
            self.__biomolecule = biomolecule.strip().upper()
        else:
            self.__biomolecule = biomolecule
        self.__alignment_data = list(AlignIO.read(self.__msa_file, 'fasta'))

        logger.info('\n\tMSA file: {0}'
            '\n\tReference sequence file: {1}'
            '\n\tbiomolecule: {2}'
            ''.format(self.__msa_file, self.__refseq_file,
                self.__biomolecule,
            )
        )
        return None


    @property
    def alignment_data(self):
        """
        """
        return self.__alignment_data


    def compute_msa_columns_gap_size(self):
        """Computes the gap size of each column in MSA

        Parameters
        ----------
            self : MSATrimmer
                Instance of MSATrimmer class

        Returns
        -------
            msa_columns_gap_size : tuple
                A tuple of column gap sizes. The column gap size is computed as
                the fraction of gaps in a particular MSA column.

        """
        logger.info('\n\tObtaining columns containing more than {}% of gaps'.format(
            self.__max_gap * 100)
        )
        seqs_len = len(self.__alignment_data[0].seq)
        num_seqs = len(self.__alignment_data)
        logger.info('\n\tTotal number of sequences read from MSA file:{}'
            '\n\tLength of the sequences:{}'.format(num_seqs, seqs_len)
        )
        msa_columns_gap_size = list()
        for i in range(seqs_len):
            num_gaps = 0
            for record in self.__alignment_data:
                state_i = record.seq[i]
                if state_i == '.' or state_i == '-': num_gaps += 1
            gap_fraction_i = float(num_gaps)/float(num_seqs)
            msa_columns_gap_size.append(gap_fraction_i)
        max_gap_size = max(msa_columns_gap_size)
        min_gap_size  = min(msa_columns_gap_size)
        logger.info('\n\tMinimum and maximum gap percentages, respectively:'
            '{0:.2f}% and {1:.2f}%'.format(max_gap_size * 100, min_gap_size * 100)
        )
        return tuple(msa_columns_gap_size)


    def msa_columns_beyond_max_gap(self):
        """Obtains the columns in MSA tha contain more than the given fraction of
        gaps treshold.

        Parameters
        ----------
            self : MSATrimmer
                An instance of MSATrimmer class

        Returns
        -------
            msa_columns_beyond_max_gap : tuple
                A tuple of MSA columns that contain fraction of gaps beyond the
                max_gap
        """
        columns_gap_size = self.compute_msa_columns_gap_size()
        seqs_len = len(self.__alignment_data[0].seq)
        msa_columns_beyond_max_gap = [
            i for i in range(seqs_len) if columns_gap_size[i] > self.__max_gap
        ]
        return tuple(msa_columns_beyond_max_gap)


    def trim_by_gap_size(self):
        """Returns a tuple of MSA columns that have beyond self.__max_gap gap
        fraction.

        Parameters
        ---------
            self : MSATrimmer
                An instance of MSATrimmer class

        Returns
        -------
            columns_to_remove : tuple
                A tuple containing columns that are going to to trimmed. These
                are MSA columns that have a gap fraction beyond self.__max_gap.
        """
        columns_to_remove = self.msa_columns_beyond_max_gap()
        return tuple(columns_to_remove)


    def trim_by_refseq(self, remove_all_gaps=False):
        """Obtains columns in MSA that contain gaps more that the gap treshold
        and do not involve residues in the best matchin sequence with reference.
        If remove_all_gaps is set True, all columns involving gaps in the matching
        sequence to reference are removed.

        Parameters
        ----------
            self : MSATrimmer
                An instance of MSATrimmer
            remove_all_gaps : bool
                If set to True, all columns with gaps in the matching sequence
                with the reference are removed.

        Returns
        -------
            columns_to_remove : tuple
                A tuple of MSA column positions. These columns are going to
                be removed from the MSA.
        """
        seqbackmapper = SequenceBackmapper(msa_file = self.__msa_file,
            refseq_file = self.__refseq_file,
            biomolecule = self.__biomolecule,
        )
        matching_seqs = seqbackmapper.find_matching_seqs_from_alignment()
        logger.info('\n\tRemoving gapped columns corresponding to best'
            ' matching sequence to the reference'
        )
        first_matching_seq = matching_seqs[0]
        logger.info('\n\tSequence in MSA that matches the reference'
            '\n\t{}'.format(first_matching_seq)
        )

        gap_symbols = ['-', '.']
        if not remove_all_gaps:
            candidate_columns_to_remove = self.msa_columns_beyond_max_gap()
            # find out MSA columns that does correspond to gaps w.r.t the sequence
            # in MSA that matches with the reference
            logger.info('\n\tNumber of columns with more than {0:.2f}% gaps:{1}'
                ''.format(self.__max_gap* 100, len(candidate_columns_to_remove))
            )
            columns_to_remove = [
                i for i in candidate_columns_to_remove if first_matching_seq[i] in gap_symbols
            ]
            logger.info('\n\tNumber of columns to remove: {}'.format(len(columns_to_remove)))
        else: # if remove all gaps
            logger.info('\n\tRemoving all columns corresponding to gaps in the matching sequence')
            seqs_len = len(self.__alignment_data[0].seq)
            columns_to_remove = [
                i for i in range(seqs_len) if first_matching_seq[i] in gap_symbols
            ]
            logger.info('\n\tNumber of columns to be removed from MSA:{}'.format(
                len(columns_to_remove))
            )

        return tuple(columns_to_remove)

    
    def get_msa_trimmed_by_refseq(self, remove_all_gaps=False):
        """
        """
        columns_to_remove = self.trim_by_refseq(remove_all_gaps=remove_all_gaps)
        trimmed_msa = list()
        for record in self.__alignment_data:
            seq, seqid = record.seq, record.id
            trimmed_seq = [seq[i] for i in range(len(seq)) if i not in columns_to_remove]
            id_seq_pair = seqid, ''.join(trimmed_seq) 
            trimmed_msa.append(id_seq_pair)
        return trimmed_msa

