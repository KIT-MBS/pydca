from pydca.fasta_reader import fasta_reader
from . import scoring_matrix
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import logging
import os

"""Performs sequence back-mapping of a reference sequence to an MSA sequeces.
The back-mapping is done by searching the best matching sequence to the reference.
The searching is carried out by pairwise local alignment of the reference seqeunce
with the sequences in the multiple sequence alignment (MSA). Among best matching
sequence (if there are more than one), the first found is taken. The reference
and best-matching sequences are then pairwise locally aligned to find out the
portion of the the reference subsequence. The residues in this reference
subsequence are then mapped to their counterparts in the MSA columns.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)

class SequenceBackmapper:
    """Defines a sequence backmapper class. Instances of SequenceBackmapper perform
    mapping of a reference sequence to the best matching sequence in the multiple
    sequence alignment (MSA) data.
    """
    def __init__(self, msa_file = None, alignment_data = None, ref_seq=None,
            refseq_file= None, biomolecule=None):
        """Initializes a SequenceBackmapper instance
        Parameters
        -----------
            msa_file : str
                Path to multiple sequence alignment file in FASTA format.
            alignment_data : list
                A list of alignned sequence in integer representation.
            ref_seq : str
                A reference sequences in character form.
            refseq_file: str
                Path to FASTA formatted file containing a reference sequence.
            biomolecule : str
                Type of the sequence data (protein or RNA)

        Returns
        -------
            None : None

        """
        self.__biomolecule = biomolecule.strip().upper()
        if msa_file:
            self.__alignment = fasta_reader.get_alignment_char_form(msa_file,
                biomolecule=self.__biomolecule
            )
        elif alignment_data:
            unique_seqs = []
            for seq in alignment_data:
                if seq not in unique_seqs: unique_seqs.append(seq)
            unique_seqs_char_form = fasta_reader.sequences_to_char_form(
                unique_seqs, self.__biomolecule)
            self.__alignment = unique_seqs_char_form
        else:
            logger.error('\n\tPlease provide alignment file or a list of alignments')
            raise ValueError

        if refseq_file:
            self.__ref_sequence = self._reference_sequence(
                refseq_file=refseq_file)
        elif ref_seq:
            self.__ref_sequence = ref_seq.strip().upper()
        else:
            logger.error('\n\tPlease provide a reference sequence or a FASTA'
                ' file containing a reference sequence')
            raise ValueError
        self._validate_refseq()

        return None


    @property
    def alignment(self):
        """SequenceBackmapper alignment getter
        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class
        Returns
        -------
            self.__alignment : list
                The list of MSA in char represenation
        """
        return self.__alignment


    @property
    def ref_sequence(self):
        """SequenceBackmapper reference sequence getter.

        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns
        -------
            self.__ref_sequence : str
                The reference sequence of a SequenceBackmapper instance.
        """
        return self.__ref_sequence

    def __str__(self):
        """Provides a human readable representation of a SequenceBackmapper
        instance.

        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns
        -------
            describe : str
                Description about a SequenceBackmapper instance.
        """
        describe = '<A sequence backmapper object of biomolecule type {}>'
        return describe.format(self.__biomolecule)


    def _validate_refseq(self):
        """Validate the reference sequence by checking if it contains gaps
        or non-standard residues

        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns
        -------
            None : None
        """
        standard_res_plus_gap = fasta_reader.RES_TO_INT_ALL[self.__biomolecule].keys()
        gap_symbols = ['-', '.', '~']
        standard_residues = [
            state for state in standard_res_plus_gap if state not in gap_symbols
        ]
        for res in self.__ref_sequence:
            if res not in standard_residues:
                logger.error('\n\tReference sequence should only contain standard residues')
                raise ValueError
        return None


    def _reference_sequence(self, refseq_file):
        """Reads a reference sequence from FASTA file. Note that if there are
        more than one sequences in the file, only the one that is found first
        taken as a reference sequence.

        Parameters
        -----------
            refseq_file: str
                Path to fasta file containing the reference sequence.
                If there are multiple sequences, the first is taken

        Returns
        -------
            ref_sequence : str
                Reference sequence
        """
        logger.info('\n\tObtaining reference sequence from file:'
            '\n\t\t{}'.format(refseq_file)
        )
        ref_seqs = fasta_reader.get_alignment_char_form(
            refseq_file,
            biomolecule = self.__biomolecule,
        )
        ref_sequence = ref_seqs[0]
        if len(ref_seqs) > 1 :
            logger.warning('\n\tFound multiple reference sequences in file {}.'
                '\n\tFirst sequence taken as reference'.format(os.path.basename(
                refseq_file)))
        if not ref_sequence:
            logger.error('\n\tNo reference sequence found')
            raise ValueError
        logger.info('\n\tReference sequence:\n\t{}'.format(ref_sequence))

        return ref_sequence.strip().upper()


    def align_pairs_local(self, ref_seq, other_seq, score_only = False):
        """Performs pairwise alignment give two sequences

        Parameters
        ----------
            ref_seq : str
                Reference sequence
            other_seq : str
                Sequence to be aligned to reference
            biomolecule : str
                Sequences type, protein or RNA

        Returns
        -------
            alignments : tuple
                A tuple of pairwise aligned sequences, alignment score and
                start and end indices of alignment
        """
        if self.__biomolecule == 'RNA':
            scoring_mat = scoring_matrix.NUC44
            GAP_OPEN_PEN = -8
            GAP_EXTEND_PEN = 0
        elif self.__biomolecule == 'PROTEIN':
            scoring_mat = blosum62
            GAP_OPEN_PEN = -10
            GAP_EXTEND_PEN = -1
        else:
            logger.error('\n\tUnknown biomolecule type.'
                ' Cannot figure out the scoring matrix.')
            raise ValueError

        alignments = pairwise2.align.localds(
            ref_seq,
            other_seq,
            scoring_mat,
            GAP_OPEN_PEN,
            GAP_EXTEND_PEN,
            score_only = score_only,
        )

        return alignments


    def find_matching_seqs_from_alignment(self):
        """Finds the best matching sequences to the reference
        sequence in the alignment. If multiple matching sequences
        are found, the first one (according to the order in the MSA)
        is taken

        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns
        -------
            best_matching_seqs : list
                A list of best matching sequences to reference
        """

        logger.info('\n\tSearching for sequence(s) that match best with the'
            ' reference sequence')
        # if the first sequence (gaps removed) in MSA matches with reference,
        # return this sequence.
        first_seq_in_alignment = self.__alignment[0]
        first_seq_in_alignment_gaps_removed = first_seq_in_alignment.replace('-','')
        if first_seq_in_alignment_gaps_removed == self.__ref_sequence:
            logger.info('\n\tFirst sequence in alignment (gaps removed) matches reference,'
                '\n\tSkipping regorous search for matching sequence'
            )
            first_seq = list()
            first_seq.append(first_seq_in_alignment)
            return first_seq
        pairwise_scores = []
        for seq_indx, seq in enumerate(self.__alignment):
            seq_gaps_removed = seq.replace('-','')

            score = self.align_pairs_local(
                self.__ref_sequence,
                seq_gaps_removed,
                score_only = True,
                )
            score_at_indx = (seq_indx, score)
            pairwise_scores.append(score_at_indx)

        seq_indx, max_score = max(pairwise_scores, key=lambda x: x[1])
        matching_seqs_indx = [
            indx  for indx, score in pairwise_scores if score == max_score
        ]

        best_matching_seqs = [
            self.__alignment[indx] for indx in matching_seqs_indx
        ]
        num_matching_seqs = len(best_matching_seqs)
        if num_matching_seqs > 1 :
            logger.warning('\n\tFound {} sequences in MSA that match the reference'
                '\n\tThe first sequence is taken as matching'.format(num_matching_seqs)
            )
        return best_matching_seqs

    @staticmethod
    def align_subsequences(ref_middle_subseq = None,
            template_subseq_in_msa = None, num_res_middle_template = None):
        """Aligns the portion of reference sequence (ref_middle_subseq) by
        scanning through the portion of template sequence (template_subseq_in_msa)
        and inserting the gaps encountered into the portion of the reference.

        Parameters
        ----------
            ref_middle_subseq : str
                The portion of the reference sequence that matched when pairwise
                aligned with the template sequence in MSA.
            template_subseq_in_msa : str
                The portion of the template sequence that matched with the
                reference sequence when pairwise aligned with the reference
                sequence.
            num_res_middle_template : int
                The number of residues (excluding gaps) in matched portion of
                the template sequence.

        Returns
        -------
            ''.join(mapped_ref_subseq) : str
                The mapped form of reference sequence. This mapped form contains
                the residues and gaps of the matching portion of the reference
                sequence as well as the newly introduced gap states from the
                tempate sequence as it appears in MSA data.
        """
        mapped_ref_subseq = []
        res_count = 0
        pos = 0
        gap_symbol = '-'.strip()
        for site in template_subseq_in_msa:
            if res_count == num_res_middle_template: break
            if site != gap_symbol:
                mapped_ref_subseq.append(ref_middle_subseq[pos])
                pos += 1
                res_count += 1
                #added after Fabrizio's bug report
                if pos == len(ref_middle_subseq): break
            else:

                if ref_middle_subseq[pos] != gap_symbol:
                    mapped_ref_subseq.append(gap_symbol)
                else:
                    mapped_ref_subseq.append(ref_middle_subseq[pos])
                    pos += 1
        mapped_ref_subseq.extend(list(ref_middle_subseq[pos:]))
        return ''.strip().join(mapped_ref_subseq)


    def map_to_reference_sequence(self):
        """Mapps the reference sequence to the template sequence in alignment
        data. The template sequence is the best scoring sequence when pairwise
        locally aligned with the reference sequence. Here are the steps for
        backmapping:

        i) Find the best matching (template) sequence from the MSA. This is done
        by using pairwise local alignment with the reference sequence.

        ii) Find the aligned portions of the reference and template sequence
        when they are aligned locally. In this step, start  and end indices
        of the matching subsequences, and the number of residues in these two
        subsequences are recorded.

        iii) Map the matching portion of the reference sequence to the matching
        portion of the template sequence in its form in the MSA data. In this
        step, the gaps that are in the template sequence but are not in the
        reference are inserted into the reference sequence at the corresponding
        positions.

        iv) Map the sites (indices) of the reference sequence. In this step,
        mapping is done using the matching start index position of the template
        sequence as it appears in the MSA. The reference sequence residues are
        counted starting from the residue that has been found at the starting
        position of the matching portion when the reference and template
        sequences were locally aligned.

        Parameters
        ----------
            self : SequenceBackmapper
                An instance of SequenceBackmapper class

        Returns
        --------
            mapped_sites : dict
                A dictionary containing the position of mapped residues as they
                appear in the reference sequence as keys and their corresponding
                mapping index in the MSA as values. E.g. {4:9, 5:10, 6:13, ..}
                mapps site 5 (index 4) in the MSA to index 9 (site 10) int the
                reference sequence, ... and so on. The non-aligned residues that might
                appear at the begining or end of the matching subsequences are
                not mapped.
        """
        logger.info('\n\tBackmapping reference sequence to MSA')
        # find best matching sequences to the ref. sequence from the alignment.
        template_sequences_in_msa = self.find_matching_seqs_from_alignment()
        # take the first matching sequence (there can be multiple matches)
        template_seq_in_msa = template_sequences_in_msa[0]
        logger.info('\n\tTemplate sequence in msa:\n{}'.format(template_seq_in_msa))
        # remove the gaps from the matching sequence so that it can be pair
        # aligned with the reference sequence.
        gap_symbol = '-'.strip()
        null_char = ''.strip()
        template_seq_in_msa_gaps_removed = template_seq_in_msa.replace(
            gap_symbol, null_char
        )
        logger.info('\n\tReference sequence and Template sequence (gaps removed)'
            ' respectively:\n{}\n{}'.format(self.__ref_sequence,
                template_seq_in_msa_gaps_removed
            ),
        )
        # pairwise locally align the reference and matching sequenece
        ref_and_template_aligned = self.align_pairs_local(self.__ref_sequence,
            template_seq_in_msa_gaps_removed,
        )
        # caputre the aligned form of the ref. and matching sequences
        # as well as the score, and start and end indices
        ref_seq_aligned = ref_and_template_aligned[0][0]
        template_seq_aligned = ref_and_template_aligned[0][1]
        the_score = ref_and_template_aligned[0][2]
        start_indx = ref_and_template_aligned[0][3]
        end_indx = ref_and_template_aligned[0][4]
        # capture the matching portions of the aligned sequences
        ref_middle_subseq = ref_seq_aligned[start_indx:end_indx]
        template_middle_subseq = template_seq_aligned[start_indx:end_indx]
        logger.info('\n\tMatching subsequences of the reference and the template'
            ' respectively:\n{}\n{}'
            '\n\tMatching start and end positions:[{}, {}]'.format(ref_middle_subseq,
            template_middle_subseq, start_indx + 1, end_indx + 1,
            ),
        )
        # capture the number of residues (excluding gaps) in each subsequence
        num_leading_res_template = len(
            template_seq_aligned[:start_indx].replace(gap_symbol, null_char)
        )
        num_leading_res_ref = len(
            ref_seq_aligned[:start_indx].replace(gap_symbol, null_char)
        )
        num_res_middle_template = len(
            template_middle_subseq.replace(gap_symbol, null_char)
        )
        num_res_middle_ref = len(
            ref_middle_subseq.replace(gap_symbol, null_char)
        )

        # find start index of the matching seq. in the alignment
        res_count = 0
        for k, site in enumerate(template_seq_in_msa, start=0):

            if res_count == num_leading_res_template:
                start_indx_in_msa = k
                break
            if site != gap_symbol: res_count += 1

        template_subseq_in_msa = template_seq_in_msa[start_indx_in_msa:] # excludes only leading sites
        backmapped_ref_subseq = self.align_subsequences(
            ref_middle_subseq = ref_middle_subseq,
            template_subseq_in_msa = template_subseq_in_msa,
            num_res_middle_template = num_res_middle_template,
        )
        logger.info('\nBackmapped ref subsequence:\n{}'.format(backmapped_ref_subseq))
        mapped_sites = dict() # keys are refseq sites and values are matching seq sites in MSA
        mapped_res_count = 0
        for k , site in enumerate(backmapped_ref_subseq, start = 0):
            #we can only map a maximum of alignment length sites
            if k == len(template_seq_in_msa) - start_indx_in_msa: break
            if site != gap_symbol:
                mapped_sites[mapped_res_count + num_leading_res_ref] = start_indx_in_msa + k
                mapped_res_count += 1
        logger.info('\n\tNumber of residues mapped: {}'
            '\n\tNumber of residues in the (original) reference sequence: {}'.format(
                len(mapped_sites), len(self.__ref_sequence))
        )
        # invert the mapping so that keys are matching seq sites in MSA and values are refseq sites.
        mapped_sites_from_msa_to_ref = {
            value:key for key, value in mapped_sites.items()
        }
        return mapped_sites_from_msa_to_ref


if __name__ == '__main__':
    """
    from pydca.config_dca.config_log import LOGGING_CONFIG
    from argparse import ArgumentParser
    import logging
    import logging.config
    logging.config.dictConfig(LOGGING_CONFIG)
    parser = ArgumentParser()
    parser.add_argument('msa_file', help = 'FASTA file containing alignment data')
    parser.add_argument('refseq_file', help = 'FASTA file containing reference sequence')
    parser.add_argument('biomolecule', choices = ['protein', 'PROTEIN', 'rna', 'RNA'])
    args = parser.parse_args()
    alignment_data = fasta_reader.get_alignment_int_form(args.msa_file,
        biomolecule=args.biomolecule)
    seq_backmapper = SequenceBackmapper(alignment_data = alignment_data,
        refseq_file = args.refseq_file , biomolecule = args.biomolecule)
    backmapped_sites =seq_backmapper.map_to_reference_sequence()
    """
    ref_middle_subseq =      'AAAAAAAA---AA' # 8 res
    template_subseq_in_msa = '---BBB--BB-BBBB----' # 9 res
    num_res_middle_template = 9
    seq_backmapper.trigger_gaps(
        ref_middle_subseq = ref_middle_subseq,
        template_subseq_in_msa = template_subseq_in_msa,
        num_res_middle_template = num_res_middle_template,
    )
    """

    """
