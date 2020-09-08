import logging
from Bio import AlignIO

"""Reads alignment data from FASTA files, can convert residue representation in
in the sequences for char to int or vice versa.

Author: Mehari B. Zerihun
"""

__all__ = [
    'get_alignment_from_fasta_file',
    'get_alignment_int_form',
    'get_alignment_char_form',
    'sequences_to_char_form',
]

logger = logging.getLogger(__name__)

#Protein residues  (numbering is mostly random, not alphabetical)
#==================================================================
#1. Alanine - ala - A, 2. Arginine - arg - R, 3. Asparagine - asn - N
#4. Aspartic acid - asp - D, 5. Cysteine -cys - C, 6. Glutamine - gln - Q
#7. Glutamic acid - glu - E, 8. Glycine - gly - G, 9. Histidine - his - H
#10. Isoleucine - ile - I, 11. Leucine - lue -L, 12. Lysine - lys - K
#13. Methionine - met - M, 14. Phenylalanine - phe - F, 15. Proline - pro - P
#16. Serine - ser - S, 17. Threonine - thr - T, 18. Tryptophan - trp - W
#19. Tyrosine - tyr - Y, 20. Valine - val - V

#RNA nucleotides
#====================================
#1. Adenine - A, 2. Cytosine - C
#3. Guanine - G, 4. Uracil - U

RES_TO_INT_ALL = {
    'PROTEIN':{
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W':19, 'Y':20,
        '-':21, '.':21, '~':21,
    },
    'RNA':{
        'A':1, 'C':2, 'G':3, 'U':4, '-':5, '.':5, '~':5,
    },
}


class FastaReaderError(Exception):
    """Raise exceptions related to reading alignment data
    """


def res_to_char(biomolecule):
    """Creates a mapping for residues from char to int to
    int to char.

    Parameters
    ----------
        biomolecule : str
            Type of biomolecule residues represent
            (must be protein or rna)

    Returns
    -------
        RES_TO_CHAR : dict
            A dictionary of residues with the keys in
            int representation and the values in char representation
    """
    biomolecule = biomolecule.strip().upper()
    RES_TO_INT = RES_TO_INT_ALL[biomolecule]
    EXCLUDE = ['.', '~']

    RES_TO_CHAR = {
        val:key for key, val in RES_TO_INT.items() if key not in EXCLUDE
    }
    return RES_TO_CHAR


def get_alignment_from_fasta_file(file_name):
    """Read sequences from FASTA file using Bio.AlignIO.read()

    Parameters
    ----------
        file_name : str
            Path to FASTA formatted file.

    Returns
    -------
        alignment : list
            A list of biomolecular sequence strings.
    """
    alignment = []
    try:
        record_iterator = AlignIO.read(file_name, 'fasta')
        #biopython just reads the records if there are tags (>some key).
        #It doesn't know if the file is really a biological sequence or not
    except Exception as expt:
        error_msg='\n\tError occured while reading from fasta file: {}.' +\
            '\n\tError type:{}\n\tArguments:{!r}'
        logger.error(error_msg.format(file_name, type(expt).__name__, expt.args))
        raise
    else:
        if any(True for _ in record_iterator):
            for record in record_iterator:
                seq = record.seq.strip()
                if seq: alignment.append(seq.upper())
            if not alignment:
                logger.error(
                    '\n\trecord_iterator returned by AlignIO.read()'
                    ' has no sequences',
                )
                raise ValueError

        else:
            logger.error(
                '\n\trecord_iterator returned by AlignIO.read() is empty',
            )
            raise ValueError
    return alignment


def alignment_letter2int(alignment, biomolecule='protein'):
    """
    Converts sequences in a multiple sequence alignment from one letter to integer representation.
    """
    biomolecule=biomolecule.strip().upper()
    if biomolecule not in ['PROTEIN','RNA']:
        logger.error(
            '\n\t{} entered. Biomolecule must be either PROTEIN or RNA'.format(
                biomolecule))
        raise ValueError
    NUM_SITE_STATES = 21 if biomolecule == 'PROTEIN' else 5
    RES_TO_INT = RES_TO_INT_ALL[biomolecule]
    alignment_int_form = []

    num_seqs_with_non_standard_res = 0
    num_non_standard_res = 0
    total_num_seqs_in_msa = 0
    for seq in alignment:
        try:
            seq_int = [RES_TO_INT[res.upper()] for res in seq]
        except KeyError:
            num_seqs_with_non_standard_res += 1
            seq_int = []
            for res in seq:
                res = res.upper()
                if res in RES_TO_INT.keys():
                    seq_int.append(RES_TO_INT[res.upper()])
                else:
                    num_non_standard_res += 1
                    seq_int.append(NUM_SITE_STATES)
        total_num_seqs_in_msa += 1
        if seq_int not in alignment_int_form:
            alignment_int_form.append(seq_int)
    if num_seqs_with_non_standard_res > 0:
        logger.info('\n\tFound {} non-standard residues in {} sequences'
            ''.format(num_non_standard_res, num_seqs_with_non_standard_res)
        )
    logger.info('\n\tTotal number of sequences read from file: {}'.format(total_num_seqs_in_msa))
    if not alignment_int_form:
        logger.error('\n\tNo data found in alignment in integer representation')
        raise ValueError
    return alignment_int_form


def get_alignment_int_form(file_name, biomolecule='protein'):
    """Converts sequences in integer representation. The sequences are
    first read by get_alignment_from_fasta_file(file_name) function that returns
    a list of sequences.

    Parameters
    ----------
        file_name : str
            Fasta file name containing the alignment data.
        biomolecule : str
            The type of biomolecule the sequence data reprsents. This can be
            either portein or RNA in lower or upper cases.

    Returns
    -------
      alignment_int_form : list
        a list of alignments, each sequence in a list of integers.
    """

    alignment = get_alignment_from_fasta_file(file_name)
    alignment_int_form = alignment_letter2int(alignment, biomolecule)

    return alignment_int_form


def get_alignment_char_form(file_name, biomolecule = 'PROTEIN'):
    """Give a list of sequcences whose residues are represented by integers,
    this function converts to a list of sequences with the residues represented
    by chars. The sequences in integer representation are obtained from
    get_alignment_int_form(file_name, biomolecule) function.

    Parameters
    ----------
        file_name : str
            FASTA file name
        biomolecule : str
            Type of biomolecule (protein or RNA)

    Returns
    -------
        seqs_char_form : list
            A list of sequences whose residues are represented by characters.
    """
    biomolecule = biomolecule.strip().upper()
    seqs_int_form = get_alignment_int_form(
        file_name,
        biomolecule=biomolecule,
    )

    msg = '\n\tConverting sequences back to character representation'
    logger.info(msg)

    RES_TO_CHAR = res_to_char(biomolecule)

    seqs_char_form = []
    for seq in seqs_int_form:
        seq_char = ''.join([RES_TO_CHAR[res] for res in seq])
        seqs_char_form.append(seq_char)
    return seqs_char_form


def sequences_to_char_form(seqs_lst, biomolecule):
    """Give a list of sequences whose residues are represented by integers, this
    function converts the sequences in which the residues are represented by
    characters. Note that this function does not read from FASTA files.

    Parameters
    ----------
        seqs_lst : list
            List of sequences each in int represenation

    Returns
    -------
        seqs_char_lst : list
            List of sequences each in char representation
    """
    biomolecule = biomolecule.strip().upper()
    RES_TO_CHAR = res_to_char(biomolecule)

    seqs_char_form = []
    for seq_int in seqs_lst:
        seq_char = [RES_TO_CHAR[res] for res in seq_int]
        seqs_char_form.append(''.join(seq_char))
    return seqs_char_form


if __name__ == '__main__':
    """
    from argparse import ArgumentParser
    import pdb
    pdb.set_trace()
    logging.basicConfig()
    parser = ArgumentParser(description='Parse args when loading fasta_reader as __main__')
    parser.add_argument('msa_file', help='alignment file in FASTA format')
    parser.add_argument('biomolecule', help='name of biomolecule', choices=['protein', 'PROTEIN', 'rna', 'RNA'])
    args = parser.parse_args()
    sequences_int_form = get_alignment_int_form(args.msa_file, args.biomolecule)
    sequences_char_form = sequences_to_char_form(sequences_int_form, args.biomolecule)

    for seq_int, seq_char in zip(sequences_int_form, sequences_char_form):
        print(seq_int, len(seq_int))
        print(seq_char)
    """
