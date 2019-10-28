import os

class InputFilesPath:
    """Defines paths to input files for testing.
    attributes:
        rna_msa_file:
            Absolute path to MSA FASTA file containg RNA sequences.
        rna_ref_file:
            Absolute path to FASTA file containing RNA reference sequence
        protein_msa_file:
            Absolute path to MSA FASTA file containing protein sequences
        protein_ref_file:

    """
    rna_msa_file = os.path.abspath("tests/tests_input/MSA_RF00059_trimmed_gap_treshold_50.fa")
    rna_ref_file = os.path.abspath("tests/tests_input/ref_seq_RF00059.faa")
    protein_msa_file = os.path.abspath("tests/tests_input/PF02826.faa")
    protein_ref_file = os.path.abspath("tests/tests_input/ref_seq_PF02826.faa")
