from __future__ import absolute_import, division
import unittest
import os
import glob
from pydca.fasta_reader import fasta_reader
from .input_files_path import InputFilesPath
class TestCase(unittest.TestCase):
    def setUp(self):
        """
        """
        self.__rna_msa_file = InputFilesPath.rna_msa_file
        self.__protein_msa_file = InputFilesPath.protein_msa_file
        self.__rna = 'rna'
        self.__protein = 'protein'



    def test_get_alignment_from_fasta_file(self):
        rna_seqs = fasta_reader.get_alignment_from_fasta_file(
            self.__rna_msa_file,
        )
        self.assertIsNotNone(rna_seqs)
        protein_seqs = fasta_reader.get_alignment_from_fasta_file(
            self.__protein_msa_file,
        )
        self.assertIsNotNone(protein_seqs)

    def test_get_alignment_int_form(self):
        rna_seqs_int_form = fasta_reader.get_alignment_int_form(
            self.__rna_msa_file,
            biomolecule = self.__rna,
        )
        self.assertIsNotNone(rna_seqs_int_form)
        protein_seqs_int_form = fasta_reader.get_alignment_int_form(
            self.__protein_msa_file,
            biomolecule = self.__protein,
        )

        self.assertIsNotNone(protein_seqs_int_form)


if __name__ == '__main__':
    unittest.main()
