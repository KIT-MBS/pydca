from __future__ import  absolute_import
import unittest
from pydca.sequence_backmapper import sequence_backmapper as seq_backmapper
from pydca.fasta_reader import fasta_reader
from .input_files_path import  InputFilesPath


class SequenceBackmapperTestCase(unittest.TestCase):

    def setUp(self):
        """
        """
        self.__rna_msa_file = InputFilesPath.rna_msa_file
        self.__rna_ref_file = InputFilesPath.rna_ref_file
        self.__protein_msa_file = InputFilesPath.protein_msa_file
        self.__protein_ref_file = InputFilesPath.protein_ref_file


        rna_alignment_int_form = fasta_reader.get_alignment_int_form(
            self.__rna_msa_file, biomolecule = 'rna')

        self.__rna_backmapper = seq_backmapper.SequenceBackmapper(
            alignment_data = rna_alignment_int_form,
            refseq_file = self.__rna_ref_file,
            biomolecule = 'rna',
        )

        protein_alignment_int_form = fasta_reader.get_alignment_int_form(
            self.__protein_msa_file, biomolecule='protein')

        self.__protein_backmapper = seq_backmapper.SequenceBackmapper(
            alignment_data = protein_alignment_int_form,
            refseq_file = self.__protein_ref_file,
            biomolecule = 'protein'
        )


    def test_map_to_reference_sequence(self):
        mapped_sites_rna = self.__rna_backmapper.map_to_reference_sequence()
        self.assertTrue(len(mapped_sites_rna) > 1)
        mapped_sites_protein = self.__protein_backmapper.map_to_reference_sequence()
        self.assertTrue(len(mapped_sites_protein) > 1)


if __name__ == '__main__':
    unittest.main()
