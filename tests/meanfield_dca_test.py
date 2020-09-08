from __future__ import absolute_import
from pydca.meanfield_dca.meanfield_dca import MeanFieldDCA
from .input_files_path import InputFilesPath
import unittest

from Bio import AlignIO

class MeanFieldDCATestCase(unittest.TestCase):
    """Test MeanFieldDCA instance behaviour
    """
    def setUp(self):
        #rna test files
        self.__rna_msa_file = InputFilesPath.rna_msa_file
        self.__rna_ref_file = InputFilesPath.rna_ref_file
        #protein test files
        self.__protein_msa_file = InputFilesPath.protein_msa_file
        self.__protein_ref_file = InputFilesPath.protein_ref_file


        self.__mfdca_instance_protein = MeanFieldDCA(
            self.__protein_msa_file,
            'protein',
        )
        self.__mfdca_instance_rna = MeanFieldDCA(
            self.__rna_msa_file,
            'rna',
        )


    def test_compute_sorted_DI_rna(self):
        """
        """
        #self.__mfdca_instance_protein.compute_sorted_DI()
        sorted_DI = self.__mfdca_instance_rna.compute_sorted_DI()

    def test_compute_sorted_DI_protein(self):
        """
        """
        sorted_DI = self.__mfdca_instance_protein.compute_sorted_DI()


class MeanFieldDCAInputTestCase(unittest.TestCase):
    """
    Consistency test for decoupling msa handling from msa file format
    """
    def setUp(self):
        self.__protein_msa_file = InputFilesPath.protein_msa_file
        self.__protein_ref_file = InputFilesPath.protein_ref_file

    def test_input(self):
        mfdca_file = MeanFieldDCA(self.__protein_msa_file, 'protein')
        fnapc_file = mfdca_file.compute_sorted_FN_APC()

        # read MSA
        msa = AlignIO.read(self.__protein_msa_file, 'fasta')


        mfdca = MeanFieldDCA(msa, 'protein')
        fnapc = mfdca.compute_sorted_FN_APC()

        self.assertEqual(fnapc, fnapc_file)


if __name__ == '__main__':
    unittest.main()
