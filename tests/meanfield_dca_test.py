from __future__ import absolute_import
from pydca.meanfield_dca.meanfield_dca import MeanFieldDCA
from .input_files_path import InputFilesPath
import unittest

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

if __name__ == '__main__':
    unittes.main()
