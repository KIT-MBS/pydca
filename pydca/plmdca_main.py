from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.dca_utilities import dca_utilities
from argparse import ArgumentParser
import ctypes 
import logging
import sys
import os
import glob 


"""Implements Python wrapper for the plmDCA implemented in c++. In addition,
it defines command line interface to carry out DCA computation using the 
pseudolikelihoom maximization algorithm.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)

class PlmDCA:
    """
    """
    plmdca_lib_path = os.path.abspath(glob.glob('pydca/plmdca/*_plmdca*')[0])
    def __init__(self):
        """
        """
        self.__plmdca = ctypes.CDLL(self.plmdca_lib_path)
        return None  


    
    def test_plmdca(self):
        """
        """
        self.__test_plmdca = self.__plmdca.test_plmdca 
        self.__test_plmdca.argtypes = None 
        self.__test_plmdca.restype = None 
        self.__test_plmdca()
        return None  



    def __del__(self):
        """Frees memory that has been used to return data from c++ extern functions
        to python
        """
        pass 


if __name__ == "__main__":
    test_plmdca = PlmDCA()
    print('shared library for plmdca: ', test_plmdca.plmdca_lib_path)
    test_plmdca.test_plmdca()
