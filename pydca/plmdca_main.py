from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.fasta_reader.fasta_reader import get_alignment_from_fasta_file
from pydca.dca_utilities import dca_utilities
import pydca.plmdca as _plmdca
from argparse import ArgumentParser
import ctypes 
import logging
import sys
import os
import glob 


"""Implements Python wrapper for the plmDCA implemented in c++. In addition,
it defines command line interface to carry out DCA computation using the 
pseudolikelihood maximization algorithm.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)


class PlmDCAException(Exception):
    """Implements exceptions related to PlmDCA computation
    """

class PlmDCA:
    """Wraps  the c++ implementation of plmDCA.

    Attributes
    ----------
        plmdca_so : str 
            Path to the plmDCA shared object created from the C++ source code
    """
    import pydca.plmdca.plmdca_location as _plmdca_loc
    plmdca_so  = glob.glob(os.path.abspath(
            os.path.join(os.path.dirname(_plmdca_loc.__file__), '_plmdca*'))
    )
    try:
        plmdca_lib_path = os.path.abspath(plmdca_so[0])
    except IndexError:
        logger.error('\n\tUnable to find plmdca dynamic library path.' 
                '\nAre you running  pydca as a module before installation?'
                '\nIn this case you need to build the plmdca shared object.'
        )
        raise 


    def __init__(self, msa_len, biomolecule, seqid=None, lambda_h=None, 
            lambda_J=None, num_iter_steps=None):
        """Initilizes plmdca instances
        """

        self.__biomolecule = biomolecule.strip().upper()
        self.__plmdca = ctypes.CDLL(self.plmdca_lib_path)
        self.__msa_len = msa_len 
        self.__num_dca_pairs = int(msa_len * (msa_len  - 1) / 2)
        self.__seqid = 0.8 if seqid is None else seqid 
        if self.__seqid <= 0 or self.__seqid > 1.0: 
            logger.error('\n\t{} is an invalid value of sequences identity (seqid) parameter'.format(self.__seqid))
            raise PlmDCAException 
        self.__lambda_h=0.01 if lambda_h is None else lambda_h
        if self.__lambda_h < 0 :
            logger.error('\n\tlambda_h must be a positive number. You passed lambda_h={}'.format(self.__lambda_h))
            raise PlmDCAException  
        self.__lambda_J=0.01 if lambda_J is None else lambda_J
        if self.__lambda_J < 0: 
            logger.error('\n\tlambda_J must be a positive number. You passed lambda_J={}'.format(self.__lambda_J))
            raise PlmDCAException


        
        log_message="""Created plmDCA instance with:
            biomolecule: {}
            MSA sequence length: {}
            sequence identity: {}
            lambda_h: {}
            lambda_J: {}
        """.format(self.__biomolecule, self.__msa_len, self.__seqid, self.__lambda_h, self.__lambda_J)
        logger.info(log_message)
        return None  

    
    def test_plmdca(self):
        """Tests plmDCA computation
        """
        self.__test_plmdca = self.__plmdca.test_plmdca 
        self.__test_plmdca.argtypes = (ctypes.c_uint, )
        self.__test_plmdca.restype = ctypes.POINTER(ctypes.c_double * self.__num_dca_pairs)
        # Get scores
        dca_scores_ptr = self.__test_plmdca(self.__msa_len)
        dca_scores  = [ score for score in dca_scores_ptr.contents ]
        counter = 0
        for i in range(self.__msa_len -1):
            for j in range(i + 1, self.__msa_len):
                print('pair:({}, {}) : {}'.format(i, j, dca_scores[counter]))
                counter += 1
        # Free memory  
        self.__free_dca_scores_cont = self.__plmdca.free_dca_scores_cont 
        self.__free_dca_scores_cont.argtypes = (ctypes.POINTER(ctypes.c_double * self.__num_dca_pairs), )
        self.__free_dca_scores_cont.restypes = None 
        self.__free_dca_scores_cont(dca_scores_ptr)
        return None
# End of PlmDCA class definition 


def configure_logging():
    """Configures logging. When configured, the logging level is INFO and
    messages are logged to stream handler. Log level name are colored whenever
    the terminal supports that. INFO level is Green, WARNING level is Yellow and
    ERROR level is Red.
    """
    from pydca.config_dca.config_log import LOGGING_CONFIG
    from pydca.config_dca.config_log import ConsoleColor as c_color
    import logging.config

    logging.config.dictConfig(LOGGING_CONFIG)
    logging.addLevelName(logging.INFO, '{}{}{}'.format(
        c_color.green, logging.getLevelName(logging.INFO), c_color.nocolor))
    logging.addLevelName(logging.WARNING, '{}{}{}'.format(
        c_color.yellow, logging.getLevelName(logging.WARNING), c_color.nocolor))
    logging.addLevelName(logging.ERROR, '{}{}{}'.format(
        c_color.red, logging.getLevelName(logging.ERROR), c_color.nocolor))
    return None


class CmdArgs:
    """Defines command line argument variables for plmDCA.
    """
    subcommand_name = 'subcommand_name'
    msa_file = 'msa_file'
    msa_file_help = 'Multiple sequence alignment (MSA) file in FASTA format'
    biomolecule = 'biomolecule'
    biomolecule_help = """Type of biomolecule.
    It should be either protein or RNA in lower or upper case letters.
    """
    verbose_optional = '--verbose'
    verbose_optional_help = 'Show logging information on the terminal.'
    subcommand_name_help = 'Subcommands destination'
    seqid_optional = '--seqid'
    seqid_optional_help = """Cut-off value of sequences similarity above which they
    are lumped together.
    """
    lambda_h_optional = '--lambda_h'
    lambda_h_optional_help = 'Value of fields penalizing constant'
    lambda_J_optional = '--lambda_J'
    lambda_J_optional_help = 'Value of couplings penalizing constant'
    num_iter_steps_optional = '--num_iter_steps'
    num_iter_steps_optional_help = 'Number of pseudolikelihood maximization iteration steps'
# end of class CmdArgs 


def get_plmdca_inst(msa_file, biomolecule, seqid=None, lambda_h=None, lambda_J=None, 
        num_iter_steps=None):
    """Creates a PlmDCA instance and returns it.

    Parameters
    ----------
        msa_file : str 
            Path to FASTA formatted MSA file.
        biomolecule : str
            Type of biomolecule the MSA data represents.
        seqid : float
            Sequences identity cut-off value.
        lambda_h : float 
            Value of fileds penalizing constant. 
        lambda_J : float
            Value of couplings penalizing constant. 
        num_iter_steps : int 
            Number of iteration for gradient decent.

    Returns
    ------- 
        plmdca_inst : PlmDCA 
            An instance of PlmDCA class
    """
    # get alignment length from MSA
    alignment = get_alignment_from_fasta_file(msa_file)
    msa_len = len(alignment[0])
    plmdca_inst = PlmDCA(msa_len, biomolecule, 
        seqid=seqid, lambda_h=lambda_h, 
        lambda_J=lambda_J, num_iter_steps=num_iter_steps,
    )
    return plmdca_inst 


def execute_from_command_line(msa_file, biomolecule, seqid=None, lambda_h=None, lambda_J=None,
        num_iter_steps=None, verbose=False):
    """Runs plmdca computation from the command line.

    Parameters
    ----------
        msa_file : str 
            Path to FASTA formatted MSA file.
        biomolecule : str
            Type of biomolecule the MSA data represents.
        seqid : float
            Sequences identity cut-off value.
        lambda_h : float 
            Value of fileds penalizing constant. 
        lambda_J : float
            Value of couplings penalizing constant. 
        num_iter_steps : int 
            Number of iteration for gradient decent. 
        verbose : bool 
            True or False. Determines if plmdca computation is done in verbose mode or not. 
    """
    if verbose : configure_logging()
    plmdca_inst = get_plmdca_inst(msa_file, biomolecule, seqid=seqid, 
        lambda_h=lambda_h, lambda_J=lambda_J, num_iter_steps=num_iter_steps
    )
    plmdca_inst.test_plmdca()
    return None 


def run_plm_dca():
    """Performs plmDCA computation based on argument passed from the command line. 
    """
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest = CmdArgs.subcommand_name)
    parser_compute_fn = subparsers.add_parser('compute_fn', help='computes the Frobenius norm of couplings.')
    parser_compute_fn.add_argument(CmdArgs.biomolecule, help=CmdArgs.biomolecule_help)
    parser_compute_fn.add_argument(CmdArgs.msa_file, help=CmdArgs.msa_file_help)
    parser_compute_fn.add_argument(CmdArgs.seqid_optional, help=CmdArgs.seqid_optional_help, type=float)
    parser_compute_fn.add_argument(CmdArgs.lambda_h_optional, help=CmdArgs.lambda_h_optional_help, type=float)
    parser_compute_fn.add_argument(CmdArgs.lambda_J_optional, help=CmdArgs.lambda_J_optional_help, type=float)
    parser_compute_fn.add_argument(CmdArgs.num_iter_steps_optional, help=CmdArgs.num_iter_steps_optional_help, type=int)
    parser_compute_fn.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_optional_help, action='store_true')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    args_dict = vars(args)

    execute_from_command_line(args_dict.get('msa_file'), args_dict.get('biomolecule'),
        seqid=args_dict.get('seqid'),
        lambda_h=args_dict.get('lambda_h'),
        lambda_J = args_dict.get('lambda_J'),
        num_iter_steps = args_dict.get('num_iter_steps'),
        verbose = args_dict.get('verbose'),
    )
    return None 

if __name__ == "__main__":
    """
    """
    run_plm_dca()