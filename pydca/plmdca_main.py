from pydca.plmdca import plmdca 
from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.fasta_reader.fasta_reader import get_alignment_from_fasta_file
from pydca.dca_utilities import dca_utilities
from argparse import ArgumentParser
import logging
import sys


"""Top level module for plmDCA. Defines command line interface and 
configures logging.

Author: Mehari B. Zerihun
"""

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
    num_iter_steps_optional = '--num_iterations'
    num_iter_steps_optional_help = 'Number of pseudolikelihood maximization iteration steps'
# end of class CmdArgs 


def get_plmdca_inst(biomolecule, msa_file, seqid=None, lambda_h=None, lambda_J=None, 
        num_iterations=None):
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
    plmdca_inst = plmdca.PlmDCA(biomolecule, msa_file, 
        seqid=seqid, lambda_h=lambda_h, 
        lambda_J=lambda_J, num_iterations = num_iterations,
    )
    return plmdca_inst 


def execute_from_command_line(biomolecule, msa_file, seqid=None, lambda_h=None, lambda_J=None,
        num_iterations=None, verbose=False):
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
    plmdca_inst = get_plmdca_inst(biomolecule, msa_file, seqid=seqid, 
        lambda_h=lambda_h, lambda_J=lambda_J, num_iterations = num_iterations
    )
    plmdca_inst.compute_params()
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

    execute_from_command_line(args_dict.get('biomolecule'),args_dict.get('msa_file'), 
        seqid=args_dict.get('seqid'),
        lambda_h=args_dict.get('lambda_h'),
        lambda_J = args_dict.get('lambda_J'),
        num_iterations = args_dict.get('num_iterations'),
        verbose = args_dict.get('verbose'),
    )
    return None 

if __name__ == "__main__":
    """
    """
    run_plm_dca()