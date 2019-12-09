from pydca.meanfield_dca import meanfield_dca
from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.dca_utilities import dca_utilities
from argparse import ArgumentParser
import logging
import sys
import os

"""Top level module for the pydca package. It implements command line
intefaces, including logging configuration, command line arguments and help
messages.

Author: Mehari B. Zerihun
"""

logger = logging.getLogger(__name__)

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
    """Defines variables related to command line parsing.
    """
    subcommand_name = 'subcommand_name'
    #variables for command line use
    output_dir_optional = '--output_dir'
    output_dir_help = """Directory path to which output results are written.
    If the directory is not existing, it will be created provided that the user
    has a privilege to do so. If this path is not provided, an output directory
    is created using the base name of the MSA file, with a prefix and/or postfix
    added to it.
    """
    verbose_optional = '--verbose'
    verbose_help = 'Show logging information on the terminal.'
    apc_optional = '--apc'
    apc_help = """Compute the average product corrected (APC) DCA score.
    """

    msa_file = 'msa_file'
    msa_file_help = 'Multiple sequence alignment (MSA) file in FASTA format'

    biomolecule = 'biomolecule'
    biomolecule_help = """Type of biomolecule.
    It should be either protein or RNA in lower or upper case letters.
    """
    refseq_file = 'refseq_file'
    refseq_file_optional = '--refseq_file'
    refseq_file_help = """FASTA formatted file containing a reference sequence.
    The reference sequence should not contain gaps or non-standard residues.
    """

    pseudocount_optional = '--pseudocount'
    pseudocount_help = """Relative value of the pseudocount so as to regularize
    emperical single-site and pair-site frequencies obtained from alignment data.
    Note that this is not the actual value of the pseudocount, instead, it is
    the ratio X/(X + Meff) where X is the actual pseudocount and Meff is the
    effective number of sequences.
    """

    seqid_optional = '--seqid'
    seqid_help = """Cut-off value of sequences similarity above which they
    are lumped together.
    """

    ranked_by_optional = '--ranked_by'
    ranked_by_optional_help="""Method in which DCA scores are calculated. There are
    four options: direct information (DI), Frobenius norm (FN) and their average
    product corrected forms (DI_APC, FN_APC).
    """

    linear_dist_optional = '--linear_dist'
    linear_dist_help="""Minimum separation beteween site pairs in sequence. 
    """
    num_site_pairs_optional = '--num_site_pairs'
    num_site_pairs_help = """The maximum number of site pairs whose couplings are 
    to be extracted.
    """
## end of class CmdArgs

DCA_COMPUTATION_SUBCOMMANDS = ['compute_di', 'compute_fn','compute_couplings',
    'compute_fields', 'compute_params','compute_fi', 'compute_fij',
]


def get_mfdca_instance(msa_file, biomolecule, force_seq_type=False, **kwargs):
    """Inititalizes the MeanFieldDCA instance based on the arguments passed
    on the command line

    Parameters
    ----------
        msa_file : str
            Path to FASTA file containing MSA data.
        biomolecule : str
            Type of biomolecule (protein or RNA)
        **kwargs : dict
            This list of keyword arguments contains:
            i)      seqid : The sequence identity.
            ii)     pseudocount : The relative pseudo count.

    Returns
    -------
        mfdca_inst : MeanFieldDCA
            A MeanFieldDCA instance.
    """
    #mfdca_inst = meanfield_dca.MeanFieldDCA(msa_file, biomolecule)
    seqid = kwargs.get('seqid')
    pseudocount = kwargs.get('pseudocount')
    mfdca_inst = meanfield_dca.MeanFieldDCA(msa_file, biomolecule,
        pseudocount=pseudocount,seqid=seqid
    )
    return mfdca_inst


def execute_from_command_line(msa_file=None, biomolecule=None, seqid=None,
        pseudocount=None, the_command=None, refseq_file = None,
        verbose=False, output_dir=None, apc=False, ranked_by=None,
        linear_dist=None, num_site_pairs=None):
    """Do computations according to the parameters passed on the command line

    Parameters
    ----------
        msa_file : str
            Path to MSA file.
        biomolecule : str
            Name of biomolecule (protein or RNA, case insensitive).
        seqid : float
            Value of sequence identity.
        pseudocount : float
            Value of relative pseudo count.
        the_command : str
            Name of the command passed from the command line.
        refseq_file : str
            Path to reference sequence file.
        verbsoe : bool
            Display logging message to the screen.
        ouput_dir : str
            Output directory to write ouput files.
        apc : bool
            Perform average product correction.

    Returns
    -------
        None : None 
    """


    if verbose: configure_logging()

    if the_command.strip() in DCA_COMPUTATION_SUBCOMMANDS:
        mfdca_instance = get_mfdca_instance(msa_file, biomolecule, seqid=seqid,
            pseudocount=pseudocount,
        )
        seqbackmapper = None
        # update mapped_sites  when refseq is provided
        if refseq_file:# do backmapping when reference sequence file is provided
            seqbackmapper = SequenceBackmapper(
                alignment_data=mfdca_instance.alignment,
                refseq_file = refseq_file,
                biomolecule = mfdca_instance.biomolecule
            )
        param_metadata = dca_utilities.mfdca_param_metadata(mfdca_instance)
        #create path to output directory is not supplied by user
        if not output_dir:
            msa_file_base_name, ext = os.path.splitext(os.path.basename(msa_file))
            output_dir = 'MFDCA_output_' + msa_file_base_name
        # create dca coutput directory
        dca_utilities.create_directories(output_dir)
        # compute DCA score
        if the_command.strip()=='compute_di':
            if apc: # do average product correction if apc is passed from the command line
                sorted_DI = mfdca_instance.compute_sorted_DI_APC(seqbackmapper=seqbackmapper)
                score_type = ' MF DI average product corrected (APC)'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='MFDCA_apc_di_scores_', postfix='.txt'
                )
            else: # compute raw DCA score if apc is not asked
                sorted_DI = mfdca_instance.compute_sorted_DI(seqbackmapper=seqbackmapper)
                score_type = 'raw DI'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='MFDCA_raw_di_scores_', postfix='.txt'
                )
            
            dca_utilities.write_sorted_dca_scores(di_file_path,sorted_DI,
                metadata=param_metadata,
                score_type = score_type
            )
        # compute Frobenius norm of couplings
        if the_command.strip()=='compute_fn':
            if apc:
                score_type = 'MFDCA Frobenius norm, average product corrected (APC)'
                sorted_FN = mfdca_instance.compute_sorted_FN_APC(seqbackmapper=seqbackmapper)
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'MFDCA_apc_fn_scores_', postfix='.txt'
                )
            else:
                score_type = 'MFDCA raw Frobenius norm'
                sorted_FN = mfdca_instance.compute_sorted_FN(seqbackmapper=seqbackmapper)
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'MFDCA_raw_fn_scores_', postfix='.txt'
                )
            dca_utilities.write_sorted_dca_scores(fn_file_path, sorted_FN,
                metadata = param_metadata,
                score_type = score_type
            )
        # compute global probability local fields
        if the_command.strip()=='compute_fields':
            fields = mfdca_instance.compute_fields()
            residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                mfdca_instance.biomolecule
            )
            metadata = param_metadata + residue_repr_metadata
            fields_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                msa_file, prefix = 'fields_', postfix='.txt'
            )
            dca_utilities.write_fields_csv(fields_file_path, fields, 
                metadata=metadata,
            )
        # compute fields and couplings
        if the_command.strip() == 'compute_params':
            fields, couplings = mfdca_instance.compute_params(
                seqbackmapper = seqbackmapper,
                ranked_by = ranked_by,
                linear_dist = linear_dist,
                num_site_pairs = num_site_pairs,
            )
            residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                mfdca_instance.biomolecule
            )
            metadata = param_metadata + residue_repr_metadata
            # write fields to text file
            fields_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'fields_', postfix='.txt'
                )
            param_metadata.append('#\tTotal number of sites whose fields are extracted: {}'.format(len(fields))) 
            dca_utilities.write_fields_csv(fields_file_path, fields, metadata=param_metadata)
            couplings_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'couplings_', postfix='.txt'
                )
            param_metadata.pop() 
            param_metadata.append('#\tTotal number of site pairs whose couplings are extracted: {}'.format(len(couplings)))
            if ranked_by is None: # the default is FN_APC
                ranked_by = 'FN_APC'
            param_metadata.append('#\tDCA ranking method used: {}'.format(ranked_by))
            if linear_dist is None: # default is |i - j| > 4
                linear_dist = 4
            param_metadata.append('#\tMinimum separation beteween site pairs in sequence: |i - j| > {}'.format(linear_dist))
            dca_utilities.write_couplings_csv(couplings_file_path, couplings, metadata=param_metadata)

        #Compute single site frequencies
        if the_command.strip() == 'compute_fi':
            #pass --pseudocount 0.0 if raw frequencies are desired
            fi = mfdca_instance.get_reg_single_site_freqs()
            residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                mfdca_instance.biomolecule)
            metadata = param_metadata + residue_repr_metadata
            fi_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                msa_file, prefix='fi_', postfix='.txt')
            dca_utilities.write_single_site_freqs(fi_file_path, fi,
                seqs_len = mfdca_instance.sequences_len,
                num_site_states = mfdca_instance.num_site_states,
                metadata = metadata)

        #Compute pair site frequencies
        if the_command.strip() == 'compute_fij':
            # pass --pseudocount 0.0 to compute raw fij
            file_path = dca_utilities.get_dca_output_file_path(output_dir,
                msa_file, prefix='fij_', postfix='.txt',
            )
            residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                mfdca_instance.biomolecule,
            )
            metadata = param_metadata + residue_repr_metadata
            fij = mfdca_instance.get_reg_pair_site_freqs()
            dca_utilities.write_pair_site_freqs(file_path, fij,
                seqs_len = mfdca_instance.sequences_len,
                num_site_states = mfdca_instance.num_site_states,
                metadata = metadata,
            )

    return None


def run_meanfield_dca():
    """Entry point for DCA computations.

    Parameters
    ----------
        All arguments to be used are captured from the command line.

    Returns
    -------
        None.
    """
    parser = ArgumentParser()

    #Create subparsers
    subparsers = parser.add_subparsers(dest = CmdArgs.subcommand_name)

    # Mean-field DI computation parser
    parser_compute_di = subparsers.add_parser('compute_di',
        help = 'Computes the direct information.'
        ' Example: mfdca compute_di <biomolecule> <MSA> --verbose, where'
        ' <biomolecule> takes a value protein or rna (case insensitive)'
        ' <MSA> takes path to the MSA file',
    )
    
    parser_compute_di.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
    parser_compute_di.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
    parser_compute_di.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help, type = float)
    parser_compute_di.add_argument(CmdArgs.pseudocount_optional, help=CmdArgs.pseudocount_help, type = float)
    parser_compute_di.add_argument(CmdArgs.refseq_file_optional, help = CmdArgs.refseq_file_help)
    parser_compute_di.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)
    parser_compute_di.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help, action='store_true')
    parser_compute_di.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help, action='store_true')

    # Mean-field FN computation parser
    parser_compute_fn = subparsers.add_parser('compute_fn',
        help = 'Compute the Frobenius norm of couplings.'
            ' Example: see compute_di',
    )
    
    parser_compute_fn.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
    parser_compute_fn.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
    parser_compute_fn.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help, type = float)
    parser_compute_fn.add_argument(CmdArgs.pseudocount_optional, help=CmdArgs.pseudocount_help, type = float)
    parser_compute_fn.add_argument(CmdArgs.refseq_file_optional, help = CmdArgs.refseq_file_help)
    parser_compute_fn.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)
    parser_compute_fn.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help, action='store_true')
    parser_compute_fn.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help, action='store_true')

    # parameters (fields and couplings) computation parser
    parser_compute_params = subparsers.add_parser('compute_params',
        help = 'Computes the parameters of global probability model, i.e., '
            ' couplings and fields in one run.'
    )
    parser_compute_params.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
    parser_compute_params.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
    parser_compute_params.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help, type = float)
    parser_compute_params.add_argument(CmdArgs.pseudocount_optional, help=CmdArgs.pseudocount_help, type = float)
    parser_compute_params.add_argument(CmdArgs.refseq_file_optional, help = CmdArgs.refseq_file_help)
    parser_compute_params.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)
    parser_compute_params.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help, action='store_true')
    parser_compute_params.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help, action='store_true')
    parser_compute_params.add_argument(CmdArgs.ranked_by_optional, help=CmdArgs.ranked_by_optional_help, 
        choices= ('FN', 'FN_APC', 'DI', 'DI_APC', 'fn', 'fn_apc', 'di', 'di_apc')
    )
    parser_compute_params.add_argument(CmdArgs.linear_dist_optional, help=CmdArgs.linear_dist_help, type=int)
    parser_compute_params.add_argument(CmdArgs.num_site_pairs_optional, help=CmdArgs.num_site_pairs_help, type=int)

    
    #Single site frequencies computation parser
    parser_compute_fi = subparsers.add_parser('compute_fi',
        help = 'Computes regularized single-site frequencies from MSA.'
            ' If raw frequencies are desired, use --pseudocount 0'
    )
    parser_compute_fi.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
    parser_compute_fi.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
    parser_compute_fi.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help, type = float)
    parser_compute_fi.add_argument(CmdArgs.pseudocount_optional, help=CmdArgs.pseudocount_help, type = float)
    parser_compute_fi.add_argument(CmdArgs.refseq_file_optional, help = CmdArgs.refseq_file_help)
    parser_compute_fi.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)
    parser_compute_fi.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help, action='store_true')
    parser_compute_fi.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help, action='store_true')
    
    #pair site frequencies computation parser
    parser_compute_fij = subparsers.add_parser('compute_fij',
        help = 'Computes regularized pair-site frequencies from MSA. If raw'
        ' frequenceis are desired, set the pseudocount to zero. Use --help'
        ' for more information.'
    )
    parser_compute_fij.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
    parser_compute_fij.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
    parser_compute_fij.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help, type = float)
    parser_compute_fij.add_argument(CmdArgs.pseudocount_optional, help=CmdArgs.pseudocount_help, type = float)
    parser_compute_fij.add_argument(CmdArgs.refseq_file_optional, help = CmdArgs.refseq_file_help)
    parser_compute_fij.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)
    parser_compute_fij.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help, action='store_true')
    parser_compute_fij.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help, action='store_true')

    #display help if no argument is passed
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    args_dict = vars(args)
    # Do the computations based on the arguments passed from the command line
    execute_from_command_line(
        biomolecule = args_dict.get('biomolecule'),
        msa_file = args_dict.get('msa_file'),
        seqid = args_dict.get('seqid'),
        pseudocount = args_dict.get('pseudocount'),
        refseq_file = args_dict.get('refseq_file'),
        the_command = args_dict.get('subcommand_name'),
        verbose = args_dict.get('verbose'),
        output_dir = args_dict.get('output_dir'),
        apc = args_dict.get('apc'),
        ranked_by = args_dict.get('ranked_by'),
        linear_dist = args_dict.get('linear_dist'),
        num_site_pairs = args_dict.get('num_site_pairs'),
    )
    logger.info('\n\tDONE')
    return None


if __name__ == '__main__':
    #run DCA computation when this file is loaded as __main__
    run_meanfield_dca()
