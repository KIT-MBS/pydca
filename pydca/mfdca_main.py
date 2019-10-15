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
    force_seq_type_optional = '--force_seq_type'
    force_seq_type_help = """Typically the anticipated number of residue types
    plus a gap is 21 for protein, and 5 for RNA sequences. If there is a significant
    deviation from these values, it is interpreted as the user has inadvertently
    entered a biomolecule type mismatch and an error may be raised.
    This can happen when the alignment data contains too many/few non-standared
    residues or when a wrong biomolecule type is entered. If you are sure about
    the biomolecule type the MSA data represents, use --force_seq_type to bypass
    this error.
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
    pdb_file = 'pdb_file'
    pdb_file_help = """Path to a PDB file.
    """
    dca_file = 'dca_file'
    dca_file_help = """File containing the result of DCA computation. The first
    and second columns should contain site pairs (i, j) such that i < j. Optionally,
    a third column can contain the DCA scores. The DCA scores are not mandatory
    as the site-pairs are assumed to be sorted, in descending order of DCA score.
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
        pseudocount=pseudocount,seqid=seqid,
        force_seq_type=force_seq_type)
    return mfdca_inst


def add_args_to_subparser(the_parser, subcommand_name):
    """Adds arguments to a pasrser. These added arguments are common to all
    sub-commands that are used for mean-field DCA computation.

    Parameters
    ----------
        the_parser : ArgumentParser
            The (sub) parser that to which the arguments are added.

    Returns
    -------
        None
    """

    the_parser.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_help,
        action='store_true',
    )
    if subcommand_name in DCA_COMPUTATION_SUBCOMMANDS:
        the_parser.add_argument(CmdArgs.biomolecule, help = CmdArgs.biomolecule_help)
        the_parser.add_argument(CmdArgs.msa_file, help = CmdArgs.msa_file_help)
        the_parser.add_argument(CmdArgs.apc_optional, help=CmdArgs.apc_help,
            action = 'store_true',
        )
        the_parser.add_argument(CmdArgs.seqid_optional, help = CmdArgs.seqid_help,
            type = float,
        )
        the_parser.add_argument(CmdArgs.pseudocount_optional,
            help=CmdArgs.pseudocount_help,
            type = float,
        )
        the_parser.add_argument(CmdArgs.refseq_file_optional,
            help = CmdArgs.refseq_file_help
        )
        the_parser.add_argument(CmdArgs.output_dir_optional,
            help=CmdArgs.output_dir_help,
        )
        the_parser.add_argument(CmdArgs.force_seq_type_optional,
            help = CmdArgs.force_seq_type_help,
            action = 'store_true',
        )


def execute_from_command_line(msa_file=None, biomolecule=None, seqid=None,
        pseudocount=None, the_command=None, refseq_file = None,
        force_seq_type=False, verbose=False, output_dir=None, apc=False):
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
        force_seq_type : bool
            Force computation although there might be a biomolecule type  mismatch.
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

    ######  start dca computation #############################################
    if the_command.strip() in DCA_COMPUTATION_SUBCOMMANDS:
        mfdca_instance = get_mfdca_instance(msa_file, biomolecule, seqid=seqid,
            pseudocount=pseudocount, force_seq_type=force_seq_type,
        )

        param_metadata = dca_utilities.mfdca_param_metadata(mfdca_instance)
        #create path to output directory is not supplied by user
        if not output_dir:
            msa_file_base_name, ext = os.path.splitext(os.path.basename(msa_file))
            output_dir = 'MFDCA_output_' + msa_file_base_name
        #create dca coutput directory
        dca_utilities.create_directories(output_dir)
        #### Compute DCA score
        if the_command.strip()=='compute_di' or the_command.strip()=='compute_couplings':
            mapped_sites = None
            if refseq_file:# do backmapping when reference sequence file is provided
                seq_backmapper = SequenceBackmapper(
                    alignment_data=mfdca_instance.alignment,
                    refseq_file = refseq_file,
                    biomolecule = mfdca_instance.biomolecule)
                mapped_sites = seq_backmapper.map_to_reference_sequence()
            if apc: # do average product correction if apc is passed from the command line
                sorted_DI, couplings = mfdca_instance.compute_sorted_DI_APC()
                score_type = ' MF DI average product corrected (APC)'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='MFDCA_apc_di_scores_', postfix='.txt'
                )
            else: # compute raw DCA score if apc is not asked
                sorted_DI, couplings = mfdca_instance.compute_sorted_DI()
                score_type = 'raw DI'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='MFDCA_raw_di_scores_', postfix='.txt'
                )
            if the_command.strip() == 'compute_di':# write DI scores
                dca_utilities.write_sorted_dca_scores(di_file_path,sorted_DI,
                    site_mapping=mapped_sites, metadata=param_metadata,
                    score_type = score_type
                )
            else: #  if the command is compute couplings
                ranked_pairs_list =  dca_utilities.get_ranked_pairs(
                    sorted_DI, site_mapping = mapped_sites
                )
                couplings_file_path = dca_utilities.get_dca_output_file_path(
                    output_dir, msa_file, prefix='couplings_', postfix='.txt'
                )
                residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                    mfdca_instance.biomolecule
                )
                #save couplings to file
                metadata  = param_metadata + residue_repr_metadata
                dca_utilities.write_couplings(
                    couplings_file_path, couplings,
                    num_site_states = mfdca_instance.num_site_states,
                    ranked_site_pairs = ranked_pairs_list,
                    metadata = metadata,
                )
        # compute Frobenius norm of couplings
        if the_command.strip()=='compute_fn':
            mapped_sites = None
            if refseq_file:# do backmapping when reference sequence file is provided
                seq_backmapper = SequenceBackmapper(
                    alignment_data=mfdca_instance.alignment,
                    refseq_file = refseq_file,
                    biomolecule = mfdca_instance.biomolecule)
                mapped_sites = seq_backmapper.map_to_reference_sequence()
            if apc:
                score_type = 'MFDCA Frobenius norm, average product corrected (APC)'
                sorted_FN = mfdca_instance.compute_sorted_FN_APC()
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'MFDCA_apc_fn_scores_', postfix='.txt'
                )
            else:
                score_type = 'MFDCA raw Frobenius norm'
                sorted_FN = mfdca_instance.compute_sorted_FN()
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'MFDCA_raw_fn_scores_', postfix='.txt'
                )
            dca_utilities.write_sorted_dca_scores(fn_file_path, sorted_FN,
                site_mapping = mapped_sites, metadata = param_metadata,
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
            dca_utilities.write_fields(fields_file_path, fields, metadata=metadata,
                num_site_states = mfdca_instance.num_site_states,
            )
        # compute fields and couplings
        if the_command.strip() == 'compute_params':
            fields, couplings = mfdca_instance.compute_hamiltonian()
            residue_repr_metadata = dca_utilities.mfdca_residue_repr_metadata(
                mfdca_instance.biomolecule
            )
            metadata = param_metadata + residue_repr_metadata
            fields_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                msa_file, prefix='fields_', postfix='.txt',
            )
            couplings_file_path = dca_utilities.get_dca_output_file_path(
                output_dir, msa_file, prefix='couplings_', postfix='.txt',
            )
            dca_utilities.write_params(
                fields_file_path=fields_file_path,
                couplings_file_path=couplings_file_path,
                fields=fields,
                couplings=couplings,
                metadata=metadata,
                num_site_states = mfdca_instance.num_site_states,
            )

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
            # pass --pseudocount 0.0 to compute raw fijs
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

    #direct information computation parser
    parser_compute_di = subparsers.add_parser('compute_di',
        help = 'Computes the direct information.'
        ' Example: mfdca compute_di <biomolecule> <MSA> --verbose, where'
        ' <biomolecule> takes a value protein or rna (case insensitive)'
        ' <MSA> takes path to the MSA file',
    )
    add_args_to_subparser(parser_compute_di, 'compute_di')
    parser_compute_fn = subparsers.add_parser('compute_fn',
        help = 'Compute the Frobenius norm of couplings.'
            ' Example: see compute_di',
    )
    add_args_to_subparser(parser_compute_fn, 'compute_fn')

    #couplings computation parser
    parser_compute_couplings = subparsers.add_parser('compute_couplings',
        help = 'Computes the couplings for a protein or RNA from MSA input data.'
            ' The biomolecule type (protein/RNA) must be specified. In addition,'
            ' a file containing a reference sequence can be supplied. For more '
            ' information about this command use --help.'
    )
    add_args_to_subparser(parser_compute_couplings, 'compute_couplings')

    # fields computation parser
    parser_compute_fields = subparsers.add_parser('compute_fields',
        help = 'Computes the global probability fields'
    )
    add_args_to_subparser(parser_compute_fields, 'compute_fields')

    # parameters (fields and couplings) computation parser
    parser_compute_params = subparsers.add_parser('compute_params',
        help = 'Computes the parameters of global probability model, i.e., '
            ' couplings and fields in one run.'
    )
    add_args_to_subparser(parser_compute_params, 'compute_params')

    #Single site frequencies computation parser
    parser_compute_fi = subparsers.add_parser('compute_fi',
        help = 'Computes regularized single-site frequencies from MSA.'
            ' If raw frequencies are desired, use --pseudocount 0'
    )
    add_args_to_subparser(parser_compute_fi, 'compute_fi')

    #pair site frequencies computation parser
    parser_compute_fij = subparsers.add_parser('compute_fij',
        help = 'Computes regularized pair-site frequencies from MSA. If raw'
        ' frequenceis are desired, set the pseudocount to zero. Use --help'
        ' for more information.'
    )
    add_args_to_subparser(parser_compute_fij, 'compute_fij')

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
        force_seq_type=args_dict.get('force_seq_type'),
        verbose = args_dict.get('verbose'),
        output_dir = args_dict.get('output_dir'),
        apc = args_dict.get('apc'),
    )
    logger.info('\n\tDONE')
    return None


if __name__ == '__main__':
    #run DCA computation when this file is loaded as __main__
    run_meanfield_dca()
