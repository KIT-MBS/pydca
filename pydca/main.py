from pydca.meanfield_dca import meanfield_dca
from pydca.sequence_backmapper.sequence_backmapper import SequenceBackmapper
from pydca.dca_utilities import dca_utilities
from pydca.contact_visualizer.contact_visualizer import DCAVisualizer, PDBContent
from pydca.msa_trimmer.msa_trimmer import MSATrimmer
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
    rna_secstruct_file = 'rna_secstruct_file'
    rna_secstruct_file_optional = '--rna_secstruct_file'
    rna_secstruct_file_help = """File containing an RNA secondary structure. The
    structure should be in one line and in dot bracket notation, i.e., the allowed
    symbols are '.', '(', and ')'. Comments should start with # sign.
    """
    wc_neighbor_dist_optional = '--wc_neighbor_dist'
    wc_neighbor_dist_help = """For RNAs, if two residues (i, j) form
    secondary structure pair and the wc_neighbor_dist is n,  all (2n + 1)(2n + 1)
    (neighboring plus WC pairs) are excluded. That is all pairs (-n + i, -n + j),
    (-n + 1 + i, -n + j), ..., (i, j), ..., (i + n, j + n). Using
    wc_neighbor_dist = 2 excludes 25 pairs.
    """
    pdb_chain_id = 'pdb_chain_id'
    pdb_chain_id_help = """ID of a PDB chain. This helps to identify the desired
    chain since PDB files can contain multiple protein/RNA chains or a
    mix of protein and RNA chains in one file.
    """
    pdb_id_optional = '--pdb_id'
    pdb_id_help = """The PDB ID in the PDB database.
    """

    linear_dist = 'linear_dist'
    linear_dist_optional = '--linear_dist'
    linear_dist_help = """The distance between two residues in sequence above
    which they are considered to be tertiary contacts. For RNAs see also
    wc_neighbor_dist parameter.
    """
    contact_dist = 'contact_dist'
    contact_dist_optional = '--contact_dist'
    contact_dist_help = """Cut-off distance between two residues in a PDB structure
    below which they are considered to be in contact. This distance is between two
    heavy atoms of the residues.
    """
    num_dca_contacts = 'num_dca_contacts'
    num_dca_contacts_optional = '--num_dca_contacts'
    num_dca_contacts_help = """Number of DCA contacts to be taken among all DCA
    ranked contacts. The counting is done after appropriate filtering of the DCA
    ranked contacts is done if a filtering criteria (e.g., a non-zero linear distance)
    is set.
    """
    max_gap_optional = '--max_gap'
    max_gap_help = """The maximum fraction of gaps in MSA column. When an MSA
    data is trimmed, columns containing gap fraction more than this value are
    removed.
    """
    remove_all_gaps_optional = '--remove_all_gaps'
    remove_all_gaps_help = """Removes columns in the MSA correponding to
    the matching sequence's gap positions.
    """
## end of class CmdArgs

DCA_COMPUTATION_SUBCOMMANDS = ['compute_di', 'compute_fn','compute_couplings',
    'compute_fields', 'compute_params','compute_fi', 'compute_fij',
]

DCA_VISUALIZATION_SUBCOMMANDS = ['plot_contact_map', 'plot_tp_rate']

FILE_CONTENT_SUBCOMMANDS = ['pdb_content', 'refseq_content', 'dca_content',
    'rna_secstruct_content',
]

MSA_TRIMMING_SUBCOMMANDS = ['trim_by_refseq', 'trim_by_gap_size']
# all subcommands
ALL_SUBCOMMANDS = list()
ALL_SUBCOMMANDS.extend(DCA_COMPUTATION_SUBCOMMANDS)
ALL_SUBCOMMANDS.extend(DCA_VISUALIZATION_SUBCOMMANDS)
ALL_SUBCOMMANDS.extend(FILE_CONTENT_SUBCOMMANDS)
ALL_SUBCOMMANDS.extend(MSA_TRIMMING_SUBCOMMANDS)

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
        pseudocount=pseudocount,sequence_identity=seqid,
        force_seq_type=force_seq_type)
    return mfdca_inst


def get_dcavisualizer_instance(biomolecule, pdb_chain_id, pdb_file, refseq_file,
        dca_file = None, rna_secstruct_file=None, linear_dist=None,
        contact_dist=None, num_dca_contacts=None,
        wc_neighbor_dist=None, pdb_id=None):
    """Creates and returns a DCAVisualizer instance

    Parameters
    ----------
        biomolecule : str
            The biomolecule type (protein or RNA)
        pdb_chain_id : str
            Chain ID of a PDB chain in PDB file.
        pdb_file : str
            Path to PDB file
        refseq_file : str
            Path to FASTA formated reference sequence file.
        dca_file : str
            Path to text file containing DCA ranked pairs.
        rna_secstruct_file : str
            Path to text file containing an RNA secondary structure.
        linear_dis : int
            Distance between two residues in a sequence.
        contact_dist : float
            Distance between two residues below which they are considered to be
            contacts.
        pdb_id : str
            The PDB ID as obtained from PDB database.
    Returns
    -------
        dcavisualizer_inst : DCAVisualizer
            An instance of DCAVisualizer class.
    """
    dcavisualizer_inst = DCAVisualizer(biomolecule, pdb_chain_id, pdb_file,
        refseq_file=refseq_file, dca_file=dca_file,
        rna_secstruct_file=rna_secstruct_file, linear_dist=linear_dist,
        contact_dist=contact_dist, num_dca_contacts=num_dca_contacts,
        wc_neighbor_dist=wc_neighbor_dist, pdb_id=pdb_id,
    )
    return dcavisualizer_inst


def get_pdb_content_instance(pdb_file):
    """Creates a PDBContent instance and returns it.

    Parameters
    ----------
        pdb_file : str
            Path to a PDB file.

    Returns
    -------
        pdb_content_instance : PDBContent
            An instance of PDBContent.
    """
    pdb_content = PDBContent(pdb_file)
    return pdb_content


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
    if subcommand_name in DCA_VISUALIZATION_SUBCOMMANDS:
        the_parser.add_argument(CmdArgs.biomolecule, help=CmdArgs.biomolecule_help)
        the_parser.add_argument(CmdArgs.pdb_chain_id, help=CmdArgs.pdb_chain_id_help)
        the_parser.add_argument(CmdArgs.pdb_file, help=CmdArgs.pdb_file_help)
        the_parser.add_argument(CmdArgs.refseq_file, help=CmdArgs.refseq_file_help)
        the_parser.add_argument(CmdArgs.dca_file, help=CmdArgs.dca_file_help)
        the_parser.add_argument(CmdArgs.rna_secstruct_file_optional,
                help=CmdArgs.rna_secstruct_file_help,
            )
        the_parser.add_argument(CmdArgs.linear_dist_optional,
            help=CmdArgs.linear_dist_help, type = int,
        )
        the_parser.add_argument(CmdArgs.contact_dist_optional,
            help=CmdArgs.contact_dist_help, type = float,
        )
        the_parser.add_argument(CmdArgs.num_dca_contacts_optional,
            help = CmdArgs.num_dca_contacts_help, type = int,
        )
        the_parser.add_argument(CmdArgs.wc_neighbor_dist_optional, type= int,
            help = CmdArgs.wc_neighbor_dist_help,
        )
        the_parser.add_argument(CmdArgs.pdb_id_optional, help = CmdArgs.pdb_id_help)

    if subcommand_name in FILE_CONTENT_SUBCOMMANDS:
        if subcommand_name == 'pdb_content':
            the_parser.add_argument(CmdArgs.pdb_file, help = CmdArgs.pdb_file_help)
    if subcommand_name in MSA_TRIMMING_SUBCOMMANDS:
        the_parser.add_argument(CmdArgs.max_gap_optional,
            type = float, help = CmdArgs.max_gap_help,
        )
        if subcommand_name == 'trim_by_refseq':
            the_parser.add_argument(CmdArgs.biomolecule, help=CmdArgs.biomolecule_help)
            the_parser.add_argument(CmdArgs.msa_file, help=CmdArgs.msa_file_help)
            the_parser.add_argument(CmdArgs.refseq_file, help=CmdArgs.refseq_file_help)
            the_parser.add_argument(CmdArgs.remove_all_gaps_optional,
                help= CmdArgs.remove_all_gaps_help, action='store_true',
            )
        if subcommand_name == 'trim_by_gap_size':
            the_parser.add_argument(CmdArgs.msa_file, help=CmdArgs.msa_file_help)
    return None


def execute_from_command_line(msa_file=None, biomolecule=None, seqid=None,
        pseudocount=None, the_command=None, refseq_file = None,
        force_seq_type=False, verbose=False, output_dir=None, apc=False,
        pdb_file=None, pdb_chain_id=None, dca_file=None, rna_secstruct_file=None,
        linear_dist=None, contact_dist=None, num_dca_contacts=None,
        wc_neighbor_dist=None, pdb_id=None, max_gap=None, remove_all_gaps=False):
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
        pdb_file : str
            PDB file path.
        pdb_chain_id : str
            The name of a PDB chain in PDB file.
        dca_file : str
            DCA file path.
        rna_secstruct_file : str
            RNA secondary structure file path.
        linear_dist : int
            Minimum separation between two residues in sequence.
        contact_dist : float
            Maximum distance between two residues in PDB to be considered as contacts.
        num_dca_contacts : int
            Number of DCA contacts to be taken from ranked pairs.
        wc_neighbor_dist : int
            Maximum radius for a residue to be considered as a neighbor of WC pairs.
        pdb_id : str
            The ID of a PDB as it is in PDB repository.
        max_gap : float
            Maximum fraction of gaps in an MSA column.
        remove_all_gaps : float
            Remove all gaps in MSA correponding to the best matching sequence to
            a reference.


    Returns
    -------
        None : None 
    """


    if verbose: configure_logging()


    if the_command not in ALL_SUBCOMMANDS:
        logger.error('\n\t{} is unknown command.'.format(the_command))
        raise ValueError
    ######  start dca computation #############################################
    if the_command.strip() in DCA_COMPUTATION_SUBCOMMANDS:
        mfdca_instance = get_mfdca_instance(msa_file, biomolecule, seqid=seqid,
            pseudocount=pseudocount, force_seq_type=force_seq_type,
        )

        param_metadata = dca_utilities.mfdca_param_metadata(mfdca_instance)
        #create path to output directory is not supplied by user
        if not output_dir:
            msa_file_base_name, ext = os.path.splitext(os.path.basename(msa_file))
            output_dir = 'DCA_output_' + msa_file_base_name
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
                score_type = ' DI average product corrected (APC)'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='DCA_apc_di_scores_', postfix='.txt'
                )
            else: # compute raw DCA score if apc is not asked
                sorted_DI, couplings = mfdca_instance.compute_sorted_DI()
                score_type = 'raw DI'
                di_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix='DCA_raw_di_scores_', postfix='.txt'
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
                score_type = 'Frobenius norm, average product corrected (APC)'
                sorted_FN = mfdca_instance.compute_sorted_FN_APC()
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'DCA_apc_fn_scores_', postfix='.txt'
                )
            else:
                score_type = 'raw Frobenius norm'
                sorted_FN = mfdca_instance.compute_sorted_FN()
                fn_file_path = dca_utilities.get_dca_output_file_path(output_dir,
                    msa_file, prefix = 'DCA_raw_fn_scores_', postfix='.txt'
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
    ###---- end of if the_command.strip() in DCA_COMPUTAION_SUBCOMMANDS -------
    if the_command.strip() in DCA_VISUALIZATION_SUBCOMMANDS:
        dca_visualizing_commands = ['plot_contact_map', 'plot_tp_rate']
        if the_command.strip() in dca_visualizing_commands:
            # get an instance of DCAVisualizer
            dcavisualizer = get_dcavisualizer_instance(
                biomolecule, pdb_chain_id, pdb_file, refseq_file=refseq_file,
                dca_file = dca_file, rna_secstruct_file=rna_secstruct_file,
                linear_dist = linear_dist, contact_dist=contact_dist,
                num_dca_contacts = num_dca_contacts,
                wc_neighbor_dist = wc_neighbor_dist,
                pdb_id = pdb_id,
            )
            # get metadata about dcavisualizer
            dcavisualizer_metadata = dca_utilities.get_dcavisualizer_metadata(
                dcavisualizer
            )
        if the_command.strip() == 'plot_contact_map':
            contact_categories_dict = dcavisualizer.plot_contact_map()
            if not output_dir:
                pdb_file_basename, ext = os.path.splitext(os.path.basename(pdb_file))
                output_dir = 'contact_map_' + pdb_file_basename
            output_file_path = dca_utilities.get_dca_output_file_path(
                output_dir, pdb_file, prefix= 'contact_map', postfix='.txt'
            )
            dca_utilities.create_directories(output_dir)
            dca_utilities.write_contact_map(output_file_path,
                contact_categories_dict, metadata = dcavisualizer_metadata
            )

        if the_command.strip() == 'plot_tp_rate':
            true_positive_rates_dict = dcavisualizer.plot_true_positive_rates()
            if not output_dir:
                pdb_file_basename, ext = os.path.splitext(os.path.basename(pdb_file))
                output_dir = 'TPR_' + pdb_file_basename
            output_file_path = dca_utilities.get_dca_output_file_path(
                output_dir, pdb_file, prefix= 'TPR_', postfix='.txt'
            )
            tpr_metadata = ['\n# First column is DCA true positive rate per rank'
                '\n# Second column is the PDB true positive rate per rank'
            ]
            metadata = dcavisualizer_metadata[:6] + tpr_metadata
            dca_utilities.create_directories(output_dir)
            dca_utilities.write_tp_rate(output_file_path,
                true_positive_rates_dict=true_positive_rates_dict,
                metadata=metadata,
            )

    ### if the_command.strip() in FILE_CONTENT_SUBCOMANDS
    if the_command.strip() in FILE_CONTENT_SUBCOMMANDS:
        if the_command.strip() == 'pdb_content':
            pdb_content = get_pdb_content_instance(pdb_file)
            pdb_content.show_struct_info()

    ### MSA trimming commands
    if the_command.strip() in MSA_TRIMMING_SUBCOMMANDS:
        if the_command.strip() == 'trim_by_refseq':
            msa_trimmer = MSATrimmer(msa_file,
                biomolecule=biomolecule,
                refseq_file=refseq_file,
                max_gap= max_gap,
            )
            columns_to_remove = msa_trimmer.trim_by_refseq(
                remove_all_gaps=remove_all_gaps,
            )
        if the_command.strip() == 'trim_by_gap_size':
            msa_trimmer = MSATrimmer(msa_file, max_gap=max_gap)
            columns_to_remove = msa_trimmer.trim_by_gap_size()
        if not output_dir:
            msa_file_basename, ext = os.path.splitext(os.path.basename(msa_file))
            output_dir = 'Trimmed_' +  msa_file_basename
            dca_utilities.create_directories(output_dir)
        output_file_path = dca_utilities.get_dca_output_file_path(
            output_dir, msa_file, prefix= 'Trimmed_', postfix='.fa'
        )
        dca_utilities.write_trimmed_msa(output_file_path, msa_trimmer=msa_trimmer,
            columns_to_remove=columns_to_remove,
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

    parser_dca_visualizer_contact_map = subparsers.add_parser('plot_contact_map',
        help = 'Provides a quick contact map comparison of DCA computation results.'
            ' For a given PDB chain ID, the PDB contacts are extracted from PDB'
            ' file. Then top N ranked DCA pairs are used for the contact map'
            ' comparison. Use --help for more information.'
    )
    add_args_to_subparser(parser_dca_visualizer_contact_map, 'plot_contact_map')
    parser_dca_visualizer_tp_rate = subparsers.add_parser('plot_tp_rate',
        help = ' Plots the true positive rate per rank of a DCA computation'
            ' result. The DCA file should contain ranked pairs (i, j) such that'
            ' i < j. If the biomolecule is RNA, a secondary structure file should'
            ' be provided if one is interested on tertiary contacts. Use --help'
            ' for more information.'

    )
    add_args_to_subparser(parser_dca_visualizer_tp_rate, 'plot_tp_rate')

    parser_pdb_content = subparsers.add_parser('pdb_content',
        help = 'Displays information about the contents of a PDB file.'
            ' Use --verbose optional argument to display the PDB summary'
            ' on the terminal.',
    )
    add_args_to_subparser(parser_pdb_content, 'pdb_content')

    parser_trim_by_refseq = subparsers.add_parser('trim_by_refseq',
        help='Removes MSA columns containing fraction of gaps more than'
        ' the value specified by {} (default 0.5) if these columns'
        ' do not correspond to residues of the sequence in MSA that matches'
        ' with the reference. Setting max_gap to zero removes all columns'
        ' except those corresponding to the residues of the matching sequence'
        ' to the reference.'.format(CmdArgs.max_gap_optional)
    )
    add_args_to_subparser(parser_trim_by_refseq, 'trim_by_refseq')

    parser_trim_by_gap_size = subparsers.add_parser('trim_by_gap_size',
        help = 'Removes MSA columns containing gap fraction more than the value'
            ' specified by {} (default 0.5)'.format(CmdArgs.max_gap_optional)
    )
    add_args_to_subparser(parser_trim_by_gap_size, 'trim_by_gap_size')

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
        pdb_file = args_dict.get('pdb_file'),
        pdb_chain_id = args_dict.get('pdb_chain_id'),
        dca_file = args_dict.get('dca_file'),
        rna_secstruct_file = args_dict.get('rna_secstruct_file'),
        linear_dist = args_dict.get('linear_dist'),
        contact_dist = args_dict.get('contact_dist'),
        num_dca_contacts = args_dict.get('num_dca_contacts'),
        wc_neighbor_dist = args_dict.get('wc_neighbor_dist'),
        pdb_id = args_dict.get('pdb_id'),
        max_gap = args_dict.get('max_gap'),
        remove_all_gaps = args_dict.get('remove_all_gaps'),
    )
    logger.info('\n\tDONE')
    return None


if __name__ == '__main__':
    #run DCA computation when this file is loaded as __main__
    run_meanfield_dca()
