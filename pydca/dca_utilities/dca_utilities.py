from pydca.fasta_reader import fasta_reader
import numpy as np
import os, errno
import logging
import zipfile

logger = logging.getLogger(__name__)

def create_directories(the_path):
    """Creat directory (directories) given a path.
    Parameters
    ----------
        the_path : str
            A string for the name of directory (directories if nested)
    Returns
    -------
        None
    """
    try:
        os.makedirs(the_path)
    except OSError as e:
        if e.errno !=errno.EEXIST:
            logger.error('Unable to create directory using path {}'.format(
            the_path))
            raise
    return None


def get_dca_output_file_path(output_dir, msa_file_name, prefix='', postfix=''):
    """Locates the file path to which DCA ouput can be written.

    Parameters
    -----------
        output_dir : str
            DCA computation related output directory.
        msa_file_name : str
            Name of the alignment file. This is used to create output file by
            adding relevant prefix and/or postfix to it.
        prefix : str
            A string that will be leading the base  name of msa file during
            output file creation.
        postfix : str
            A string that will be trailing the base name of msa file during
            output file creation.

    Returns
    -------
        output_file_path : str
            relative path to the output file name

    """
    msa_file_base_name = os.path.basename(msa_file_name)
    msa_file_root, ext = os.path.splitext(msa_file_base_name)
    output_file_base_name = prefix.strip() + msa_file_root.strip() + postfix.strip()
    output_file_path = os.path.join(output_dir, output_file_base_name)
    return output_file_path


def make_archive(root_dir, dest_dir = None):
    """Given path of directory root_dir  with with base name foo, i.e., root_dir
    has absolute path of the form /path/to/my/foo, all the files in foo are
    archived. By default the archived file will be /path/to/my/foo.zip or if
    the destination directory dest_dir is provided, it will be
    /path/to/destination/foo.zip.

    Parameters
    ----------
        root_dir : str
            The starting directory of files to be archived.
        dest_dir : str
            Directory to save the archive

    Returns
    -------
        archive_name : str
            The absolute path of archived file
    """

    root_dir = os.path.abspath(root_dir)
    if not os.path.isdir(root_dir):
        logger.error('\n\tUnable to find the following directory'
            ' to make archive:\n\t{}'.format(root_dir))
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_dir)
    root_dir_basename = os.path.basename(root_dir)
    dest_dir =  os.path.dirname(root_dir) if not dest_dir else os.path.abspath(dest_dir)
    archive_name = os.path.join(dest_dir, root_dir_basename + '.zip')
    logger.info('\n\tArchiving files from directory:{}'
        '\n\tBase name of root directory:{}'
        '\n\tArchive destination:{}'
        '\n\tArchive file name:{}'.format(root_dir, root_dir_basename,
            dest_dir, archive_name)
    )
    with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip_h:
    #start archiving
        for dir_name, sub_dirs, the_files in os.walk(root_dir):
            #dir_name has absolute path
            ##relative_path is the relative path between parent of root_dir and
            ##dir_name
            relative_path = os.path.relpath(dir_name, os.path.dirname(root_dir))
            for file_ in the_files:
                file_abs_path = os.path.abspath(os.path.join(dir_name, file_))
                name_in_archive = os.path.join(relative_path, file_)
                logger.info('\n\tFile added to archive:{}'.format(name_in_archive))
                zip_h.write(file_abs_path,  name_in_archive)

    return archive_name


def mfdca_param_metadata(mfdca_instance):
    """Generates a list of mean field DCA metadata that can be added as
    an output file header.

    Parameters
    ----------
        mfdca_instance : MeanfieldDCA
            An instance of MeanFieldDCA class.

    Returns
    -------
        metadata : list
            A list of metada related to DCA parameters used.
    """
    metadata = [
        '# PARAMETERS USED FOR THIS COMPUTATION: ',
        '#      Sequence type: {}'.format(mfdca_instance.biomolecule),
        '#      Total number of sequences in alignment data: {}'.format(
            mfdca_instance.num_sequences),
        '#      Length of sequences in alignment data: {}'.format(
            mfdca_instance.sequences_len),
        '#      Effective number of sequences: {}'.format(
            mfdca_instance.effective_num_sequences),
        '#      Value of sequence identity: {}'.format(
            mfdca_instance.sequence_identity),
        '#      Value of relative pseudocount: {}'.format(
            mfdca_instance.pseudocount),
    ]
    return metadata


def mfdca_residue_repr_metadata(biomolecule):
    """Generates single-character: interger mapping of residues

    Parameters
    ----------
        biomolecule : str
            The biomolecule type (protein or RNA, case insensitive).

    Returns
    -------
        metadata_list : list
            A list containing a string of the residue mapping. The string is put
            in a list so that it can be concatenated with other meta data which
            are list of data.

    """
    metadata_list = ['# RESIDUES IDENTIFICATION']
    mapping_dict = fasta_reader.res_to_char(biomolecule)
    char_int_pair_list = mapping_dict.items()
    int_char_pair_list_sorted = sorted(char_int_pair_list, key=lambda k:k[0])

    #put five (int, char) res pairs per line
    num_rows = int(len(int_char_pair_list_sorted) / 5)
    for i in range(num_rows + 1):
        #since slicing a list out of bound is silently truncated, we do not
        #need to check the upper bound, we get the behaviour we wanted here.
        row = int_char_pair_list_sorted[i * 5:(i + 1) * 5]
        row.insert(0, '# ')
        metadata_list.append(''.join(map(str, row)))
    return metadata_list


def get_ranked_pairs(sorted_DI, site_mapping = None):
    """Puts ranked site pairs in a list. In particular, when a site mapping is
    provided, the site pairs are put in a list according to that mapping.

    Parameters
    ----------
        sorted_DI : list
            A list of tuples, the tuples containing the site pairs and score.
            It is sorted according to DCA score (in descending order).

    Returns
    -------
        ranked_pairs : list
            A list containing tuples of site pairs.
    """
    ranked_pairs_list = []
    if site_mapping is not None:
        for pair, score in sorted_DI:
            try:
                mapped_pair = tuple(
                    [site_mapping[pair[0]], site_mapping[pair[1]]]
                )
            except KeyError:
                pass
            else:
                ranked_pairs_list.append(mapped_pair)
    else:# if no mapping is done, get all pairs from sorted DI.
        for pair, score in sorted_DI:
            ranked_pairs_list.append(pair)

    return ranked_pairs_list

def write_sorted_dca_scores(file_name, sorted_DI, site_mapping=None,
    metadata=None, score_type = None):
    """Writes sorted direct information to file.

    Parameters
    ----------
        file_name : str
            Path to file where DCA data is going to be saved.
        sorted_DI : list
            A list containing site pair
        site_mapping : dict
            A dictionary of mappings between the sites in the reference sequence
            and that of the MSA.
        metadata : list
            A list of strings that describes the metadata to be written to ouput
            file.
        score_type : str
            DCA score type. Average prodcut corrected or raw DCA score.

    Returns
    -------
        None
    """
    logger.info('\n\tWriting DCA scores to file {}'.format(file_name))
    mapped_DI = list()
    if site_mapping: # Map sites when reference sequence is supplied
        mapping_key_file = os.path.join(os.path.dirname(file_name), 'mapping_key.txt')
        logger.info('\n\tWriting site mapping key in file: {}'.format(mapping_key_file))
        mapping_key_pairs = list(site_mapping.items())
        with open(mapping_key_file, 'w') as mpk_fh:
            mpk_fh.write('# The first column represents sites in the MSA.'
                '\n# The second is the corresponding position in '
                '\n# the reference sequence. For example 9--->1 means the Ninth site in'
                '\n# the MSA is mapped to the First in the reference.\n'
            )
            for pair in mapping_key_pairs:
                mpk_fh.write('{}--->{}\n'.format(pair[0] + 1, pair[1] + 1))
        num_skipped_pairs = 0
        for pair, score in sorted_DI: # iterate over all unique pairs (i, j) for all j > i
                try:                  # and corresponding DCA scores.
                    mapped_pair = tuple([site_mapping[pair[0]], site_mapping[pair[1]]])
                except KeyError:
                    num_skipped_pairs += 1
                else:
                    if mapped_pair[0] > mapped_pair[1]:
                        logger.error('\n\tIncorrect ordering of mapped residues.'
                        '\n\tfirst residue position ({}) is greater than'
                        '\n\tsecond residue position ({})'.format(mapped_pair[0],
                                mapped_pair[1]
                            )
                        )
                    pair_and_score = tuple([mapped_pair, score])
                    mapped_DI.append(pair_and_score)
        logger.info('\n\tTotal number of site pairs in MSA not mapped to'
            ' reference: {}'.format(num_skipped_pairs))
    if site_mapping and not mapped_DI:
        logger.warn('\n\tSearch for mapping sites based on site mapping'
            ' dictionary failed.\n Writing to file skipped')
        return None

    if mapped_DI: sorted_DI = mapped_DI

    with open(file_name, 'w') as fh:
        fh.write('#' + '='*70 + '\n')
        if metadata:
            for line in metadata: fh.write('{}\n'.format(line))
        fh.write('# The First and Second columns represent sites and the'
            '\n# Third column is {} DCA score\n'.format(score_type))
        fh.write('#' + '='*70 + '\n')
        for pair, score in sorted_DI:
                i, j = pair
                fh.write('{0:<7} {1:<14} {2:<35}\n'.format(i+1, j+1, score))
    return None


def get_couplings_for_pair(couplings = None, pair = None, num_site_states = None):
    """Extracts the couplings corresponding to a give site pair.

    Parameters
    ----------
        couplings : np.array(dtype=np.float64)
            A 2d numpy array of the couplings
        pair : tuple
            A tuple of two sites forming a pair.

    Returns
    -------
        couplings_ij : np.array(dtype=np.float64)
            A 2d numpy array of the couplings corresponding sites i,j forming
            a pair.
    """
    q = num_site_states - 1 # here q is number of residues (gaps execluded)
    row_start, column_start = pair[0]*q, pair[1]*q
    row_end, column_end = row_start + q, column_start + q
    couplings_ij = couplings[row_start:row_end, column_start:column_end]

    return couplings_ij


def write_couplings(file_name, couplings, num_site_states = None,
        ranked_site_pairs = None, metadata = None, seqs_len=None):
    """Writes the couplings to file.

    Parameters
    -----------
        file_name : str
            Path to the output file to which the couplings are to be saved.
        couplings : np.array(dtype=np.float64)
            A 2d nupy array containing the couplings.
        ranked_site_pairs : list
            A list of tuples, the tuples representing ranked site pairs by DCA
            score. This argument is used to write the couplings in sorted oreder
            of the DCA score of the pairs, but the couplings themselves are not
            sorted within a site pair.
        metadata : str
            Data for file header.

    Returns
    -------
            None.
    """

    logger.info('\n\tSaving couplings to file:\n\t{}'.format(file_name))
    q = num_site_states
    with open(file_name, 'w') as fh:
        #write header meta data
        fh.write('#'+ '='*70 + '\n')
        if metadata:
            for data in metadata:
                fh.write('{}\n'.format(data))
            fh.write('# In the data below, First and Second integers refer to \n'
                '# sites and the Third and Fourth integers refer to residue.\n'
                '# The last values are the couplings\n'
                '# Residues are mapped as shown above\n')
            fh.write('#' + '='*70 + '\n')
        #write the couplings
        if ranked_site_pairs is not None:
            for pair in ranked_site_pairs:
                couplings_ij = get_couplings_for_pair(couplings = couplings,
                    pair = pair, num_site_states = num_site_states,
                )
                for a in range(num_site_states - 1):
                    for b in range(num_site_states - 1):
                        couplings_ij_ab = couplings_ij[a,b]
                        fh.write('{},{},{},{},{}\n'.format(
                            pair[0] + 1, pair[1] + 1, a + 1, b + 1, couplings_ij_ab)
                        )
        else: # write couplings in order of sites (no ranking by DCA score)
            for i in range(seqs_len - 1):
                for j in range(i + 1, seqs_len):
                    pair = (i, j)
                    couplings_ij = get_couplings_for_pair(couplings=couplings,
                        pair = pair, num_site_states= num_site_states,
                    )
                    for a in range(num_site_states - 1):
                        for b in range(num_site_states - 1):
                            couplings_ij_ab = couplings_ij[a,b]
                            fh.write('{},{},{},{},{}\n'.format(
                                pair[0] + 1, pair[1] + 1, a + 1, b + 1, couplings_ij_ab)
                            )
    return None


def write_fields(file_name, fields, num_site_states=None, metadata=None):
    """Writes local fields to file.

    Parameters
    ----------
        file_name : str
            Path to output file name to write local fields.
        fields : dict
            A dictionary whose keys are sites and values are a 1d numpy array of
            fields corresponding to states of that particular site.
        num_site_states : int
            Total number of states at sites. It takes either five or 21 for RNAs
            and proteins, respectively.
        site_mapping: dict
            A dictionary that provides mapping between the sites of the MSA data
            and that of a reference sequence.
        metadata : list
            A list containing metadata to be written to the header of fields
            output file.
    """
    logger.info('\n\tWriting fields to file {}'.format(file_name))
    num_sites = len(fields.keys())
    with open(file_name, 'w') as fh:
        fh.write('#{}\n'.format(70*'='))
        if metadata is not None:
            for data in metadata:
                fh.write('{}\n'.format(data))
            fh.write('# In the data below, the first integers referes to sites,\n'
                '# the Second integer to residue type and the last data\n'
                '# is local field. Residues are mapped as shown above\n'
                '#{}\n'.format(70*'=')
            )
            for i in range(num_sites):
                for a in range(num_site_states-1):
                    fh.write('{},{},{}\n'.format(i + 1, a + 1, fields[i][a]))
    return None


def write_single_site_freqs(file_name, fi, seqs_len = None,
        num_site_states = None, metadata = None):
    """Writes single site frequencies to file.

    Parameters
    ----------
        file_name : str
            Path to file to which single site frequencies are to be saved.
        fi : np.array(dtype=np.float64)
            A 2d numpy array of single site frequencies.
        metadata : str
            Data for file header.

    Returns
    -------
        None.
    """

    logger.info('\n\tSaving single site frequencies to file:\n\t{}'.format(
        file_name))
    with open(file_name, 'w') as fh:
        fh.write('#' + '=' * 70 + '\n')
        if metadata:
            for data in metadata:
                fh.write('{}\n'.format(data))
            fh.write('# Below, the First integer refers to the site, the \n'
                '# Second the residue at that site, and the Third is the \n'
                '# frequency. Residue numbers are mapped as shown above.\n')
            fh.write('#' + '='*70 + '\n')
        for i in range(seqs_len):
            for a in range(num_site_states):
                fh.write('{},{},{}\n'.format(i + 1, a + 1, fi[i,a]))

    return None


def write_pair_site_freqs(file_name, fij, seqs_len = None,
    num_site_states = None, metadata = None):
    """Writes pair site frequenceis to file.

    Parameters
    ----------
        file_name : str
            Name of file to which pair site frequencies are to be saved.
        fij : np.array(dtype=np.float64)
            A 2d numpy array of pair site frequencies.
        metadata : str
            Data for file header.

    Returns
    -------
        None
    """

    logger.info('\n\tSaving pair site frequencies (gaps are excluded)'
        ' to file: \n\t{}'.format(file_name))
    with open(file_name, 'w') as fh:
        fh.write('#' + '=' * 70 + '\n')
        if metadata:
            for data in metadata:
                fh.write('{}\n'.format(data))
            fh.write('# Below, the First and Second integers refer to sites, the \n'
                '# Third and Fourth residues, and the Last one is frequency for pairs.\n'
                '# Residue numbers are mapped as shown above.\n')
            fh.write('#' + '=' * 70 + '\n')
        pair_counter = 0
        for i in range(seqs_len -1):
            for j in range(i + 1, seqs_len):
                for a in range(num_site_states - 1):
                    for b in range(num_site_states -1):
                        fh.write('{},{},{},{},{}\n'.format(
                            i + 1, j + 1, a + 1, b + 1, fij[pair_counter, a, b]))
                pair_counter += 1

    return None


def write_params(fields_file_path=None, couplings_file_path=None, fields=None,
        couplings=None, num_site_states=None, metadata=None):
    """Writes the couplings and fields to their corresponding output files.

    Parameters
    ----------
        fields_file_path : str
            File to write fields (of the global probability model)
        couplings_file_path : str
            File to write couplings
        fields : np.array
            Numpy 2d array containing the fields.
        coupling : np.array
            Numpy 2d array of the couplings.
        num_site_states : int
            The total number of site states. 5 and 21 for RNAs and proteins,
            respectively.
        metadata : list
            A list of metadata to be written to output files.

    Returns
    -------
        None : None
    """
    # write fields
    write_fields(
        fields_file_path, fields,
        num_site_states=num_site_states,
        metadata = metadata,
    )
    #write couplings
    seqs_len = len(fields.keys())
    write_couplings(
        couplings_file_path, couplings,
        num_site_states=num_site_states,
        metadata = metadata, seqs_len=seqs_len
    )
    return None


def get_dcavisualizer_metadata(dcavisualizer_inst):
    """Tabulates a list of information about the parameters used
    in a DCAVisualizer instance.

    Parameters
    ----------
        dcavisualizer_inst : DCAVisualizer
            An instacne of a DCAVisualizer class.

    Returns
    -------
        dcavisualizer_metadata : list(strs)
            A list of strings, containing information about some of the
            dcavisualizer_inst parameter.
    """

    contact_dist = dcavisualizer_inst.contact_dist
    linear_dist = dcavisualizer_inst.linear_dist
    pdb_id = dcavisualizer_inst.pdb_id
    pdb_chain_id = dcavisualizer_inst.pdb_chain_id
    wc_neighbor_dist = dcavisualizer_inst.wc_neighbor_dist
    biomolecule = dcavisualizer_inst.biomolecule
    dcavisualizer_metadata = [
        '# PARAMETES USED FOR THIS COMPUTATION',
        '#\tMinimum PDB contact distance : {}'.format(contact_dist),
        '#\tLinear distance between residues in chain > : {}'.format(linear_dist),
        '#\tWC neighbor distance (if RNA) : {}'.format(wc_neighbor_dist),
        '#\tBIOMOLECULE : {}'.format(biomolecule),
        '#\tPDB-ID : {}'.format(pdb_id),
        '#\tPDB-CHAIN-ID : {}'.format(pdb_chain_id),
        '# First and Second columns are the positions of contacting residues in',
        '# referece sequence. The Third column is an annotation of contact',
        '# category. The categories can be:',
        '# tp->true posiitve, fp->false positives, pdb->PDB contacts,',
        '# missing->missing in PDB chain, tp-wc->true positive and WC pair (RNA)',
        '# tp-nwc->true positive and non-WC (RNA)'
    ]
    return dcavisualizer_metadata


def write_tp_rate(file_name, true_positive_rates_dict=None, metadata=None):
    """Writes the true positive rate to output file

    Parameters
    ----------
        file_name : str
            path to output file where the the TP rate is going to be written.
        tp_rates : list(tuples)
            A list of tuples containing the rank and corresponding true positive
            rate.
        metadata : list(strs)
            A list of strings for containing information about the parameters
            used for true positive rate computation. This metadata will be writen
            at the header of the output file (fila_name).

    Returns
    -------
        None : None
    """
    dca_tp_rates = true_positive_rates_dict['dca']
    pdb_tp_rates = true_positive_rates_dict['pdb']
    with open(file_name, 'w') as fh:
        fh.write('#' + '='*70 + '\n')
        for data in metadata:
            fh.write('{}\n'.format(data))
        fh.write('#' + '='*70 + '\n')
        for dca_tpr, pdb_tpr in zip(dca_tp_rates, pdb_tp_rates):
            fh.write('{0:.6f}\t{1:.6f}\n'.format(dca_tpr, pdb_tpr))

    return None


def write_contact_map(file_name, contact_categories_dict, metadata=None):
    """Writes the contact map to output file.

    Parameters
    ----------
        file_name : str
            Path to output file
        contact_categories_dict : dict
            A dictionary of categorized contacts whose keys are contact type,
            i.e., true positives, false positives, pdb, or missing. The values
            are lists of tuples of site pairs.
    Returns
    -------
        None : None
    """
    describe_columns = [
        '# Column-1 :  contact category',
        '# Column-2 : site-number in sequence (first pairing site)',
        '# Column-3 : site-number in sequence (second pairing site)',
        '# Column-4 : closest atom pairs for residue pairs',
        '# Column-5 : site-number in PDB (first pairing site)',
        '# Column-6 : site-number in PDB (second pairing site)',
        '# Column-7 : distance between pairing atoms (column-4) in Angstrom',

    ]
    metadata.extend(describe_columns)
    with open(file_name, 'w') as fh:
        fh.write('#' + '='*70 + '\n')
        for data in metadata:
            fh.write('{}\n'.format(data))
        fh.write('#' + '='*70 + '\n')
        for category in contact_categories_dict:
            for pair in contact_categories_dict[category]:
                p = list(pair)
                pdb_metadata = list(contact_categories_dict[category][pair])
                line = p + pdb_metadata
                line.insert(0, category)
                fmted_line = '\t\t'.join(str(elem) for elem in line)
                fh.write(fmted_line + '\n')

    return None


def write_trimmed_msa(file_name, msa_trimmer=None, columns_to_remove=None, metadata=None):
    """Writes trimmed MSA file in FASTA format

    Parameters
    ----------
        file_name : str
            Path to file where trimmed MSA data is going to be written
        msa_trimmer : MSATrimmer
            An instance of MSATrimmer class
        columns_to_remove : list/tuple
            A list of columns to be trimmed from the MSA
        metadata : list/tuple
            A list/tuple of metadata about trimming parameters used to trim the
            MSA

    Returns
    -------
        None : None
    """
    logger.info('\n\tWritting trimmed MSA in to file {}'.format(file_name))
    with open(file_name, 'w') as fh:
        for record in msa_trimmer.alignment_data:
            trimmed_seq = [
                record.seq[i] for i in range(len(record.seq)) if i not in columns_to_remove
            ]
            fh.write('>{}\n{}\n'.format(record.id, ''.join(trimmed_seq)))
    return None


if __name__  == '__main__':
    """
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('root_dir', help='root di')
    parser.add_argument('--dest_dir', help='destination directory for the archive')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    make_archive(args.root_dir, dest_dir = args.dest_dir)
    """
