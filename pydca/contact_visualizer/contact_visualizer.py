import matplotlib.pyplot as plt
import Bio.PDB as biopdb
import Bio.SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio import pairwise2
from argparse import ArgumentParser
import logging
from ..sequence_backmapper import scoring_matrix
import os
from collections import OrderedDict
import itertools
import requests 

logger = logging.getLogger(__name__)

STANDARD_RESIDUES = {
    'RNA' : ('A', 'C', 'G', 'U'),

    'PROTEIN':('ALA', 'ARG', 'ASN', 'ASP', 'CYS',
        'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO',
        'SER', 'THR', 'TRP', 'TYR', 'VAL')
}

RES_THREE_CHAR_TO_ONE = {
    'PROTEIN': {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
        'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V',
    },
}

STANDARD_RESIDUES['PROTEIN_ONE_CHAR'] = tuple(
    RES_THREE_CHAR_TO_ONE['PROTEIN'].values()
)


def is_protein_sequence(the_sequence):
    """Identifies wheather a sequence is protein or not. The decision is made
    based only on standared residues.

    Parameters
    ----------
        the_sequence : str
            A biological sequence.
    Returns
    -------
        True or False : bool
            Returns True if the parameter, the_sequence, contains only
            protein standared residues, False otherwise.
    """
    the_sequence = the_sequence.strip().upper()
    for res in the_sequence:
        if res not in STANDARD_RESIDUES['PROTEIN_ONE_CHAR']:
            return False
    return True


def is_rna_sequence(the_sequence):
    """Identiy wheather a sequence is RNA or not. The decision is made
    base only on standared residues.

    Parameters
    ----------
        the_sequence : str
            The sequence to be investigated.

    Returns
    -------
        True or False :  bool
            Return True if the parameter, the_sequence, contains only RNA
            standared residues, False otherwise.
    """
    the_sequence = the_sequence.strip().upper()
    for res in the_sequence:
        if res not in STANDARD_RESIDUES['RNA']:
            return False
    return True


def sequence_matches_biomolecule(the_sequence, biomolecule):
    """Verifies if a given sequence corresponds to the biomolecule type.

    Parameters
    -----------
        the_sequence : str
            The biomolecular sequence to be verified
        biomolecule : str
            The type of biomolecule (protein or RNA)
    """
    biomolecule = biomolecule.strip().upper()
    if biomolecule == 'PROTEIN':
        return is_protein_sequence(the_sequence)
    elif biomolecule == 'RNA':
        return is_rna_sequence(the_sequence)
    else :
        logger.error('\n\tUnknown biomolecule type {}'.format(biomolecule))
        raise ValueError


class PDBContentException(Exception):
    """Used for raising PDB analysis related exceptions.
    """

class PDBContent:
    """Parses PDB files and stores useful information about content of a given
    PDB file.
    """

    def __init__(self, pdb_file, biomolecule = None):
        """Initializes PDBcontent with PDB file and biomolecule type.

        Parameters
        ----------
            pdb_file :  str
                File path to the PDB file
            biomolecule :  str
                Type of the biomolecule (protein or RNA), if not given it takes
                None value
        Returns
        -------
            None
        """
        # pdb_file may be just a PDB ID. We download it if it is a valid PDB ID.
        if os.path.isfile(pdb_file):
            self.__pdb_file = pdb_file 
        else:
            pdb_file_basename = os.path.basename(pdb_file)
            if pdb_file_basename[0].isdigit() and len(pdb_file_basename)==4:
                pdb_id = pdb_file_basename.upper()
                self.__pdb_file = self.download_pdb(pdb_id)
            else:
                logger.error('\n\t{} represents neither a PDB file nor a valid PDB ID.'
                    ' Please provide a valid PDB file or PDB ID'.format(pdb_file)
                )
                raise PDBContentException
        if biomolecule is not None: biomolecule = biomolecule.strip().upper()
        if biomolecule:
            if biomolecule not in ('PROTEIN', 'RNA'):
                logger.error(
                    '\n\t{} is invalid choice. Must be protein or RNA.'.format(
                        biomolecule),
                )
                raise PDBContentException
        self.__biomolecule = biomolecule
        self.__pdb_structure =  self.get_pdb_structure()
        self.__pdb_chain_sequences = self.collect_chain_sequences()
        return None


    @property
    def pdb_chain_sequences(self):
        """sequences getter property

        Parameters
        ----------
            self : PDBContent
                An instance of PDBContent class.

        Returns
        -------
            self.__pdb_chain_sequences : OrderdDict
                A dictionary containing chain sequences. The keys are chain IDs
                and the values are tuples of sequence type and the sequence itself
                respectively.
        """
        return self.__pdb_chain_sequences


    @property
    def pdb_structure(self):
        """pdb structure getter property.

        Parameters
        ----------
            self : PDBContent
                An instance of PDBContent class.
        Returns
        -------
            self.__pdb_structure : Bio.PDB.Structure.Structure
                The PDB structure obtained from self.__pdb_file. This structure
                is the return value of Bio.PDB.PDBParser.get_structure(struct_id,
                self.__pdb_file).
        """
        return self.__pdb_structure


    def get_pdb_structure(self, pdb_name = None, show_warning = False):
        """Reads residue objects from PDB chain from the PDB structure file pdb_file.

        Parameters
        ----------
            pdb_name:
                Name of the PDB structure as given by the user.

        Returns
        -------
            pdb_structure:
                The PDB structure created by Biopython.
        """
        #if pdb_name is None:
        #    pdb_name = os.path.splitext(os.path.basename(self.__pdb_file))[0]
        PERMISSIVE = 1 if show_warning else 0
        logger.info('\n\tReading PDB structure from file: {}'.format(self.__pdb_file))
        pdb_parser = Bio.PDB.PDBParser(PERMISSIVE=PERMISSIVE) # makes no difference
        # unless in newest version of biopython (default shows warning in older versions)
        structure = pdb_parser.get_structure(pdb_name, self.__pdb_file)
        num_models = len(structure)
        logger.info('\n\tObtained PDB structure')
        logger.info('\n\tNumber of models in PDB structure {}: {}'.format(
            pdb_name, num_models)
        )
        return structure


    @staticmethod
    def download_pdb(pdb_id):
        """Downloads a PDB file

        Parameters
        ----------
            pdb_id : str
                The PDB identifier

        Returns
        -------
            pdb_file_local_path : str 
                Path to downloaded PDB file local path. Currently it is 
                in the working directory. #TODO add a destination directory?
            None
        """
        
        pdb_file_url = 'https://files.rcsb.org/view/{}.pdb'.format(pdb_id)
        logger.info('\n\tDownloading PDB file from: {}'.format(pdb_file_url))
        try:
            r = requests.get(pdb_file_url)
        except Exception:
            logger.error('\n\tProblem while downloading PDB file from: {}'.format(pdb_file_url)
            )
            raise
        else:
            pdb_file_local_path = '{}_downloaded.pdb'.format(pdb_id)
            with open(pdb_file_local_path, 'wb') as pdb_h:
                pdb_h.write(r.content)
        return pdb_file_local_path


    def extract_structure_info(self):
        """Extracts information from PDB header.

        Parameters
        ----------
            self : PDBAnalyzer()

        Returns
        -------
            struct_info : dict
                A dictionary of PDB header metadata.
        """
        structure = self.get_pdb_structure()
        info_keys = ['resolution', 'structure_method', 'name', 'head',
            'deposition_date', 'release_date',
            'compound', 'author', 'journal_reference',
        ]
        logger.info('\n\tExtracting header info from PDB file:'
            '\n\t{}'.format(self.__pdb_file)
        )
        keys_and_info = [(key, structure.header[key]) for key in info_keys]
        structure_info = OrderedDict(keys_and_info)
        return structure_info


    def show_struct_info(self):
        """Displays PDB structure info extracted from PDB file header metadata.

        Parameters
        ----------

        Returns
        -------
            None
        """
        header_info = self.extract_structure_info()
        fmt_info = ''
        for key, info in header_info.items():
            current_info = '\n\t{}: {}'.format(key, info)
            fmt_info += current_info
        underline = '-'*15
        logger.info('\n\tPDB header info:\n\t' + underline + fmt_info)
        return None


    @staticmethod
    def get_chain_id_list(structure):
        """Finds the list of chains in a PDB structures and returns a list of ids
        of the chains.

        Parameters
        ----------
            structure : Biopython PDB structure object.

        Returns
        -------
            chain_id_list : list
                A list of chain ids.
        """
        list_of_models = []
        first_model = structure[0]
        chains = first_model.get_list()
        num_chains = len(chains)
        logger.info('\n\tFound {} chains from PDB model'.format(num_chains))
        chain_id_list = list()
        for chain in chains:
            chain_id_list.append(chain.id)
        logger.info('\n\tList of chains IDs: {}'.format(chain_id_list))
        return chain_id_list


    def filter_residues(self, residues, biomolecule):
        """Filters the standared residues from others (e.g. hetatom residues).

        Parameters
        ----------
            residues : list
                A list of Biopython PDB structure residues objects.

        Returns
        -------
            standard_residues : list
                A lif of Biopython PDB structure standared residues (after hetro
                atom residues are filtered).
        """
        biomolecule = biomolecule.strip().upper()
        standard_residues = []
        for res in residues:
            if res.get_resname().strip() in STANDARD_RESIDUES[biomolecule]:
                if not res.id[0].strip(): standard_residues.append(res) # filter out hetro residues
        return standard_residues


    @staticmethod
    def to_sequence(residue_name_list, biomolecule):
        """Converts a list of residue ids to a sequence of one char residues ids.

        Parameters
        ----------
            residue_name_list :  list
                A list of residue names ordered according to the PDB chain this list
                of residues is obtained from.

        Returns
        -------
            sequence : str
                A sequence of one char represetation of residue IDs.
        """
        biomolecule = biomolecule.strip().upper()
        if biomolecule == 'PROTEIN':
            sequence = ''.join(
                [RES_THREE_CHAR_TO_ONE[biomolecule][res] for res in residue_name_list]
            )
        elif biomolecule == 'RNA':
            sequence = ''.join(residue_name_list)
        else:
            logger.error('\n\tUnknown value for biomolecule: {}'.format(biomolecule))
            raise PDBContentException
        return sequence


    def collect_chain_sequences(self):
        """Collects chain ID and corresponding sequences from PDB files

        Parameters
        ----------
            self :  PDBContent
                PDBContent instance
        Returns
        -------
            chain_seqs : OrderedDict[chain_id : tuple(biomolecule, sequence)]
                A dictionary with chain IDs as key and a tuple of biomolecule type
                and sequence as values.
        """
        structure = self.__pdb_structure
        chain_ids = self.get_chain_id_list(structure)
        model = structure[0]
        chain_seqs = OrderedDict()

        for chain_id in chain_ids:
            all_residues =  model[chain_id].get_list()
            biomolecule = 'PROTEIN'
            standard_residues = self.filter_residues(all_residues, biomolecule)
            if not standard_residues:
                biomolecule = 'RNA'
                standard_residues = self.filter_residues(all_residues, biomolecule)
            if not standard_residues :
                logger.error('\n\tUnable to obtain standared residues'
                    ' for chain {} of the PDB file {}'.format(chain_id, self.__pdb_file)
                )
                raise PDBContentException
            res_id_list =  list()
            for res in standard_residues:
                res_id_list.append(res.get_resname().strip())
            current_seq = self.to_sequence(res_id_list, biomolecule)
            chain_seqs[chain_id] = (biomolecule, current_seq)
        return chain_seqs


    def display_chain_sequences(self):
        """Displays chains IDs and corresponding sequence.

        Parameter
        ---------
            self :  PDBContent
                An instance of PDBContent class
        Returns
        -------
            None
        """
        chain_seqs = self.__pdb_chain_sequences
        formatted_msg = ''
        for chain_id in chain_seqs:
            formatted_msg += '\n\tBiomolecule: {}\n\tChain ID: {}'\
                    '\n\tSequence: {}'.format(chain_seqs[chain_id][0], chain_id,
                    chain_seqs[chain_id][1]
                )
        logger.info(formatted_msg)
        return None


class RefSeqContentException(Exception):
    """Exception to be raised when an error is occuring in the
    RefSeqContent class.
    """

class  RefSeqContent:
    """Stores and analyzes sequence data from reference sequence fasta files.

    Attributes
    ----------
        self.__refseq_file : str
            Path to reference sequence file.
        self.__biomoleule : str
            Type of the biomolecule, (protein or RNA).
        self.__sequences : list(str)
            A list of reference sequences as read from self.__refseq_file.
    """

    def __init__(self, refseq_file, biomolecule=None):
        """RefSeqContent initializer.

        Parameters
        ----------
            refseq_file : str
                Path to a text file containing reference sequence.
            biomolecule : str
                The type of biomolecule the reference sequence is intended to
                represent.

        Returns
        -------
            None
        """
        self.__refseq_file = refseq_file
        if biomolecule is not None:
            self.__biomolecule = biomolecule.strip().upper()
            if self.__biomolecule not in ('PROTEIN', 'RNA'):
                logger.error('\n\t {} is not valid biomolecule.'
                    ' must be protein (PROTEIN) or rna (RNA)'.format(biomolecule),
                )
                raise RefSeqContentException
        self.__sequences = self.get_ref_seqs_from_fasta_file()
        return None


    @property
    def ref_sequences(self):
        """sequences getter

        Parameters
        ----------
            self : RefSeqContent
                Instance of RefSeqContent
        """
        return self.__sequences


    def get_ref_seqs_from_fasta_file(self):
        """Read sequences from PDB fasta file and return the first sequence.

        Parameters
        ----------
            pdb_seq_fasta_file : str
                Fasta file containing sequences of PDB as obtained from PDB
                database.

        Returns
        -------
            ref_sequences : OrderedDict
                A dictionary of sequences with the chain IDs as keys and the
                a tuple of sequence type and the sequence as values.
        """
        ref_sequences = OrderedDict()
        seqs_read =  Bio.SeqIO.parse(self.__refseq_file, 'fasta')
        for k, record in enumerate(seqs_read, start=1):
            sequence = record.seq.upper()
            if sequence:
                seq_type = self.identify_seq_type(sequence)
                seq_str =  ''.join(list(sequence))
                ref_sequences[k]= (seq_type, seq_str)
        if not ref_sequences.values():
            logger.error('\n\tUnable to find sequences from file {}'.format(
                self.__refseq_file)
            )
            raise RefSeqContentException
        logger.info('\n\tTotal number of sequences found in'
            ' reference fasta file: {}'.format(len(ref_sequences))
        )
        seqs_list = ref_sequences.values()
        formatted_seqs = ['\n\t{}'.format(seq) for seq in seqs_list]

        return  ref_sequences


    def display_reference_sequences(self):
        """Displays the list of reference sequences recorded from input file.

        Parameters
        ----------
            self : RefSeqContent
                An instance of RefSeqContent
        Returns
        -------
            None
        """
        fmtd_msg = ''
        for id, seq in self.__sequences.items():
            fmtd_msg += '\n\tSequence {}, {} : {}'.format(id, seq[0], seq[1])
        logger.info(fmtd_msg)
        return None


    @staticmethod
    def identify_seq_type(the_sequence):
        """Identifies wheather a sequence type is protein or RNA. Only standared
        residues are considered valid.

        Parameters
        ----------
            the_sequence : str
                Protein or RNA sequence
        Returns
        -------
            seq_type : str
                The sequence type that parameter the_sequence represents. It
                takes only either of the two values: 'RNA' or 'PROTEIN'
        """
        the_sequence = the_sequence.strip().upper()

        if is_rna_sequence(the_sequence):
            seq_type = 'RNA'
        elif is_protein_sequence(the_sequence):
            seq_type = 'PROTEIN'
        else:
            logger.error('\n\tThe sequence {}'
                '\n\tis neither protein nor RNA. Does it contain non-standared residues?'
                '\n\tPlease note that reference sequences are assumed to contain'
                '\n\tonly standared residues.'.format(the_sequence)
            )
            raise RefSeqContentException

        return seq_type


class RNASecStructContentException(Exception):
    """Used for raising exceptions related to RNA secondary structure.
    """

class RNASecStructContent:
    """Stores information about RNA secondary structure
    """
    def __init__(self, secstruct_file):
        """Initializes an RNASecStructContent instance.
        Parameters
        ----------
            secstruct_file : str
                Path to RNA secondary structure file.
        Attributes
        ----------
            self.__secstruct_file : str
                Path to the location of text file containing RNA secondary structure.
            self.__left_brackets : tuple(str)
                A tuple containing symbols used as left brackets to represent
                the RNA secondary structure.
            self.__right_brackets : tuple(str)
                A tuple containing symbols used as right brackets to represent
                the RNA secondary structure.
            self.__nonwc_symbols : tuple(str)
                A tuple containing symbols used to represent non-WC nucleotides
                in the RNA secondary structure.
            self.__wcpairs : tuple(tuples)
                A tuple of WC pairs. The pairs are represented by their index
                position in the secondary structure data.
        """
        self.__secstruct_file = secstruct_file
        self.__left_brackets = ('(')
        self.__right_brackets = (')')
        self.__nonwc_symbols = ('.')
        self.__secstruct = self.read_rna_secstruct()
        self.__wcpairs = self.get_wcpairs()
        return None


    @property
    def secstruct_file(self):
        """RNA secondray structure file path getter

        Returns
        -------
            self.__secstruct_file : str
                Path to the RNA secondary structure text file.
        """
        return self.__secstruct_file


    @property
    def secstruct(self):
        """RNA secondary structure getter

        Parameters
        ----------
            self : RNASecStructContent
                Instance of RNASecStructContent class

        Returns
        -------
            self.__secstruct : tuple
                A tuple of RNA secondary structure as read from input file.
        """
        return self.__secstruct


    @property
    def wcpairs(self):
        """WC pairs getter.

        Parameters
        ----------
            self : SecStructContent
                An instance of SecStructContent class

        Returns
        -------
            self.__wcpairs : tuple(tuples)
        """

        return self.__wcpairs


    def get_wcpairs(self):
        """Extracts WC pairs from the entire RNA secondary structure information.

        Parameters
        ----------
            self : SecStructContent
                An instance of SecStructContent class

        Returns
        -------
            wcpairs : tuple(tuples)
                A tuple containing WC pairs as tuples.
        """
        wcpairs = self.get_wcpair_indices(self.__secstruct)
        return wcpairs


    def read_rna_secstruct(self):
        """Extracts RNA secondary structure  pairs from a dot bracket represen-
        taion of the secondary structure.
        Parameters
        ----------
            self : RNASecStructContent

        Returns
        -------
            secstruct_data : tuple
                A tuple of secondary structure information as read from file.

        """
        with open(self.__secstruct_file) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('#'): continue
                if line:
                    secstruct_str = line
                    break
        if not secstruct_str:
            logger.error('\n\tUnable to get secondary structure data from'
                ' file {}'.format(self.__secstruct_file)
            )
            raise RNASecStructContentException
        brackets = list(self.__left_brackets) + list(self.__right_brackets)
        allowed_symbols = brackets + list(self.__nonwc_symbols)
        for val in secstruct_str:
            if val not in allowed_symbols:
                logger.error('\n\t{} is invalid symbol to represent RNA secondary'
                    ' structure'.format(val)
                )
                raise RNASecStructContentException
        logger.info('\n\tSecondary structure:'
            '\n\t\t{}'.format(secstruct_str),
        )
        #put the secondary structure in a list
        secstruct_data = list()
        for val in secstruct_str : secstruct_data.append(val)
        secstruct_data = tuple(secstruct_data)
        return secstruct_data


    def get_wcpair_indices(self, secstruct_data):
        """Extracts WC pairs from the secondary structure data. The pairing is
        done by taking the indices of WC pairs.

        Parameters
        ----------
            self : SecStructContent
                An instancce of SecStructContent class

            secstruct_data : list/tuple
                A list or tuple containing RNA secondary structure data. This
                data is expected to containing only symbols that are in
                self.__left_brackets, self.__right_brackets and self.__nonwc_symbols.

        Returns
        -------
            secstruct_pairs : tuple(tuples)
                A tuple containing WC pairs. The pairs are represented using
                their index position in the secstruct_data parameter. The indexing
                starts from zero.
        """
        left_bracket_indx = list()
        secstruct_pairs = list()
        for k, symbol in enumerate(secstruct_data, start=0):
            if symbol in self.__left_brackets: left_bracket_indx.append(k)
            if symbol in self.__right_brackets:
                try:
                    wc_pairs = (left_bracket_indx[-1], k)
                    secstruct_pairs.append(wc_pairs)
                    #remove the paired left bracket
                    del left_bracket_indx[-1]
                except IndexError:
                    logger.error('\n\tInvalid RNA secondary structure.'
                        ' Make sure that there are no missing brackets.'
                    )
                    raise
        if left_bracket_indx:
            logger.error('\n\tInvalid RNA secondary structure. Make sure that there'
                ' are no missing brackets.'
            )
            raise ValueError
        secstruct_pairs.sort(key=lambda x: x[0])
        logger.info('\n\tNumber of RNA secondary structure pairs: {}'.format(
            len(secstruct_pairs)),
        )
        secstruct_pairs = tuple(secstruct_pairs)
        return secstruct_pairs


class DCAContentException(Exception):
    """Used for raising exceptions related to data obtained from DCA file
    """


class DCAContent:
    """Stores data obtained from DCA compuation output files. The DCA file
    is assumed to contain site pairs in the first and second column. The DCA
    score can be omitted, but site pairs are assumed to be ranked according to
    DCA score in descending DCA score.

    Attributes
    ----------
        self.__dca_file : str
            Path to the DCA file
        self.__dca_ranked_pairs : list(tuples)
            A list containing tuples of DCA ranked pairs.
    """
    def __init__(self, dca_file = None, sorted_dca_scores=None):
        """Initializes DCAContent objects

        Parameters
        ----------
            dca_file : str
                Path to DCA computation output file.

        Returns
        -------
            None : None
        """
        self.__dca_file = dca_file
        if dca_file is not None:
            self.__dca_ranked_pairs = self.shift_dca_ranked_pair_indices()
        elif sorted_dca_scores is not None:
            self.__dca_ranked_pairs = self.get_dca_pairs_from_score_list(sorted_dca_scores)
        if dca_file is None and sorted_dca_scores is None:
            logger.error('\n\tPlease provide a DCA file or a list of ranked site pairs')
            raise DCAContentException
        self.__num_dca_ranked_pairs = len(self.__dca_ranked_pairs)
        return None


    @property
    def dca_ranked_pairs(self):
        """Property for accessing DCA ranked site pairs.

        Parameters
        ----------
            self : DCAContent

        Returns
        -------
            tuple(self.__dca_ranked_pairs) :  tuple(tuples)
                A tuple of DCA ranked pair tuples.
        """
        return tuple(self.__dca_ranked_pairs)


    @property
    def num_dca_ranked_pairs(self):
        """Property of accessing the total number of DCA ranked pairs as obtained
        from DCA file.

        Parameters
        ----------
            self : DCAContent

        Returns
        -------
            self.__num_dca_ranked_paris
                Total number of DCA ranked site pairs as read from DCA computation
                output file.
        """

        return self.__num_dca_ranked_pairs


    def get_dca_pairs_from_score_list(self, sorted_dca_scores):
        """Obtains DCA ranked site pairs from a list of tuples site-pair and score.

        Parameters
        ----------
            self : DCAContent
                An instance of DCAContent class
            sorted_dca_scores : list of tuples
                A lisf to tuple of site pair and scores sorted by score in reverse
                order
        Returns
        -------
            ranked_site_pairs : list 
                A lisf of site pairs 
        """
        ranked_site_pairs = [pair for pair, score in sorted_dca_scores]
        return ranked_site_pairs

    def shift_dca_ranked_pair_indices(self):
        """DCA output files contain pairs of sites that are labeled starting
        from 1. Here they are shifted so that they are counted from 0 so as to
        be cosistent with Python's index counting.

        Parameters
        ----------
            self : DCAContent

        Returns
        -------
            shifted_dca_ranked_pairs : list(tuples)
                A list containing DCA pairs as tuples. The pairs are shfted by
                one as opposed to they are writen in the DCA results file.
        """
        dca_ranked_pairs = self.read_dca_ranked_pairs()
        shifted_dca_ranked_pairs = [
            (pair[0] - 1, pair[1] - 1) for pair in dca_ranked_pairs
        ]

        if not all(pair[0] >= 0 for pair in shifted_dca_ranked_pairs):
            logger.error('\n\tFound negative value from DCA ranked pairs')
            raise DCAContentException
        if not all(pair[1] >= 0 for pair in shifted_dca_ranked_pairs):
            logger.error('\n\tFound negative value from DCA ranked pairs')
            raise DCAContentException
        return shifted_dca_ranked_pairs


    def read_dca_ranked_pairs(self):
        """Reads DCA pairs from DCA score file

        Parameters
        ----------
            self : DCAContent

        Returns
        -------
            dca_ranked_pairs : list(tuples)
                A list containing tuples of DCA pairs as read from DCA computation
                output file.
        """
        dca_ranked_pairs = list()
        logger.info('\n\tReading ranked residue pairs from DCA file {}'.format(
            self.__dca_file)
        )
        with open(self.__dca_file) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('#'): continue
                line = line.split()
                res_pairs = (int(line[0]), int(line[1]))
                dca_ranked_pairs.append(res_pairs)
        logger.info('\n\tNumber of DCA ranked pairs found: {}'.format(
            len(dca_ranked_pairs))
        )
        return dca_ranked_pairs



class DCAVisualizerException(Exception):
    """Used for raising exceptions related to the DCAVisualizer class.
    """

class DCAVisualizer:
    """Used to visualize contact map of a PDB structure.

    Attributes
    ----------
        self.__biomolecule : str
            The biomolecule that the PDB structure is intended to represent. This
            should reflect the type of PDB chain intended to be used. For example,
            if a PDB file contains both protein and RNA chains,  the biomolecule
            chosen should be compliant with the chain biomolecule type.

        self.__pdb_file : str
            Path to PDD file

        self.__pdb_content : PDBContent
            An instance of PDBContent class
    """

    def __init__(self, biomolecule, pdb_chain_id, pdb_file, refseq_file = None,
            dca_file = None, sorted_dca_scores = None, rna_secstruct_file = None, 
            linear_dist = None, contact_dist = None, num_dca_contacts = None, 
            wc_neighbor_dist = None,
            pdb_id = None):
        """Initializes ContactVisualizer instance with attributes:
        self.__pdb_content and self.__refseq_content of type PDBContent and
        RefSeqContent, respectively.

        Parameters
        ----------
            pdb_file : str
                Path to the PDB file
            biomolecule : str
                Type of biomolecule (protein or RNA), case insensitive.
        """
        self.__biomolecule = biomolecule.strip().upper()
        if self.__biomolecule not in ('PROTEIN', 'RNA'):
            logger.error('\n\t{} is invalid input. Must be protein (PROTEIN) or'
                ' rna (RNA)'.format(self.__biomolecule),
            )
            raise DCAVisualizerException
        self.__pdb_file = pdb_file
        self.__pdb_content = PDBContent(self.__pdb_file,
            biomolecule = self.__biomolecule,
        )
        self.__pdb_chain_id = pdb_chain_id.strip().upper()
        # Distance between residues in sequence
        if linear_dist is not None:
            self.__linear_dist = linear_dist
        else:
            self.__linear_dist = 4
        if self.__linear_dist < 0:
            logger.error(
                '\n\tDistance between two residues in a sequence cannot be negative',
            )
            raise DCAVisualizerException
        # Distance between two residues to be considered contacts in PDB structure.
        if contact_dist is not None:
            self.__contact_dist = contact_dist
        else:
            self.__contact_dist = 8.0
        if self.__contact_dist < 0:
            logger.error('\n\tPDB contact distance cannot be negative')
            raise DCAVisualizerException
        #initialize a RefSeqContent
        if refseq_file is not None:
            self.__refseq_content = RefSeqContent(refseq_file,
                biomolecule = self.__biomolecule,
            )
        else:
            self.__refseq_content = None
        #initialize a DCAContent
        if dca_file is not None:
            self.__dca_content = DCAContent(dca_file=dca_file)
        elif sorted_dca_scores is not None:
            self.__dca_content = DCAContent(sorted_dca_scores=sorted_dca_scores)
        else:
            self.__dca_content = None
        #initialize a SecStructContent
        if rna_secstruct_file is not None:
            self.__rna_secstruct_content = RNASecStructContent(rna_secstruct_file)
            if wc_neighbor_dist is not None:
                self.__wc_neighbor_dist = wc_neighbor_dist
            else:
                self.__wc_neighbor_dist = 0
            if self.__wc_neighbor_dist < 0:
                logger.error('\n\tWC neighbour distance cannot be negative.')
                raise DCAVisualizerException
        else:
            self.__rna_secstruct_content = None
            self.__wc_neighbor_dist = None
            if self.__biomolecule == 'RNA':
                logger.warning('\n\tYou didn\'t supply RNA secondary structure file.')
        self.__refseq_len = len(self.get_matching_refseq_to_biomolecule())
        # initialize number of top DCA ranked contacts to be taken
        if num_dca_contacts is None:
            self.__num_dca_contacts = self.__refseq_len
        else:
            if num_dca_contacts <= self.__dca_content.num_dca_ranked_pairs:
                self.__num_dca_contacts = num_dca_contacts
            else:
                logger.error('\n\tThe maximum number of DCA ranked site pairs'
                    ' found in DCA file is {}. But you supplied {}'.format(
                        self.__dca_content.num_dca_ranked_pairs,
                        num_dca_contacts,
                    )
                )
                raise DCAVisualizerException
            if self.__num_dca_contacts < 0:
                logger.error('\n\tThe number of DCA contacts cannot be negative')
                raise DCAVisualizerException
        self.__pdb_id = pdb_id if pdb_id is not None else None

        # make sure that the reference sequence and secondary structure data
        # have equal lengths in the case or RNAs
        if self.__biomolecule == 'RNA':
            if self.__refseq_content and self.__rna_secstruct_content:
                secstruct_len = len(self.__rna_secstruct_content.secstruct)
                if self.__refseq_len != secstruct_len:
                    logger.error('\n\tThe length of RNA secondary structure and '
                        ' reference sequence are not equal.'
                        '\n\treference sequence length: {}'
                        '\n\tsecondary structure length: {}'
                        '\n\tDid you give a wrong input file?'.format(
                            self.__refseq_len, secstruct_len)
                    )
                    raise DCAVisualizerException
        logg_mssg = '''
            \tCreated {} with:
            \tBiomolecule: {}
            \tPDB chain ID: {}
            \tMinimum linear distance between residues: {}
            \tMaximum contact distance : {}
            \tReference sequence length: {}'''

        logg_mssg = logg_mssg.format( __class__.__name__, self.__biomolecule,
            self.__pdb_chain_id, self.__linear_dist, self.__contact_dist,
            self.__refseq_len,
        )
        if self.__biomolecule == 'RNA':
            rna_logg_add = '''
            \tRNA secstruct file: {}
            \tWC neighbour distance : {}
            '''.format(rna_secstruct_file, self.__wc_neighbor_dist)
            logg_mssg = logg_mssg + rna_logg_add
        logger.info(logg_mssg)
        return None

    @property
    def biomolecule(self):
        """
        """
        return self.__biomolecule

    @property
    def contact_dist(self):
        """
        """
        return self.__contact_dist


    @property
    def linear_dist(self):
        """
        """
        return self.__linear_dist


    @property
    def wc_neighbor_dist(self):
        """
        """
        return self.__wc_neighbor_dist


    @property
    def pdb_id(self):
        """
        """
        return self.__pdb_id


    @property
    def pdb_chain_id(self):
        """
        """
        return self.__pdb_chain_id


    @property
    def pdb_content(self):
        """
        """
        return self.__pdb_content


    @property
    def refseq_content(self):
        """
        """
        return self.__refseq_content


    @property
    def rna_secstruct_content(self):
        """
        """
        return self.__rna_secstruct_content


    @property
    def dca_content(self):
        """
        """
        return self.__dca_content


    def get_matching_refseq_to_biomolecule(self):
        """Obtains the reference sequence that matches the chains biomolecule type

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            refseq : str
                The reference sequence that has a biomolecule type of
                self.__biomolecule.
        """
        refseq = None
        for key, val in self.__refseq_content.ref_sequences.items():
            if self.__biomolecule == val[0]:
                refseq = val[1]
                break
        if refseq is None:
            logger.error('\n\tUnable to find a reference sequence that matches'
                ' biomolecule type {}'.format(self.__biomolecule)
            )
            raise DCAVisualizerException
        return refseq


    def align_refseq_and_pdbseq(self):
        """Align the standared residue sequence obtained from PDB file with that
        of the pdb sequence in the corresponding PDB fasta file (reference sequence file).

        Parameters
        ----------
            self : DCAVisualizer
                Instance of DCAVisualizer class

        Returns
        -------
            alignment : list(tuples)
                A list containing tuples of alignmed sequences, score, start and end indices
                of the pairwise alignment. The content of the tuples is as follows:
                (refseq_aligned, pdbseq_aligned, score, start_index, end_index). If there
                are multiple high scoring pairwise alignments, the list will have multiple
                aforementioned tuples. Typically, the reference and PDB chain sequences
                differ by a few missing residues at the begining and end of the PDB chain
                sequence.
        """

        ref_seq = self.get_matching_refseq_to_biomolecule()
        biomol_info, pdb_seq = self.__pdb_content.pdb_chain_sequences[self.__pdb_chain_id]
        if len(ref_seq) < len(pdb_seq):
            logger.warning('\n\tThe reference sequence contains less number of'
            ' residues compared to PDB chain sequence. \n\tMake sure that'
            ' you have provided correct input data.'
            )
        logger.info('\n\tReference sequence:\n\t{}'.format(ref_seq))
        logger.info('\n\tSequence from PDB chain {}:\n\t{}'.format(
            self.__pdb_chain_id, pdb_seq)
        )

        if self.__biomolecule == 'PROTEIN':
            scoring_mat = blosum62
            GAP_OPEN_PEN = -10
            GAP_EXTEND_PEN = -1
        if self.__biomolecule == 'RNA':
            scoring_mat = scoring_matrix.NUC44
            GAP_OPEN_PEN = -8
            GAP_EXTEND_PEN = 0
        alignment = pairwise2.align.localds(
            ref_seq,
            pdb_seq,
            scoring_mat,
            GAP_OPEN_PEN,
            GAP_EXTEND_PEN,
            score_only=False,
        )
        seq_ref_aligned = alignment[0][0]
        seq_pdb_aligned = alignment[0][1]
        alignment_starts_at = alignment[0][3]
        alignment_ends_at = alignment[0][4]
        alignment_score = alignment[0][2]
        logger.info('\n\tPDB fasta seq. and residue seq., respectively, aligned:'
            '\n\t{}'
            '\n\t{}'
            '\n\tAlignment score: {}'
            '\n\tAlignment starting and ending positions: {}'.format(
                seq_ref_aligned,
                seq_pdb_aligned,
                alignment_score,
                [alignment_starts_at + 1, alignment_ends_at], # biopython's pairwise2 returns [start, end), counting from 0.
            )
        )

        seq_pdb_aligned_subseq = seq_pdb_aligned[alignment_starts_at:alignment_ends_at]
        gap_char = '-'
        if gap_char in seq_pdb_aligned_subseq:
            logger.warning(
                '\n\tPDB sequence has got gaps at the middle when aligned with reference.'
                '\n\tUsually residues are missing at the beginning and end of the chain.'
            )
        return alignment


    def map_pdbseq_to_refseq(self):
        """Mappes PDB chain sequence residues to that of reference sequence.
        The reference and PDB chain sequence to be mapped are the
        results of pairwise alignment of these sequences.

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            mapped_res_pos : OrderedDict
                A dictionary  whose keys are the postion (site) of the pairwise
                aligned PDB chain sequence and whose values are the corresponding
                aligned reference sequence residue positions.
            residues_not_found_in_pdb : list
                A list of sites that are in the reference sequence but not in
                the PDB chain sequence.
        """
        # Take the first alignment result beteween reference sequence and PDB
        # chain sequence
        alignment = self.align_refseq_and_pdbseq()[0]
        logger.info('\n\tMapping PDB chain residues to reference sequence residues')
        refseq_aligned = alignment[0]
        pdbseq_aligned = alignment[1]
        score = alignment[2]
        start_indx = alignment[3]
        end_indx = alignment[4]
        ref_res_pos = -1
        pdb_res_pos = -1
        mapped_res_pos = OrderedDict()
        residues_not_found_in_pdb = list()
        gap_char = '-'
        for res_ref, res_pdb in zip(refseq_aligned, pdbseq_aligned):
            if res_ref != gap_char: ref_res_pos += 1
            if res_pdb != gap_char: pdb_res_pos += 1
            if res_ref != gap_char and res_pdb == gap_char:
                residues_not_found_in_pdb.append(ref_res_pos)
            if res_ref !=gap_char and res_pdb != gap_char:
                mapped_res_pos[pdb_res_pos] = ref_res_pos
        logger.info('\n\tTotal number of PDB chain residues mapped: {}'.format(
            pdb_res_pos + 1)
        )
        residues_not_found_logged = [res + 1 for res in residues_not_found_in_pdb]
        logger.info(
            '\n\tResidue sites that are not found in PDB chain: {}'.format(
                residues_not_found_logged
            )
        )
        return mapped_res_pos, residues_not_found_in_pdb


    def get_mapped_pdb_contacts(self):
        """Collects all resisdue contacts that are within a cut-off distance
        as specified by self.__contact_dist parameter

        Parameters
        ----------
            self : DCAVisualizer
                An instance of DCAVisualizer class.
        Returns
        -------
            mapped_residues : dict
                A dictionary whose keys are mapped sites of pairs in reference
                sequence and whose values are a tuple of metadata about the
                corresponding sites in PDB standared chain sequence. The
                metadata is closest atom name pair, residue id (site-number only)
                in the PDB and minimum atom pair distance.
            residues_not_found_in_pdb : list
                List of reference sequence sites that are not found in PDB chain.
                This list comes from map_pdbseq_to_refseq() method.
        """

        try:
            chain_biomolecule = self.__pdb_content.pdb_chain_sequences[self.__pdb_chain_id][0]
        except KeyError:
            logger.error('\n\tUnable to find a PDB chain of ID {}'
                '\n\tFrom PDB file {}'.format(self.__pdb_chain_id, self.__pdb_file)
            )
            raise
        else:

            if self.__biomolecule != chain_biomolecule:
                logger.error("\n\tChain {} does not contain {} residues.".format(
                    self.__pdb_chain_id, self.__biomolecule)
                )
                raise DCAVisualizerException
        logger.info('\n\tObtaining PDB contacts for {} PDB chain {}'
            ' using cut-off distance of {} Angstrom'.format(self.__biomolecule,
                self.__pdb_chain_id, self.__contact_dist
            )
        )
        structure = self.__pdb_content.pdb_structure
        model = structure[0]
        chain = model[self.__pdb_chain_id]
        all_residues = chain.get_list()
        standard_residues = self.__pdb_content.filter_residues(
            all_residues,
            self.__biomolecule,
        )
        mapping_key, residues_not_found_in_pdb = self.map_pdbseq_to_refseq()
        mapped_residues = dict()
        num_residues = len(standard_residues)
        for i in range(num_residues - 1):
            res_1 = standard_residues[i]
            for j in range(i + 1, num_residues):
                res_2 = standard_residues[j]
                min_atom_dist = 100000.0
                # find the closest atom pair for the two residues
                for atom_1, atom_2 in itertools.product(res_1, res_2):
                    atom_1_name = atom_1.get_name().strip()
                    atom_2_name = atom_2.get_name().strip()
                    if atom_1_name.startswith('H') or atom_2.name.startswith('H'):continue
                    dist = atom_1 - atom_2
                    if dist  < min_atom_dist:
                        min_atom_dist = dist
                        atom_pair = atom_1_name + '-' + atom_2_name
                try:
                    mapped_pair = (mapping_key[i], mapping_key[j])
                    metadata = (atom_pair, res_1.id[1], res_2.id[1], min_atom_dist)
                except KeyError: # this happens when not all PDB residues are mapped
                    pass
                else:
                    mapped_residues[mapped_pair] = metadata
        return mapped_residues, residues_not_found_in_pdb


    def get_wc_pairs_and_neighbors(self):
        """When RNA secondary structure file is supplied, the WC pairs are returned.
        If the WC neighbor distance is not Zero a total number of (2d + 1)(2d + 1)
        pairs are formed for each WC pair and added to the list of WC pairs forming
        a WC and their neighbors list of pairs. If there is no secondary structure
        information, empty list is returned.

        Parameters:
        -----------
            self : DCAVisualizer
                An instance of DCAVisualizer class.

        Returns
        -------
            wc_pairs_and_neighbors : list
                A tuple of pair tuples or None (if there is no secondary structure
                information).
        """
        wc_pairs_and_neighbors  = list()
        if not self.__rna_secstruct_content:
            logger.warning('\n\tRNA secondary structure information is not found,'
                ' cannot obtain WC pairs and their neighbors'
            )
            return  wc_pairs_and_neighbors
        else:
            #first pass logging message based on the value of WC neighbor distance.
            if self.__wc_neighbor_dist > 0:
                logger.info('\n\tObtaining WC pairs and their neighbors with'
                    ' distance: {}'.format(self.__wc_neighbor_dist)
                )
            elif self.__wc_neighbor_dist == 0:
                logger.warning('\n\tWC neighbor distance is Zero, obtaining WC pairs only')
            else:
                logger.error('\n\tInvalid value of WC neighbor distance')
                raise DCAVisualizerException
            # extract WC pairs (and neighbors, if neighbor dist > 0)
        wc_pairs = self.__rna_secstruct_content.wcpairs
        logger.info('\n\tNumber of WC pairs found: {}'.format(len(wc_pairs)))
        out_of_bounds_pairs = list()
        for first, second in wc_pairs:
            first_subsites = list()
            second_subsites = list()
            for i in range(-self.__wc_neighbor_dist, self.__wc_neighbor_dist + 1):
                left = first + i
                right = second + i
                if left < 0 or left >= self.__refseq_len:
                    # add + 1 for out of bounds pairs since they are used only for
                    # logging messages
                    out_of_bounds_pairs.append((left + 1, right + 1))
                    continue
                if right < 0 or right >= self.__refseq_len:
                    out_of_bounds_pairs.append((left + 1, right + 1))
                    continue
                first_subsites.append(left)
                second_subsites.append(right)
            for pair in itertools.product(first_subsites, second_subsites):
                wc_pairs_and_neighbors.append(pair)
        logger.info('\n\tOut of bound pairs: {}'.format(out_of_bounds_pairs))
        logger.info('\n\tTotal number of WC pairs and neighbors: {}'.format(
            len(wc_pairs_and_neighbors))
        )
        return wc_pairs_and_neighbors


    def select_top_dca_ranked_contacts(self, num_dca_contacts=None):# may be not necessary for now
        """Returns top N DCA ranked contacts. If the linear distance is set to
        non-zero, only tertiary contacts are returned. If the linear distance
        is zero, no filtering of tertiary contacts is carried out. For RNAs, WC
        pairs and their neighbors are filtered if RNA secondary structure file
        is supplied. To get top N top DCA ranked pairs, set the linear distance
        to zero for proteins; for RNAs in addition to setting the linear distance
        to zero, do not supply a file containing RNA secondary structure.

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            top_n_contacts : tuple
        """
        if num_dca_contacts is None: num_dca_contacts = self.__num_dca_contacts
        all_dca_contacts = self.__dca_content.dca_ranked_pairs
        if self.__biomolecule == 'RNA':
            wc_and_neighbors = self.get_wc_pairs_and_neighbors()
            if not wc_and_neighbors:
                logger.warning('\n\tWC and neighbors are not filtered from DCA contacts')
            dca_less_wc_and_neighbors = [
                p for p in all_dca_contacts if p not in wc_and_neighbors
            ]
            logger.info('\n\tRemoving DCA contacts that are more than {} residues'
                ' apart in chain'.format(self.__linear_dist)
            )
            dca_less_wc_and_neighbors_linear_dist = [
                p for p in dca_less_wc_and_neighbors if abs(p[0]-p[1]) > self.__linear_dist
            ]
            top_n_contacts = dca_less_wc_and_neighbors_linear_dist[:num_dca_contacts]
            return top_n_contacts
        elif self.__biomolecule == 'PROTEIN':
            logger.info('\n\tRemoving DCA contacts that are more than {} residues'
                ' apart in chain'.format(self.__linear_dist)
            )
            dca_less_linear_dist = [
                p for p in all_dca_contacts if abs(p[0] - p[1]) > self.__linear_dist
            ]
            top_n_contacts = dca_less_linear_dist[:num_dca_contacts]
            logger.info('\n\tObtaining top {} DCA contacts'.format(num_dca_contacts))
            return top_n_contacts
        else:
            logger.error('\n\t{} is invalid biomolecule type. Must be protein'
                ' or RNA'.format(self.__biomolecule)
            )
            raise DCAVisualizerException
    #end of select_top_dca_ranked_contacts(self, num_dca_contacts=None)

    def dca_ranked_pairs_filtered_by_linear_dist(self, num_dca_contacts =None):
        """Filters DCA ranked contacts by chain (linear distance).
        If the linear distance is set to Zero, no filtering is done.

        Parameters
        ----------
            self : DCAVisualizer
                An instance of DCAVisualizer class.

        Returns
        -------
            all_dca_pairs or filtered_dca_pairs :  tuple(tuples)
                A tuple containing tuples of all or filtered DCA ranked pairs.
        """
        # get all DCA ranked pairs
        if num_dca_contacts is None: num_dca_contacts = self.__num_dca_contacts
        all_dca_pairs = self.__dca_content.dca_ranked_pairs
        if self.__linear_dist == 0:
            logger.warning('\n\tChain distance between residues is Zero,'
                '\n\treturning top {} unfiltered DCA contacts'.format(num_dca_contacts)
            )
            return tuple(all_dca_pairs[:num_dca_contacts])
        elif self.__linear_dist > 0:
            logger.info('\n\tFiltering DCA ranked pairs by chain distance of :'
                ' {}'.format(self.__linear_dist)
            )
            # get filtered DCA ranked pairs as filtered by chain (linear) distance
            filtered_dca_pairs = [
                p for p in all_dca_pairs if abs(p[0] - p[1]) > self.__linear_dist
            ]
            return tuple(filtered_dca_pairs[:num_dca_contacts])
        else:
            logger.error('\n\tInvalid value of linear distance: {}'.format(
                self.__linear_dist)
            )
            raise DCAVisualizerException
    #### end of def dca_ranked_pairs_filtered_by_linear_dist(self)

    def split_and_shift_contact_pairs(self, list_of_contacts):
        """Separating contacting site pairs into two lists: One containing all
        first entries (xdata), and the other containing all second entries (ydata).

        Note that, each site pair is shifted by 1 so as to make them ready for
        visualization (since the sites have been counted starting from 0)

        Parameters
        ----------
            self : DCAVisualizer

            list_of_contacts : list/tuple
                A list or tuple containing tuples of contacting site pairs.

        Returns
        -------
            xdata : list
                List containing shifted first entries of site pairs as obtained
                from list_of_contacts parameter.
            ydata : list
                List containing shifted second entries of site pairs as obtained
                from list_of_contacts parameter
        """
        xdata = list()
        ydata = list()
        for first, second in list_of_contacts:
            xdata.append(first + 1) # shift indexing by 1 for visualization (output)
            ydata.append(second + 1)

        return xdata, ydata


    def contact_categories(self):
        """Categorizes contacts to true positives, false positives, missing in
        PDB or PDB contacts.

        Parameters
        ----------
            self : DCAVisualizer
                An instance of DCAVisualizer class

        Returns
        -------
            contact_categories_dict : dict
                A dictionary consisting of contact categories. The keys
                are tp for true positives, fp for false positives, missing for
                those misssing in PDB but are part of the DCA prediction or pdb
                for those that are obtained by mapped reference sequence to the
                corresponding PDB chain.
        """
        mapped_pdb_contacts, missing_residues = self.get_mapped_pdb_contacts()
        # take the top L (refseq_len) DCA ranked pairs
        top_dca_ranked_pairs = self.dca_ranked_pairs_filtered_by_linear_dist()

        logger.info('\n\tCategorizing contacts')
        logger.info('\n\tTaking top {} DCA ranked pairs for contact'
            ' map comparison'.format(len(top_dca_ranked_pairs))
        )
        #filtered_mapped_pdb_contacts = [
        #    p for p in mapped_pdb_contacts if abs(p[0] - p[1]) > self.__linear_dist
        #]
        missing_dca_contacts = list()
        if missing_residues: # if there are missing residues in PDB chain
            for  pair in top_dca_ranked_pairs:
                if pair[0] in missing_residues or pair[1] in missing_residues:
                    missing_dca_contacts.append(pair)
        # for contact map plot we need only contacting PDB residue pairs.
        contacts_in_pdb = OrderedDict()
        for pair in mapped_pdb_contacts:
            metadata = mapped_pdb_contacts[pair]
            if metadata[-1] < self.__contact_dist:
                contacts_in_pdb[pair] = metadata
        # get true positives and false positives among all top DCA ranked pairs
        true_positives = OrderedDict()
        false_positives = OrderedDict()
        for p1 in top_dca_ranked_pairs:
            try:
                metadata = mapped_pdb_contacts[p1]
            except KeyError: # if residues are missing in PDB, they had not been mapped
                pass
            else:
                if metadata[-1] < self.__contact_dist:
                    true_positives[p1] = metadata
                else:
                    if p1 not in missing_dca_contacts:
                        false_positives[p1] = metadata
        missing_dca_contacts_filtered = OrderedDict()
        for pair in missing_dca_contacts:
            if abs(pair[0] - pair[1]) > self.__linear_dist:
                missing_dca_contacts_filtered[pair] = pair # no pdb metadata for missing residues
        contact_categories_dict = dict()
        contact_categories_dict['tp'] = true_positives
        contact_categories_dict['fp'] = false_positives
        contact_categories_dict['missing'] = missing_dca_contacts_filtered
        contact_categories_dict['pdb'] = contacts_in_pdb
        return contact_categories_dict


    def _plot_contact_map_rna(self):
        """Plots RNA contact map and returns the contact categories dict.

        Parameters
        ----------
            self : DCAVisualizer
                An instance of DCAVisualizer class.

        Returns
        -------
            contact_categories_dict : dict
                A dictionary of contact categories.
        """
        contact_categories_dict = self.contact_categories()
        true_positives = contact_categories_dict['tp']
        false_positives = contact_categories_dict['fp']
        missing_dca_contacts  = contact_categories_dict['missing']
        mapped_pdb_contacts = contact_categories_dict['pdb']
        # Raise exception if the requested number of DCA predicted site-pairs is
        # greater than that of the number of contacts in PDB 
        filtered_pdb_contacts_list = [ 
           site_pair for site_pair, metadata in mapped_pdb_contacts.items() if abs(site_pair[1] - site_pair[0]) > self.__linear_dist  
        ]
        num_filtered_pdb_contacts = len(filtered_pdb_contacts_list)
        if  self.__num_dca_contacts > num_filtered_pdb_contacts:
            logger.error('\n\tMaximum number of PDB contacts with linear distance {} is {}.' 
                '\n\tPlease set the number of DCA contacts a maximum of this value'.format(
                    self.__linear_dist, 
                    num_filtered_pdb_contacts
                ) 
            )
            raise DCAVisualizerException

        x_false_positive, y_false_positive = self.split_and_shift_contact_pairs(
            false_positives
        )

        x_pdb_contacts, y_pdb_contacts = self.split_and_shift_contact_pairs(
            mapped_pdb_contacts
        )

        num_contacts_compared = len(true_positives) + len(false_positives)
        fraction_of_true_positives = float(len(true_positives))/float(num_contacts_compared)
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,5))
        if missing_dca_contacts:
            x_missing, y_missing = self.split_and_shift_contact_pairs(
                missing_dca_contacts
            )
            ax.scatter(y_missing, x_missing,
                s=6, color='blue',label='missing in PDB'
            )
            missing_pairs_to_be_logged = [
                (first, second) for first, second in zip(x_missing, y_missing)
            ]
            logger.info('\n\tPairs in DCA ranked sites but not found in PDB'
                ' chain site pairs: {}'
                '\n\tNote that these pairs are neither TP nor FP'
                ''.format(missing_pairs_to_be_logged)
            )
        #

        ax.scatter(x_pdb_contacts, y_pdb_contacts,
            s=6, color='grey', label='PDB contacts (PDB ID : {})'.format(self.__pdb_id),
        )
        ax.scatter(y_false_positive, x_false_positive,
            s=6, color='red', label='false positives',
        )
        # set plot attributes and display contact map comparison plot.
        ax_title = '''
        Maximum PDB contact distance : {} Angstrom
        Minimum residue chain distance: {} residues
        Number of DCA contacts : {}
        Fraction of true positives : {:.3g}
        '''.format(self.__contact_dist, self.__linear_dist,
            self.__num_dca_contacts,
            fraction_of_true_positives,
        )
        if self.__rna_secstruct_content:
            wc_pairs = self.__rna_secstruct_content.wcpairs
            top_dca_ranked_pairs = OrderedDict(
                list(true_positives.items()) + list(false_positives.items())
            )
            predicted_wc_pairs = OrderedDict()
            for pair in top_dca_ranked_pairs:
                if pair in wc_pairs:
                    predicted_wc_pairs[pair] = top_dca_ranked_pairs[pair]
            predicted_non_wc_pairs = OrderedDict()
            for pair in top_dca_ranked_pairs:
                if pair not in predicted_wc_pairs:
                    predicted_non_wc_pairs[pair] = top_dca_ranked_pairs[pair]
            predicted_non_wc_true_positives = OrderedDict()
            for pair in predicted_non_wc_pairs:
                if pair not in false_positives:
                    predicted_non_wc_true_positives[pair] = predicted_non_wc_pairs[pair]
            logger.info('\n\tCategorizing true positives as WC and non-WC DCA contacts')
            contact_categories_dict['tp-wc'] = predicted_wc_pairs
            contact_categories_dict['tp-nwc'] = predicted_non_wc_true_positives
            logger.info('\n\tRemoving true positives from contact categories.'
                '\n\t(True positives are already categorized into WC and non-WC)'
            )
            contact_categories_dict.pop('tp', None)
            wc_first, wc_second = self.split_and_shift_contact_pairs(predicted_wc_pairs)
            logger.info('\n\tNumber of correctly predicted WC pairs:{}'
                ''.format(len(predicted_wc_pairs))
            )
            tp_non_wc_first, tp_non_wc_second = self.split_and_shift_contact_pairs(
                predicted_non_wc_true_positives
            )
            ax.scatter(tp_non_wc_second, tp_non_wc_first, s=6, color='green',
                label='predicted Non-WC contacts',
            )
            ax.scatter(wc_second, wc_first, s=6, color='black',
                label='predicted WC contacts',
            )
            ax_title = ax_title + 'Correctly predicted WC pairs : {}'.format(
                len(predicted_wc_pairs)
            )
            ax_title = ax_title + '\nCorrectly predicted non-WC pairs: {}'.format(
                len(predicted_non_wc_pairs) -len(false_positives)
            )
        else: # if secondary structure information is not given
            #
            x_true_positive, y_true_positive = self.split_and_shift_contact_pairs(
                true_positives
            )
            ax.scatter(y_true_positive, x_true_positive, s=6, color='green',
                label = 'true positives'
            )

        ax.set_title(ax_title)
        ax.set_xlabel('residue position', fontsize=14)
        ax.set_ylabel('residue position', fontsize=14)
        #plt.legend()
        #fig.set_dpi(600)
        plt.tight_layout()
        pdb_file_basename, ext  = os.path.splitext(os.path.basename(self.__pdb_file))
        figure_name = 'contact_map_'+ pdb_file_basename +'.png'
        #plt.savefig(figure_name, dpi=600)
        plt.show()
        return contact_categories_dict


    def _plot_contact_map_protein(self):
        """Plots contact map of a protein seqeunce. DCA predicted contacts are
        categorized as tp, fp, and missing to denote true positives, false
        positives and those missing in PDB, respectively.

        Parameters
        ----------
            self : DCAVisualizer
                An instance of DCAVisualizer class.

        Returns
        -------
            contact_categoreis_dict : dict
                A dictionary of categorized contacts.
        """
        contact_categories_dict = self.contact_categories()
        true_positives = contact_categories_dict['tp']
        false_positives = contact_categories_dict['fp']
        missing_pairs = contact_categories_dict['missing']
        pdb_contacts =  contact_categories_dict['pdb']

        filtered_pdb_contacts_list = [ 
           site_pair for site_pair, metadata in pdb_contacts.items() if abs(site_pair[1] - site_pair[0]) > self.__linear_dist  
        ]
        num_filtered_pdb_contacts = len(filtered_pdb_contacts_list)
        if  self.__num_dca_contacts > num_filtered_pdb_contacts:
            logger.error('\n\tMaximum number of PDB contacts with linear distance {} is {}.' 
                '\n\tPlease set the number of DCA contacts a maximum of this value'.format(
                    self.__linear_dist, 
                    num_filtered_pdb_contacts
                ) 
            )
            raise DCAVisualizerException

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,5))
        if missing_pairs:
            x_missing, y_missing = self.split_and_shift_contact_pairs(missing_pairs)
            ax.scatter(x_missing, y_missing, s=6, color='blue')
        x_true_positives, y_true_positives = self.split_and_shift_contact_pairs(
            true_positives,
        )
        x_false_positives, y_false_positives = self.split_and_shift_contact_pairs(
            false_positives
        )
        x_pdb, y_pdb = self.split_and_shift_contact_pairs(pdb_contacts)
        ax_title = '''
        Maximum PDB contact distance : {} Angstrom
        Minimum residue chain distance: {} residues
        Number of DCA contacts : {}
        Fraction of true positives : {:.3g}
        '''.format(self.__contact_dist, self.__linear_dist,
            self.__num_dca_contacts,
            len(true_positives)/(len(true_positives) + len(false_positives)),
        )

        ax.scatter(y_true_positives, x_true_positives, s=6, color='green')
        ax.scatter(y_false_positives, x_false_positives, s=6, color='red')
        ax.scatter(x_pdb, y_pdb, s=6, color='grey')
        ax.set_xlabel('residue position', fontsize=14)
        ax.set_ylabel('residue position', fontsize=14)
        ax.set_title(ax_title)
        plt.tight_layout()
        plt.show()
        return  contact_categories_dict


    def plot_contact_map(self):
        """Plots the contact map of an RNA or protein by calling the respective
        methods.

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            contact_categories_dict : dict
                A dictionary whose keys are contact types (dca, pdb, missing)
            and whose values are OrderedDict or list (if missing in PDB and not
            mapped).
        """
        if self.__biomolecule == 'RNA':
            contact_categories_dict = self._plot_contact_map_rna()
            return contact_categories_dict
        elif self.__biomolecule == 'PROTEIN':
            contact_categories_dict = self._plot_contact_map_protein()
        else:
            logger.error('\n\tCannot plot contact map for biomolecule {}'.format(
                self.__biomolecule)
            )
            raise DCAVisualizerException
        return contact_categories_dict
    ## end of def plot_contact_map(self)

    def compute_true_positive_rates(self):
        """Computes the true positive rate of DCA predicted contacts by comparing
        with those of in a PDB chain.

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            true_positive_rates_dict : dict
                A dictionary whose keys are the types of true positive rates (
                e.g., dca, pdb) and whose values are a list of true positive
                rates per rank.
        """
        max_num_dca_contacts = int(0.5 * self.__refseq_len*(self.__refseq_len))
        all_filtered_dca_contacts = self.dca_ranked_pairs_filtered_by_linear_dist(
            num_dca_contacts = max_num_dca_contacts # this number is larger than
            # the expected DCA pairs but we want all pairs and python truncates the list
            # when out of bound slice is requested.
        )
        pdb_content, missing_pairs = self.get_mapped_pdb_contacts()
        #remove missing contacts if from DCA contacts, since there is no way
        # to compare them for TP rate.
        dca_contacts_less_missing = [
            pair for pair in all_filtered_dca_contacts if pair not in missing_pairs
        ]
        logger.info('\n\tNumber of filtered DCA contacts after removal of possibly'
            '\n\tmissing pairs: {}'.format(len(dca_contacts_less_missing))
        )
        # mapped PDB contacts contain all contacts, we need to filter them out.
        filtered_pdb_contacts = OrderedDict()
        for pair in pdb_content:
            if abs(pair[0] - pair[1]) > self.__linear_dist:
                if pdb_content[pair][3] < self.__contact_dist:
                    filtered_pdb_contacts[pair] = pdb_content[pair]

        num_pdb_contacts = len(filtered_pdb_contacts)
        logger.info('\n\tNumber of PDB contacts: {}'.format(
            len(filtered_pdb_contacts),
            )
        )
        num_tps = 0
        dca_tp_rates = list()
        pdb_tp_rates = list()
        for counter, dca_pair in enumerate(all_filtered_dca_contacts, start=1):
            if dca_pair in filtered_pdb_contacts:
                num_tps += 1
            current_dca_tp_rate = float(num_tps)/float(counter)
            dca_tp_rates.append(current_dca_tp_rate)
            if counter <= num_pdb_contacts:
                pdb_tp_rates.append(1.0)
            else:
                current_pdb_tp_rate = float(num_pdb_contacts)/float(counter)
                pdb_tp_rates.append(current_pdb_tp_rate)
        true_positive_rates_dict = dict()
        true_positive_rates_dict['dca'] = dca_tp_rates
        true_positive_rates_dict['pdb'] = pdb_tp_rates
        return true_positive_rates_dict


    def plot_true_positive_rates(self):
        """Plotes the true positive rate per rank of DCA ranked site pairs.
        The x-axis is in log scale.

        Parameters
        ----------
            self : DCAVisualizer

        Returns
        -------
            true_positive_rates_dict : dict
                A dictionary whose keys are true positives types (pdb, or dca)
                and whose values are the corresponding true positive rates per
                rank.

        """
        true_positive_rates_dict = self.compute_true_positive_rates()
        dca_true_positive_rates = true_positive_rates_dict['dca']
        pdb_true_positive_rates = true_positive_rates_dict['pdb']
        max_rank = len(dca_true_positive_rates)
        ranks = [i + 1 for i in range(max_rank)]
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
        ax.plot(ranks, dca_true_positive_rates)
        ax.plot(ranks, pdb_true_positive_rates)
        ax.set_xscale('log')
        ax_title = '''
        True Positive Rate Per Rank
        PDB cut-off distance : {} Angstrom
        Residue chain distance : {}
        '''
        if self.__biomolecule == 'RNA':
            add_ax_title = '''WC neighbour distance : {}
            '''.format(self.__wc_neighbor_dist)
            ax_title = ax_title + add_ax_title
        ax.set_title(ax_title.format(
                self.__contact_dist, self.__linear_dist,
            )
        )
        ax.set_xlabel('rank (log scalled)', fontsize=14)
        ax.set_ylabel('true positives/rank', fontsize=14)
        plt.grid()
        plt.tight_layout()
        plt.show()
        return true_positive_rates_dict


if __name__ == '__main__':
    from argparse import ArgumentParser
    import logging.config
    import sys
    from config_log import LOGGING_CONFIG
    from config_log import ConsoleColor as ccolor

    def configure_logging():
        """
        """
        logging.config.dictConfig(LOGGING_CONFIG)
        logging.addLevelName(
            logging.INFO,
            '{}{}{}'.format(
                ccolor.green,
                logging.getLevelName(logging.INFO),
                ccolor.nocolor
            )
        )

        logging.addLevelName(
            logging.WARNING,
            '{}{}{}'.format(
                ccolor.yellow,
                logging.getLevelName(logging.WARNING),
                ccolor.nocolor
            )
        )

        logging.addLevelName(
            logging.ERROR,
            '{}{}{}'.format(
                ccolor.red,
                logging.getLevelName(logging.ERROR),
                ccolor.nocolor
            )
        )
        return None

    parser = ArgumentParser(description = "PDB analyzer parser")
    subparser = parser.add_subparsers(dest='subcommand')
    pdb_parser = subparser.add_parser(
        'pdb_content',
        help = 'subcommand for analysing contents of a PDB file',
    )
    refseq_parser = subparser.add_parser(
        'refseq_content',
        help = 'subcommand for analyzing contents of a reference sequence file',
    )

    dcavisualizer = subparser.add_parser(
        'dcavisualizer',
        help = 'subcommand for visualizing residue-residue contacts',
    )

    dcavisualizer.add_argument('biomolecule',
        help = 'Type of biomolecule PDB chain represents',
        choices = ('protein', 'PROTEIN', 'rna', 'RNA'),
    )
    dcavisualizer.add_argument('pdb_chain_id',
        help = 'The ID of PDB chain',
    )
    dcavisualizer.add_argument('pdb_file', help = 'path to a PDB file')
    dcavisualizer.add_argument('refseq_file', help = 'path to reference sequence file')
    dcavisualizer.add_argument('dca_file', help = 'path to DCA file')
    dcavisualizer.add_argument('--contact_dist',
        type = float,
        help = 'Maximum distance below which residue pairs are considered contacts',
    )
    dcavisualizer.add_argument('--num_dca_contacts', type=int,
        help = 'Number of top DCA ranked pairs to be taken',
    )
    dcavisualizer.add_argument('--rna_secstruct_file',
        help = 'File containing RNA secondary structure information',
    )
    dcavisualizer.add_argument('--linear_dist', type=int,
        help = 'Distance between two residues in the sequence above which they'
            ' are considered to be tertiary contacts. For RNAs this parameter'
            ' sets the distance from secondary structure pairs. If (k, l) form'
            ' RNA secondary structure pairs, and the value of linear distance'
            ' is n, (2n + 1)^2 pairs are excluded from the tertiary contacts list.'
            ' That is, all pairs (-n + k, -n + l), (-n + 1 + k, -n + l), ...,'
            ' (k, l), ..., (k + n, l + n) are removed from DCA predicted contacts.'
    )
    dcavisualizer.add_argument('--verbose', action='store_true')
    tprate = subparser.add_parser(
        'tprate',
        help = 'subcommand for visualizing tprate',
    )
    #add pdb_content subcommand arguments
    pdb_parser.add_argument(
        'pdb_file',
        help = 'PDB file path',
    )
    pdb_parser.add_argument('--verbose', action='store_true')
    pdb_parser.add_argument(
        '--biomolecule',
        choices = ('protein', 'PROTEIN', 'rna', 'RNA'),
        help = 'Type of PDB structure i.e., protein or RNA'
    )
    pdb_parser.add_argument(
        '--struct_info', action='store_true',
        help = 'Displays a summary of the PDB file content.'
            ' This is done by extracting metadata from the PDB file header.',
        required = False
    )
    pdb_parser.add_argument(
        '--show_chains', action = 'store_true',
        help = 'Displays chain ID found in a PDB structure.',
        required = False,
    )
    pdb_parser.add_argument(
        '--pdb_id',
        help = 'ID of PDB structure',
        required = False,

    )
    #add refseq_content subcommands
    refseq_parser.add_argument(
        'refseq_file',
        help = 'FASTA formatted text file containing reference sequence(s)'
            'that has been used for DCA computation.',
    )
    refseq_parser.add_argument('--verbose', action='store_true')
    refseq_parser.add_argument(
        '--biomolecule',
        help = 'The biomolecule type the reference sequences represent.'
            'It should be either protein or RNA, case insensitive.',
        choices = ('protein', 'PROTEIN', 'rna', 'RNA'),
    )
    refseq_parser.add_argument(
        '--show_seqs',
        action = 'store_true',
        help = 'Displays the list of sequences read from Fasta formatted file',
    )

    secstruct_parser = subparser.add_parser('secstruct_content',
        help = 'RNA secondary structure parser',
    )
    secstruct_parser.add_argument('secstruct_file',
        help = 'File containing RNA secondary structure',
    )
    secstruct_parser.add_argument('--verbose', action='store_true')

    dca_parser = subparser.add_parser('dca_content',
        help = 'DCA results subparser',
    )
    dca_parser.add_argument('dca_file',
        help = 'File containing ranked site pairs by DCA score',
    )
    dca_parser.add_argument('--verbose', action = 'store_true')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    if args.verbose: configure_logging()

    if args.subcommand == 'pdb_content':
        pdb_test_inst = PDBContent(args.pdb_file, biomolecule = args.biomolecule)
        if args.struct_info: pdb_test_inst.show_struct_info()
        if args.show_chains: pdb_test_inst.display_chain_sequences()
    if args.subcommand == 'refseq_content':
        refseq_test_inst = RefSeqContent(args.refseq_file, biomolecule = args.biomolecule)
        if args.show_seqs:
            refseq_test_inst.display_reference_sequences()
    if args.subcommand == 'dcavisualizer':
        dcavisualizer = DCAVisualizer(
            args.biomolecule,
            args.pdb_chain_id,
            args.pdb_file,
            refseq_file = args.refseq_file,
            dca_file = args.dca_file,
            rna_secstruct_file = args.rna_secstruct_file,
            linear_dist = args.linear_dist,
            contact_dist = args.contact_dist,
            num_dca_contacts = args.num_dca_contacts,
        )
        #dcavisualizer.dca_predicted_tertiary_contacts()
        #dcavisualizer.map_pdbseq_to_refseq()
        #dcavisualizer.get_mapped_pdb_contacts()
        dcavisualizer.plot_contact_map()
        #dcavisualizer.compute_true_positive_rates()
        #dcavisualizer.plot_true_positive_rate()

    if args.subcommand == 'secstruct_content':
        secstruct_content = RNASecStructContent(args.secstruct_file)
    if args.subcommand == 'dca_content':
        dca_content = DCAContent(args.dca_file)
        pairs = dca_content.dca_ranked_pairs
