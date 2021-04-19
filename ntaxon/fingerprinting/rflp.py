import pandas as pd
import numpy as np
from Bio import SeqIO
from ntaxon.nucleotide import Sequence
from ntaxon.fingerprinting import BinaryMatrix


class _RestrictionFragmentMap:
    __data_df = pd.DataFrame()
    
    def __init__(self, map: pd.DataFrame):
        self.__data_df = map
    
    def to_df(self):
        return self.__data_df
    
    

class _RestrictionFragmentTable:
    __data_df = pd.DataFrame()

    def __init__(self, fragments: pd.DataFrame):
        self.__data_df = fragments

    def to_df(self):
        return self.__data_df

    def screen_haplotypes(self, sizeonly=True, size_min = 50, similarity_max= 0.8):
        """
        Screen for Unique Haplotypes for a species
        :param sizeonly: Boolean - True: Screen based on Nucleotide fragment size only, False - Screen based on Nucleotide Sequence
        :param size_min: Int - Minimum fragment size to consider
        :param similarity_max: Maximum similarity consideration for seq based comparision

        :return: Filtered _RestrictionFragmentTable
        """
        # TODO Validate Inputs
        pass

    def compare_fragment_size(self, query_table):
        print(0)


class RestrictionDigestion:
    __restriction_enzyme = None         # Biopython RestrictionEnzyme
    __species_names = None              # PD Series
    __accessions = None                 # PD DataFrame
    __fragments = None                  # PD DataFrame
    __map = None                        # PD DataFrame of restriction site locations
    __binary_matrix_size = None         # PD DataFrame
    __binary_matrix_map = None          # PD DataFrame

    # TODO: Develop method for binary matrix based on fragment similarity
    __binary_matrix_fragments = None    # PD DataFrame of binary matrix of similar sequences

    def binary_matrix(self, of='map'):
        return self.__binary_matrix_map

    @property
    def fragments(self):
        return self.__fragments

    def __init__(self, enzyme, fasta, descriptors=None):
        if fasta is None:
            raise Exception('Invalid Input')

        seq_list = []
        for record in SeqIO.parse(fasta, "fasta"):
            seq_list.append([record.id, str(record.seq)])

        self.__accessions = pd.DataFrame(data=seq_list, columns=['sample', 'sequence'])
        self.__restriction_enzyme = enzyme

        self.__fragments = self.__perform_restriction_digestion()
        self.__map = self.__perform_restriction_search()

        self.__binary_matrix_size = self.__construct_fragment_binary_matrix()
        self.__binary_matrix_map = self.__construct_map_binary_matrix()

    def __init_from_df(self, enzyme, accessions: pd.DataFrame, label_col: str = 'sample', sequence_col: str = 'sequence'):
        """
        :param enzyme: Biopython Restriction Enzyme Class
        :param accessions: Accession having sequences
        :param label_col: str - Name of Column for Sample Labels
        :param sequence_col: str - Name of Column for Sample Nucleotide Sequence
        """
        # Rename columns
        accessions = accessions.rename(columns={label_col: "sample", sequence_col: "sequence"})

        accessions = accessions.drop_duplicates(subset=["sample"], keep="last")
        accessions[sequence_col] = accessions["sequence"].str.strip()

        if "species" in list(accessions.columns):
            accessions['species'] = accessions['species'].str.strip()
            self.__species_names = accessions.set_index("sample")['species']

        self.__accessions = accessions
        self.__restriction_enzyme = enzyme

        self.__fragments = self.__perform_restriction_digestion()
        self.__binary_matrix_size = self.__construct_fragment_binary_matrix()
        self.__binary_matrix_map = self.__construct_map_binary_matrix()

    def __perform_restriction_search(self):
        """
        Performs Restriction Location Search on Nucleotide Sequences
        """
        accessions = self.__accessions
        r_map_df = pd.DataFrame(columns=['sample', 'location'])
        for i, r in accessions.iterrows():
            s = Sequence(r['sequence'])
            r_maps = s.restriction_search(self.__restriction_enzyme)
            for m in r_maps:
                r_map_df = r_map_df.append({
                    'sample': r['sample'],
                    'location': m
                }, ignore_index=True)
        return _RestrictionFragmentMap(map=r_map_df)
        

    def __perform_restriction_digestion(self):
        """
        Performs Restriction Digestion on Nucleotide Sequences

        :param enzyme: Biopython Restriction Enzyme Class
        :param label_col: str - Name of Column for Sample Labels
        :param sequence_col: str - Name of Column for Sample DNA/RNA Sequence

        :return: _RestrictionFragmentTableProfile - Digestion Profile Matrix
        """
        accessions = self.__accessions
        digestion = pd.DataFrame(columns=['size', 'sample', 'fragment'])
        for idx in accessions.index:
            # TODO: Check if index is duplicated
            sample = accessions.loc[idx, "sample"] # accessions[label_col][idx]
            seq = Sequence(accessions.loc[idx, "sequence"])
            haplotypes = seq.restriction_digest(self.__restriction_enzyme)
            for frg in haplotypes:
                digestion = digestion.append(
                    pd.DataFrame([{"size": frg[0], "sample": sample, "fragment": str(frg[1])}]),
                    ignore_index=True
                )

        # Label Haplotypes Groupwise
        digestion['haplotype'] = digestion.groupby(['size', 'sample']).cumcount() + 1  # Assign serial number
        digestion['haplotype'] = (digestion['size'].astype(str) + '_' + digestion['haplotype'].astype(str)).astype('category')

        # Label sample/sample species
        species_names = self.__species_names

        def get_species_name(x):
            return species_names.loc[x['sample']]

        if species_names is not None:
            digestion['species'] = digestion.apply(get_species_name, axis = 1).astype('category')

        return _RestrictionFragmentTable(digestion)

    def __construct_fragment_binary_matrix(self):
        """
        Construct binary matrix from digestion profile
        :return: BinaryMatrix
        """
        digestion = self.__fragments.to_df()
        samples = digestion['sample'].unique()

        # Compute the digestion profile [index of dataframe are sizes]
        digestion_bin = digestion[['size', 'sample', 'haplotype']].pivot(index='haplotype', columns='sample',
                                                                          values='size')
        # Change NaN to 0
        digestion_bin = digestion_bin.fillna(0)

        # Change other values to 1
        digestion_bin = digestion_bin.applymap(lambda x: 1 if x > 0 else 0)

        # Sort by size
        digestion_bin['size'] = [int(i.split('_')[0]) for i in digestion_bin.index]

        digestion_bin = digestion_bin.sort_values(by=["size"], axis=0, ascending=False)

        matrix_bin = BinaryMatrix(data=digestion_bin, locus_name=self.__restriction_enzyme.__name__)

        # Fill excluded samples with value 0 in haplotypes filtered
        if matrix_bin.sample_count != len(samples):
            excluded_samples = samples[~np.isin(samples, matrix_bin.sample_names)]
            matrix_bin.add_sample(columns=excluded_samples, value=0)

        return matrix_bin
    
    def __construct_map_binary_matrix(self):
        map = self.__map.to_df()

        # Convert to binary matrix
        bin_df = pd.crosstab(map['location'], map['sample'])

        # Mimic restriction location as size column
        bin_df['size'] = list(bin_df.index)
        # return BinaryMatrix(data=bin_df.reset_index(drop=True), locus_name=self.__restriction_enzyme.__name__)
        return BinaryMatrix(data=bin_df, locus_name=self.__restriction_enzyme.__name__)
