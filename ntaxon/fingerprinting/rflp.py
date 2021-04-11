import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from ntaxon.nucleotide import Sequence
from ntaxon.fingerprinting import BinaryMatrix


class RestrictionFragmentTable:
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

        :return: Filtered RestrictionFragmentTable
        """
        # TODO Validate Inputs
        pass

    def compare_fragment_size(self, query_table):
        print(0)


class RestrictionDigestion:
    __restriction_enzyme = None     # Biopython RestrictionEnzyme
    __species_names = None          # PD Series
    __accessions = None             # PD DataFrame
    __fragments = None              # PD DataFrame
    __binary_matrix = None          # PD DataFrame

    @property
    def binary_matrix(self):
        return self.__binary_matrix

    @property
    def fragments(self):
        return self.__fragments

    def __init__(self, enzyme, accessions: pd.DataFrame, label_col: str = 'sample', sequence_col: str = 'sequence'):
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
        self.__binary_matrix = self.__construct_binary_matrix()

    def __perform_restriction_digestion(self):
        """
        Performs Restriction Digestion on Nucleotide Sequences

        :param enzyme: Biopython Restriction Enzyme Class
        :param label_col: str - Name of Column for Sample Labels
        :param sequence_col: str - Name of Column for Sample DNA/RNA Sequence

        :return: RestrictionFragmentTableProfile - Digestion Profile Matrix
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

        return RestrictionFragmentTable(digestion)

    def __construct_binary_matrix(self):
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
