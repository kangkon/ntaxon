import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ntaxon.nucleotide import Sequence
from ntaxon.fingerprinting import BinaryMatrix


class SeqDigestionTable:
    __digestion_profile = pd.DataFrame()

    def __init__(self, digestion_table: pd.DataFrame):
        self.__digestion_table = digestion_table

    def data(self):
        return self.__digestion_table

    def screen_haplotypes(self, sizeonly=True, size_min = 50, similarity_max= 0.8):
        """
        Screen for Unique Haplotypes for a species
        :param sizeonly: Boolean - True: Screen based on Nucleotide fragment size only, False - Screen based on Nucleotide Sequence
        :param size_min: Int - Minimum fragment size to consider
        :param similarity_max: Maximum similarity consideration for seq based comparision

        :return: Filtered SeqDigestionTable
        """
        # TODO Validate Inputs

        if 'species' not in self.__digestion_table.columns:
            raise Exception("Species column does not exist in data frame")


    def compare_fragment_size(self, query_table):
        print(0)


class RestrictionDigestion:
    __restriction_enzyme = None
    __species_names = None
    __accession_input = None
    __digestion_table = None
    __digestion_profile = None

    def __init__(self, enzyme, accessions: pd.DataFrame, label_col: str = 'isolate', sequence_col: str = 'sequence',
                 species_names: pd.Series = None):
        """
        :param enzyme: Biopython Restriction Enzyme Class
        :param accessions: Accession having sequences
        :param label_col: str - Name of Column for Sample Labels
        :param sequence_col: str - Name of Column for Sample DNA/RNA Sequence
        :param species_names: pd.Series - Indexed pandas series containing species name, where index being isolate/sample name
        """
        self.__accession_input = accessions
        self.__restriction_enzyme = enzyme

        if 'species' in accessions.columns:
            self.__species_names = accessions.set_index(label_col)['species']

        # If species names are manually provided the existing values will be overwritten
        if species_names is not None:
            self.__species_names = species_names

        self.__digestion_table = self.__perform_restriction_digestion(
            enzyme=enzyme,
            label_col=label_col,
            sequence_col=sequence_col,
        )

        self.__digestion_profile = self.__construct_binary_matrix()

    def __perform_restriction_digestion(self, enzyme, label_col: str = 'isolate', sequence_col: str = 'sequence'):
        """
        Performs Restriction Digestion on Nucleotide Sequences

        :param enzyme: Biopython Restriction Enzyme Class
        :param label_col: str - Name of Column for Sample Labels
        :param sequence_col: str - Name of Column for Sample DNA/RNA Sequence

        :return: SeqDigestionProfile - Digestion Profile Matrix
        """
        accessions = self.__accession_input
        digestion = pd.DataFrame(columns=['size', 'isolate', 'fragment'])
        for idx in accessions.index:
            isolate = accessions[label_col][idx]
            seq = Sequence(accessions[sequence_col][idx])
            haplotypes = seq.haplotypes(enzyme)
            for frg in haplotypes:
                digestion = digestion.append(
                    pd.DataFrame([{"size": frg[0], "isolate": isolate, "fragment": str(frg[1])}]),
                    ignore_index=True
                )

        # Label Haplotypes Groupwise
        digestion['haplotype'] = digestion.groupby(['size', 'isolate']).cumcount() + 1  # Assign serial number
        digestion['haplotype'] = (digestion['size'].astype(str) + '_' + digestion['haplotype'].astype(str)).astype('category')

        # Label sample/isolate species
        species_names = self.__species_names
        if species_names is not None:
            digestion['species'] = digestion.apply(lambda x: species_names.loc[x.isolate], axis = 1).astype('category')

        return SeqDigestionTable(digestion)

    def __construct_binary_matrix(self, min_size=1, min_frequency=1):
        """
        Construct binary matrix from digestion profile

        :param min_size: Int - Minimum fragment size
        :param min_frequency: Int - Minimum frequency

        :return: BinaryMatrix
        """
        digestion = self.__digestion_table.data()

        digestion = digestion[digestion['size'] > min_size].groupby('haplotype').filter(lambda x: len(x) > min_frequency)

        # Compute the digestion profile [index of dataframe are sizes]
        digestion_bin = digestion[['size', 'isolate', 'haplotype']].pivot(index='haplotype', columns='isolate',
                                                                          values='size')
        # Change NaN to 0
        digestion_bin = digestion_bin.fillna(0)

        # Change other values to 1
        digestion_bin = digestion_bin.applymap(lambda x: 1 if x > 0 else 0)

        # Sort by size
        digestion_bin['size'] = [int(i.split('_')[0]) for i in digestion_bin.index]

        digestion_bin = digestion_bin.sort_values(by=["size"], axis=0, ascending=False)

        return BinaryMatrix(data=digestion_bin, locus_name=self.__restriction_enzyme.__name__, size_column="size")

    def digestion_table(self) -> SeqDigestionTable:
        return self.__digestion_table

    def digestion_profile(self) -> BinaryMatrix:
        return self.__digestion_profile

    def plot_histogram(self, figsize=(14, 8)):
        """
        Plot Histogram of band size based occurance of bands

        :param figsize: Tuple - Matplotlib plot size
        """
        plt.figure(figsize=figsize)
        ax = sns.histplot(data=self.__digestion_table.digestion_table(), x="size")
        plt.show()

    def plot_electrophoretic_diagram(self, figsize=(16, 10), min_size=1, min_frequency=1):
        """
        Plot heatmap based on presence and absence of bands

        :param figsize: Tuple - Matplotlib plot size
        :param min_size: Int - Minimum fragment size
        :param min_frequency: Int - Minimum frequency

        :param min_occurrence: Int, Filter by Minimum Occurrence
        """
        matrix_bin = self.__digestion_profile.get_binary_matrix()
        if min_size > 1 or min_frequency > 1:
            matrix_bin = self.__construct_binary_matrix(min_size=min_size, min_frequency=min_frequency).get_binary_matrix()

        plt.figure(figsize=figsize)
        ax = sns.heatmap(matrix_bin, cmap="YlGnBu", cbar=False)
        plt.show()

    def screen_ref_haplotype(self):
        print('Alignment')

