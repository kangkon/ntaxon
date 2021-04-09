import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
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
            # TODO: Check if index is duplicated
            isolate = accessions.loc[idx, label_col] # accessions[label_col][idx]
            seq = Sequence(accessions.loc[idx, sequence_col])
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

    def __construct_binary_matrix(self, min_size=1, min_frequency=1, min_common_haplotye=1):
        """
        Construct binary matrix from digestion profile

        :param min_size: Int - Minimum fragment size
        :param min_frequency: Int - Minimum frequency // Column wise (i.e. minimum times that a band occurs in samples/lanes/columns)
        :param direction: Int - Direction of frequency count 0 - column/samples 1 - row/haplotypes

        :return: BinaryMatrix
        """
        digestion = self.__digestion_table.data()
        isolates = digestion['isolate'].unique()

        if min_common_haplotye > 1:
            digestion = digestion[digestion['size'] > min_size].groupby('isolate').filter(
                lambda x: len(x) > min_common_haplotye)

        if min_size > 1 or min_frequency > 1:
            digestion = digestion[digestion['size'] > min_size].groupby('haplotype').filter(
                lambda x: len(x) > min_frequency)

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

        matrix_bin = BinaryMatrix(data=digestion_bin, locus_name=self.__restriction_enzyme.__name__, size_column="size")

        # Fill excluded isolates with value 0 in haplotypes filtered
        if matrix_bin.shape()[1] != len(isolates):
            # find columns not in current binary matrix
            #excluded_isolates = np.setdiff1d(matrix_bin.get_columns(), isolates)
            excluded_isolates = isolates[~np.isin(isolates, matrix_bin.get_columns())]
            matrix_bin.add_columns(columns=excluded_isolates, value=0)

        # if lane_order:
        #    if len(lane_order) == matrix_bin.shape[1]:
        #        matrix_bin = matrix_bin[lane_order]
        #    else:
        #        raise Exception("Length of lane_order is not same as binary matrix")

        return matrix_bin

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

    def plot_electrophoretic_diagram(self, figsize=(16, 10), min_size=1, min_frequency=1, min_common_haplotye=1,
                                     lane_order=None, size_tolerance=0):
        """
        Plot heatmap based on presence and absence of bands

        :param figsize: Tuple - Matplotlib plot size
        :param min_size: Int - Minimum fragment size
        :param min_frequency: Int - Minimum frequency
        :param min_common_haplotye: Int - Filter by minimum common haplotypes

        :param size_tolerance: Int, Group by certain fragment size tolerance
        """
        matrix_bin = self.__digestion_profile.get_binary_matrix()

        if min_size > 1 or min_frequency > 1 or min_common_haplotye > 1:
            matrix_bin = self.__construct_binary_matrix(min_size=min_size, min_frequency=min_frequency,
                                                        min_common_haplotye=min_common_haplotye).get_binary_matrix()
        if size_tolerance > 0:
            matrix_bin = self.apply_tolerance(size_tolerance=size_tolerance, matrix=matrix_bin)

        if lane_order:
            if len(lane_order) == matrix_bin.shape[1]:
                matrix_bin = matrix_bin[lane_order]
            else:
                raise Exception("Length of lane_order is not same as binary matrix")

        plt.figure(figsize=figsize)
        ax = sns.heatmap(matrix_bin, cmap="YlGnBu", cbar=False)
        plt.show()

    def apply_tolerance(self, size_tolerance=1, matrix=None):
        matrix_bin = matrix if matrix is not None else self.__digestion_profile.get_binary_matrix()
        fragments = []
        for f in matrix_bin.index.values.to_list():
            # get the size by spliting across _
            splitted = f.split('_')
            fragments.append([f, splitted[0]])

        diff_matrix = []
        for f, s in fragments:
            dif = []
            for _f, _s in fragments:
                dif.append(abs(float(s) - float(_s)))  # Get absolute values only
            diff_matrix.append(dif)

        diff_matrix = np.tril(diff_matrix, -1)  # convert to lower triangular form
        diff_matrix = pd.DataFrame(diff_matrix, columns=matrix_bin.index.values.to_list(),
                                   index=matrix_bin.index.values.to_list())

        diff_matrix_stacked = diff_matrix.stack().to_frame().reset_index()
        diff_matrix_stacked.columns = ['frag_1', 'frag_2', 'difference']

        # Filter by tolerance
        tolerance_filtered = diff_matrix_stacked[diff_matrix_stacked['difference'] <= size_tolerance]
        tolerance_filtered = tolerance_filtered[tolerance_filtered['difference'] != 0]

        # Compute homogenity for a fragment
        sample_count = matrix_bin.shape[1]  # no of samples
        fragment_occurance_count = matrix_bin.sum(axis=1)

        def compute_merge_vars(r):
            # get homogenity scores (fragment count scores)
            frag_1_hom = fragment_occurance_count.loc[r['frag_1']]
            frag_2_hom = fragment_occurance_count.loc[r['frag_2']]

            frag_sum = matrix_bin.loc[[r['frag_1'], r['frag_2']], :].sum(axis=0)
            return pd.Series([
                sample_count - frag_sum[frag_sum == 0].count(),
                r['frag_1'] if frag_1_hom > frag_2_hom else r['frag_2']
            ])

        tolerance_filtered[['merged_frag_count', 'merge_to']] = tolerance_filtered.apply(compute_merge_vars, axis=1)
        # Sort to merge fragments in correct order
        tolerance_filtered = tolerance_filtered.sort_values(by='merged_frag_count')
        # Construct the merged matrix
        matrix_bin_merged = matrix_bin.reset_index()
        for idx, r in tolerance_filtered.iterrows():
            mergable_frag = r['frag_1'] if r['merge_to'] == r['frag_2'] else r['frag_2']
            # mergable_frag_iloc = matrix_bin.index.get_loc(mergable_frag)
            # matrix_bin_merged.iloc[mergable_frag_iloc, [0]] = r['merge_to']
            matrix_bin_merged['haplotype'] = matrix_bin_merged['haplotype'].replace(mergable_frag, r['merge_to'])

        # Create new fresh dataframe # otherwise groupby causes problem
        df = pd.DataFrame(matrix_bin_merged.values.tolist(), columns=matrix_bin_merged.columns)
        # Combine rows and replace value > 1 to 1
        return df.groupby('haplotype').agg('sum').apply(lambda x: [1 if y > 0 else 0 for y in x])

    def get_haplotypes(self, size_tolerance=0, matrix=None, exclude=None, allow_polymorphism=0):
        matrix_bin = matrix if matrix is not None else self.__digestion_profile.get_binary_matrix()

        if exclude is not None:
            columns = list(matrix_bin.columns)

            if type(exclude) == str:
                if exclude not in columns:
                    raise Exception(f"{exclude} does not exist")
            elif type(exclude) == list:
                for c in exclude:
                    if c not in columns:
                        raise Exception(f"{c} does not exist")
            else:
                raise Exception("Invalid Exclude Value")

            matrix_bin.drop(exclude, axis=1)

        if size_tolerance > 0:
            matrix_bin = self.apply_tolerance(size_tolerance=size_tolerance, matrix=matrix_bin)

        sample_count = matrix_bin.shape[1]

        polymorphism_scores = 1 - matrix_bin.sum(axis=1) / sample_count
        polymorphism_scores = polymorphism_scores.round(decimals=1)
        return polymorphism_scores[polymorphism_scores <= allow_polymorphism].index.tolist()