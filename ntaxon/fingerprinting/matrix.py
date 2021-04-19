import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance


class BinaryMatrix():
    """
    BinaryMatrix od 0 and 1 mostly used for gel scoring

    Attributes
    ----------
    locus_name  : str
        Name of the locus or Name of the marker i.e. current primer or restriction enzyme

    Methods
    -------
    to_df() : DataFrame
        Binary Matrix of 0 for fragment absence and 1 for fragment presence

    """
    __binary_data_matrix_df = None
    __binary_matrix_df = None
    __boolean_matrix_df = None
    __binary_matrix_stacked_df = None

    __fragment_sizes = None
    __sample_count = 0

    locus_name = "Unknown"

    @property
    def shape(self):
        # returns: total_fragments, sample_count
        return self.__binary_matrix_df.shape

    @property
    def fragment_sizes(self):
        # Returns Pandas Indexed Series
        return self.__fragment_sizes

    @property
    def fragment_count(self):
        # Returns Pandas Indexed Series
        return self.__binary_data_matrix_df.shape[0]

    @property
    def sample_count(self):
        # Returns Int
        return self.__sample_count

    @property
    def sample_names(self):
        return list(self.__binary_matrix_df.columns)

    def __repr__(self):
        return f"<BinaryMatrix {self.fragment_count}x{self.sample_count} of Locus: {self.locus_name}>"

    def __str__(self):
        return f"<BinaryMatrix {self.fragment_count}x{self.sample_count} of Locus: {self.locus_name}>"

    def __init__(self, locus_name, data: pd.DataFrame, size_column='size'):
        # Convert NA to o
        data = data.fillna(0)

        #TODO: check data to be in proper format
        self.locus_name = locus_name

        if size_column in list(data.columns):
            data = data.rename(columns={size_column: "size"})
        else:
            data['size'] = [i for i in range(data.shape[0], 0, -1)]

        # Convert index values to string
        data.index.map(str)

        # Change Index and Columns name
        data.index.names = ["haplotype"]
        data.columns.names = ["sample"]

        self.__binary_data_matrix_df = data

        self.__apply_matrix_data()

    @staticmethod
    def __unstack_matrix_data(data):
        df = data.pivot_table(index=['haplotype', 'size'], columns='sample')
        df.columns = df.columns.droplevel().rename(None)
        df2 = df.reset_index().set_index('haplotype')
        df2.columns.name = "sample"
        return df2

    def __apply_matrix_data(self):
        data_bin = self.__binary_data_matrix_df.loc[:, self.__binary_data_matrix_df.columns != "size"]
        self.__binary_matrix_df = data_bin
        self.__boolean_matrix_df = data_bin.astype('bool')
        self.__binary_matrix_stacked_df = pd.melt(self.__binary_data_matrix_df.reset_index(), id_vars=['size', 'haplotype'])

        self.__fragment_sizes = self.__binary_data_matrix_df["size"]
        self.__sample_count = self.__binary_data_matrix_df.shape[1] - 1

    def to_df(self, size=False, transpose=False):
        out = self.__binary_matrix_df
        if size:
            out = self.__binary_data_matrix_df

        if transpose:
            return out.transpose()

        return out

    def to_list(self, size=False, transpose=False):
        data = self.__binary_data_matrix_df

        if transpose:
            data = data.transpose()

        if size:
            return data.values.tolist()

        return data.values.tolist()

    def to_df_bool(self):
        return self.__boolean_matrix_df

    def to_df_stacked(self, size=True):
        if size is False:
            return self.__binary_matrix_stacked_df.drop('size', axis=1)

        return self.__binary_matrix_stacked_df

    def set_fragment_sizes(self, sizes):
        # sizes shoud be array
        if len(sizes) != self.fragment_count:
            raise Exception("Size array does not match matrix dimension")
        data = self.__binary_data_matrix_df
        data['size'] = sizes
        self.__binary_data_matrix_df = data
        self.__apply_matrix_data()

    def set_sample_names(self, names):
        # names should be array
        if len(names) != self.sample_count:
            raise Exception("Names array does not match matrix dimension")
        data = self.__binary_data_matrix_df
        data.columns = ["size"] + [n for n in names]
        self.__binary_data_matrix_df = data
        self.__apply_matrix_data()

    def sample_order(self, names):
        if len(names) != self.__sample_count:
            raise Exception("Names array does not match sample count")
        data = self.__binary_data_matrix_df
        data = data[['size'] + names]
        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_reordered")

    def add_sample(self, columns, value=0):
        # TODO: Modify function accept value as vector
        data = self.__binary_data_matrix_df
        data.loc[:, columns] = value
        self.__apply_matrix_data()

    def exclude(self, samples=None, fragments=None):
        # samples and fragments should be list
        data = self.__binary_data_matrix_df

        if type(samples) == list:
            if any(e not in list(data.columns) for e in samples):
                raise Exception("Invalid Sample Names to Exclude")
            data = data.drop(labels=samples, axis=1)

        if type(fragments) == list:
            if any(e not in list(data.index) for e in fragments):
                raise Exception("Invalid Fragment Index Names to Exclude")
            data = data.drop(labels=fragments, axis=0)

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def merge(self, matrix):
        """
        Merge Two Binary Matrix and returns a combined a BinaryMatrix

        Parameters
        ----------
        matrix : BinaryMatrix
            BinaryMatrix Instance
        """
        matrix_a = self.__binary_data_matrix_df
        matrix_b = matrix.to_df(size=True)

        return BinaryMatrix(
            data=matrix_a.append(matrix_b, ignore_index=True),
            locus_name=matrix_a.locus_name + '_' + matrix_b.locus_name
        )

    # Get the fragment or haplotype size
    def sizeof(self, fragment: str = None, ifragment: int = None):
        """
        Get the size of a fragment
        :param fragment: str
        :param ifragment: int
        :return: int
        """
        if self.__fragment_sizes is None:
            return 0

        if ifragment:
            return self.__fragment_sizes.iloc[ifragment]
        else:
            return self.__fragment_sizes.loc[fragment]

    def sort_samples(self, distmethod="jaccard", linkage_method="ward"):
        """
        Sort sample order by similarity

        :return: Binarymatrix - sorted binary matrix
        """
        # TODO: Sort by similarity
        matrix = self.__binary_matrix_df
        print("Sorted matrix")
        return 0

    def filter_size(self, min=1, max=None):
        data = self.__binary_data_matrix_df
        if min > 1:
            data = data[data['size'] > min]
        if (max is not None) and (max > min):
            data = data[data['size'] < max]

        if data.empty:
            raise Exception("Result is empty")

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def filter_frequency(self, min=1, max=None, axis=1):
        # TODO complete frequency filter
        data = self.__binary_data_matrix_df
        if min > 1:
            data = data[data['size'] > min]
        if (max is not None) and (max > min):
            data = data[data['size'] < max]

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def filter_min_common_haplotype(self, value):
        data_stacked = self.__binary_matrix_stacked_df
        data_stacked = data_stacked[data_stacked['value'] > 0].groupby('sample').filter(lambda x: len(x) > value)
        if data_stacked.empty:
            raise Exception("Result is empty")

        data = self.__unstack_matrix_data(data=data_stacked)
        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def filter_haplotype_frequency(self, frequency=0, count=1):
        if count > self.sample_count:
            raise Exception("Haplotype Frequency Count is grater than sample count")
        df_stacked = self.__binary_matrix_stacked_df
        # filter all whose haplotype value is 1 and get count of all haplotype across samples
        haplotype_counts = df_stacked[df_stacked['value'] > 0].groupby(["haplotype"])["value"].count()

        # filter haplotypes by min common occurance
        filtered_haplotypes = haplotype_counts[haplotype_counts >= count]

        # get the sample names by filtered haplotypes
        filtered_samples = df_stacked[df_stacked['haplotype'].isin(list(filtered_haplotypes.index))]

        data = self.__unstack_matrix_data(filtered_samples)
        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def filter_by_fragment(self, fragments=[], sizes=[]):
        """
        Filter samples based on fragment presence

        Parameters
        ----------
        fragments : list
            List of haplotype names (or fragment names)
        size : list
            List of haplotype sizes (or fragment sizes)
        """
        if type(fragments) != list or type(sizes) != list:
            raise Exception("Invalid input")

        # only values with 1
        data_stacked = self.__binary_matrix_stacked_df[self.__binary_matrix_stacked_df['value'] == 1]

        results = None
        q_list = []

        if len(fragments) > 0:
            results = data_stacked[data_stacked['haplotype'].isin(fragments)]
            q_list = fragments

        elif len(sizes) > 0:
            results = data_stacked[data_stacked['size'].isin(sizes)]
            q_list = sizes

        if results is not None:
            results = results.groupby('sample').filter(lambda x: len(x) == len(q_list))

        if results is None or results.empty:
            raise Exception("No results found")

        data = self.__unstack_matrix_data(results)
        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def select_samples(self, names=None):
        # names should be array
        data = self.__binary_data_matrix_df
        if type(names) == list:
            if any(e not in list(data.columns) for e in names):
                raise Exception("Invalid Sample Names to Exclude")
            data = data[['size'] + names]

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def select_fragments(self, names=None):
        # fragment index names should be array
        data = self.__binary_data_matrix_df
        if type(names) == list:
            if any(e not in list(data.index) for e in names):
                raise Exception("Invalid Fragment Index Names to Exclude")
            data = data.loc[names, :]

        if data.empty:
            raise Exception("Result is empty")

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def select_ifragments(self, indices=None):
        # fragment index names should be array
        data = self.__binary_data_matrix_df
        if type(indices) == list:
            if any(e > self.fragment_count or e <= 1 for e in indices):
                raise Exception("Invalid Fragment Index Value")
            data = data.iloc[indices, :]

        return BinaryMatrix(data=data, locus_name=f"{self.locus_name}_filtered")

    def total_fragment_count(self, presence=True, samplewise=False):
        """
        Counts fragments in a binary matrix. Returns a series with index equivalent to binary matrix

        Parameters
        ----------
        presence : bool
            Allele present or absent i.e allele 0 or allele 1
        samplewise : bool, optional
            For returning sample wise fragment Count
        """
        matrix = self.__binary_matrix_df

        if samplewise:
            matrix_t = matrix.transpose()
            if presence:
                return matrix_t.sum(axis=1)
            else:
                return (matrix_t == 0).astype(int).sum(axis=1)

        if presence:
            return matrix[matrix.eq(1)].count(axis=1)
        else:
            return matrix[matrix.eq(0)].count(axis=1)

    def total_fragments(self, unique=True):
        """
        Counts total fragments present in a binary matrix.
        Parameters
        ----------
        unique : bool
            Sum of Unique fragments or Repeating fragments
        """
        if not unique:
            return self.total_fragment_count().sum()
        return self.__binary_matrix_df.shape[0]

    def fragment_frequency(self, singular=True, presence=True):
        """
        Calculate Allele frequency

        Parameters
        ----------
        singular : bool
            True: Calculate frequency of individual fragment based on total occurrence of that particular fragment
            False: Calculate frequency of individual fragment based on total occurrence of all fragments in matrix
        """

        if singular:
            return self.total_fragment_count(presence=presence) / self.__sample_count

        return self.total_fragment_count(presence=True)/self.total_fragments(unique=False)

    def apply_tolerance(self, size_tolerance=1):
        ## Aply size tolerance and returns a dataframe
        #matrix_bin = matrix if matrix is not None else self.__binary_matrix_df
        fragments = list(self.fragment_sizes.reset_index().values)
        fragment_index_names = list(self.__binary_data_matrix_df.index.values)

        diff_matrix = []
        for f, s in fragments:
            dif = []
            for _f, _s in fragments:
                dif.append(abs(float(s) - float(_s)))  # Get absolute values only
            diff_matrix.append(dif)

        diff_matrix = np.tril(diff_matrix, -1)  # convert to lower triangular form

        diff_matrix = pd.DataFrame(diff_matrix, columns=fragment_index_names,
                                   index=fragment_index_names)

        diff_matrix_stacked = diff_matrix.stack().to_frame().reset_index()
        diff_matrix_stacked.columns = ['frag_1', 'frag_2', 'difference']

        # Filter by tolerance
        tolerance_filtered = diff_matrix_stacked[diff_matrix_stacked['difference'] <= size_tolerance]
        tolerance_filtered = tolerance_filtered[tolerance_filtered['difference'] != 0]

        # Compute homogenity for a fragment
        sample_count = self.__sample_count  # no of samples
        each_fragment_count = self.__binary_matrix_df.sum(axis=1)
        fragment_sizes = self.__fragment_sizes

        matrix_bin = self.__binary_matrix_df

        def compute_merge_vars(r):
            # get homogenity scores (fragment count scores)
            frag_1_hom = each_fragment_count.loc[r['frag_1']]
            frag_2_hom = each_fragment_count.loc[r['frag_2']]

            frag_sum = matrix_bin.loc[[r['frag_1'], r['frag_2']], :].sum(axis=0)
            merge_to_frag = r['frag_1'] if frag_1_hom > frag_2_hom else r['frag_2']
            return pd.Series([
                sample_count - frag_sum[frag_sum == 0].count(),
                merge_to_frag,
                fragment_sizes[merge_to_frag],
            ])

        tolerance_filtered[['merged_frag_count', 'merge_to', 'merge_size']] = tolerance_filtered.apply(compute_merge_vars, axis=1)
        # Sort to merge fragments in correct order
        tolerance_filtered = tolerance_filtered.sort_values(by='merged_frag_count')

        # Construct the merged matrix
        matrix_bin_merged = self.__binary_data_matrix_df.reset_index()
        for idx, r in tolerance_filtered.iterrows():
            mergable_frag = r['frag_1'] if r['merge_to'] == r['frag_2'] else r['frag_2']
            # mergable_frag_iloc = matrix_bin.index.get_loc(mergable_frag)
            # matrix_bin_merged.iloc[mergable_frag_iloc, [0]] = r['merge_to']
            matrix_bin_merged['haplotype'] = matrix_bin_merged['haplotype'].replace(mergable_frag, r['merge_to'])
            matrix_bin_merged.loc[matrix_bin_merged['haplotype'] == r['merge_to'], 'size'] = fragment_sizes.loc[r['merge_to']]


        # Create new fresh dataframe # otherwise groupby causes problem
        df = pd.DataFrame(matrix_bin_merged.values.tolist(), columns=matrix_bin_merged.columns)
        # Combine rows and replace value > 1 to 1
        df_fixed = df.groupby(['haplotype', 'size']).agg('sum')\
            .apply(lambda x: [1 if y > 0 else 0 for y in x])\
            .reset_index()\
            .set_index('haplotype')

        return BinaryMatrix(data=df_fixed, locus_name=f"{self.locus_name}_Tolarated")
    
    @property
    def monomorphic_fragments(self):
        # Row wise Total fragment count == sample count
        return self.__binary_matrix_df.loc[self.total_fragment_count(presence=True) == self.shape[1]]

    @property
    def monomorphic_fragment_count(self):
        # Row wise Total fragment count == sample count
        return len(self.__binary_matrix_df.loc[self.total_fragment_count(presence=True) == self.shape[1]].index)

    @property
    def polymorphic_fragments(self):
        return self.__binary_matrix_df.loc[self.total_fragment_count(presence=True) != self.shape[1]]

    @property
    def polymorphic_fragment_count(self):
        # Row wise Total fragment count == sample count
        return len(self.__binary_matrix_df.loc[self.total_fragment_count(presence=True) != self.shape[1]].index)

    @property
    def percent_polymorphism(self):
        return self.polymorphic_fragment_count * 100/self.total_fragments()

    @property
    def marker_frequency(self):
        """
        Calculates Marker Frequency

        :return: Marker Frequency
        """
        return self.total_fragments(unique=False) / (self.total_fragments(unique=True) * self.__sample_count)

    def heterozygosity(self, singular=False):
        """
        Calculate Heterozygosity for current marker depending on frequency of homozygous individuals
        H = 1 - SUM(p_i^2), where p = frequency of i th allele (Nei, 1978)

        :param singular: (Bool) True - Calculate frequency of individual fragment based on total occurrence of that particular fragment;
                        False -  Calculate frequency of individual fragment based on total occurrence of all fragments in matrix

        :return: PIC value
        """
        allele_frequency = self.fragment_frequency(singular=singular)
        sigma_pi = 0
        for fragment in list(self.__binary_matrix_df.index):
            sigma_pi = sigma_pi + np.square(allele_frequency.loc[fragment])

        return 1 - sigma_pi

    def heterozygosity_2(self, singular=False):
        """
        Calculate Heterozygosity for current marker depending on frequency of heterozygous individuals
        H = 1 - PRODUCT(2 * p_i), where p = frequency of i th allele

        :param singular: (Bool) True - Calculate frequency of individual fragment based on total occurrence of that particular fragment;
                        False -  Calculate frequency of individual fragment based on total occurrence of all fragments in matrix

        :return: PIC value
        """
        allele_frequency = self.fragment_frequency(singular=singular)
        product_pi = 1
        for fragment in list(self.__binary_matrix_df.index):
            product_pi = product_pi * 2 * allele_frequency.loc[fragment]

        return product_pi

    def pic_dominant_1(self):
        """
        Calculate Polymorphic Information Content for a Dominant marker (Roldan-Ruiz et al., 2000)
        PIC = 2f(1-f)

        :return: A list of PIC for individual fragments/allele or average PIC value
        """
        marker_frequency = self.marker_frequency
        return 2 * marker_frequency * (1 - marker_frequency)

    def pic_dominant_2(self):
        """
        Calculate Polymorphic Information Content for a Dominant marker (De Riek et al., 2001)
        PIC = 1 - [f^2 + (1 - f)^2]

        :return: A list of PIC for individual fragments/allele or average PIC value
        """
        marker_frequency = self.marker_frequency
        return 1 - (np.square(marker_frequency) + np.square((1 - marker_frequency)))

    def pic_dominant_3(self, mean=True):
        """
        Calculate Polymorphic Information Content for a Dominant marker
        PIC = 1 - (p^2 + q^2), where p = frequency of presence and q = frequency of absence in a particular allele

        :param mean: (boolean) Output The average PIC

        :return: A list of PIC for individual fragments/allele or average PIC value
        """
        p = self.fragment_frequency(presence=True)
        q = self.fragment_frequency(presence=False)

        pic = 1 - (np.square(p) + np.square(q))

        if mean:
            return pic.mean()

        return pic

    def pic_codominant_1(self, singular=False):
        """
        Calculate Polymorphic Information Content for a Dominant marker
        PIC = 1 - SUM(p_i^2), where p = frequency of i th allele

        :param singular: (Bool) True - Calculate frequency of individual fragment based on total occurrence of that particular fragment;
                        False -  Calculate frequency of individual fragment based on total occurrence of all fragments in matrix

        :return: PIC value
        """
        allele_frequency = self.fragment_frequency(singular=singular)
        sigma_pi = 0
        for fragment in list(self.__binary_matrix_df.index):
            sigma_pi = sigma_pi + np.square(allele_frequency.loc[fragment])

        return 1 - sigma_pi

    def pic_codominant_2(self, singular=False):
        """
        Calculate Polymorphic Information Content for a Codominani marker
        PIC = 1 - SUM(P_i^2) - (SUM(P_i^2])^2 + SUM(P_i^4), where p = frequency of i th allele

        :param singular: (Bool) True - Calculate frequency of individual fragment based on total occurrence of that particular fragment;
                        False -  Calculate frequency of individual fragment based on total occurrence of all fragments in matrix

        :return: PIC value
        """
        allele_frequency = self.fragment_frequency(singular=singular)
        sigma_pi_sq = 0
        sigma_pi_4 = 0
        for fragment in list(self.__binary_matrix_df.index):
            sigma_pi_sq = sigma_pi_sq + np.square(allele_frequency.loc[fragment])
            sigma_pi_4 = sigma_pi_4 + np.power(allele_frequency.loc[fragment], 4)

        return 1 - sigma_pi_sq - np.square(sigma_pi_sq) + sigma_pi_4

    def pic_codominant_3(self, singular=False):
        # TODO - Fix the BUG
        """
        Calculate Polymorphic Information Content for a Codominant marker for population i and j (Botstein, 1980)

        :param singular: (Bool) True - Calculate frequency of individual fragment based on total occurrence of that particular fragment;
                        False -  Calculate frequency of individual fragment based on total occurrence of all fragments in matrix

        :return: PIC value
        """
        alleles = list(self.__binary_matrix_df.index)
        allele_frequency = self.fragment_frequency(singular=singular)
        sigma_pi_sq = 0

        # fragment here represents an allele
        for allele in alleles:
            sigma_pi_sq = sigma_pi_sq + np.square(allele_frequency.loc[allele])

            sigma_pi_pj = 0
            for i in range(len(allele) - 1):
                for j in range(i + 1, len(allele)):
                    sigma_pi_pj = sigma_pi_pj + 2 * np.square(allele_frequency.iloc[i]) * np.square(allele_frequency.iloc[j])

        return 1 - sigma_pi_sq + sigma_pi_pj
    
    def get_haplotypes(self, size_tolerance=0, exclude=None, allow_polymorphism=0):
        matrix_bin = self

        if exclude is not None:
            matrix_bin = matrix_bin.exclude(samples=exclude)

        if size_tolerance > 0:
            matrix_bin = matrix_bin.apply_tolerance(size_tolerance=size_tolerance)

        sample_count = matrix_bin.sample_count

        polymorphism_scores = 1 - matrix_bin.to_df().sum(axis=1) / sample_count
        polymorphism_scores = polymorphism_scores.round(decimals=1)
        return polymorphism_scores[polymorphism_scores <= allow_polymorphism].index.tolist()
    
    def plot_histogram(self, figsize=(14, 8)):
        """
        Plot Histogram of fragment size based occurance of fragments

        :param figsize: Tuple - Matplotlib plot size
        """
        plt.figure(figsize=figsize)
        ax = sns.histplot(data=self.__binary_matrix_stacked_df, x="size")
        plt.show()

    def plot_electrophoretic_diagram(self, figsize=(16, 10), min_size=1, min_frequency=1, min_common_haplotye=1,
                                     lane_order=None, size_tolerance=0, exclude=None):
        """
        Plot heatmap based on presence and absence of fragments

        :param figsize: Tuple - Matplotlib plot size
        :param min_size: Int - Minimum fragment size
        :param min_frequency: Int - Minimum frequency
        :param min_common_haplotye: Int - Filter by minimum common haplotypes
        :param exclude: list - Exclude samples on list

        :param size_tolerance: Int, Group by certain fragment size tolerance
        """
        matrix_bin = self

        if exclude is not None:
            matrix_bin = matrix_bin.exclude(samples=exclude)
        
        if lane_order:
            matrix_bin = matrix_bin.sample_order(lane_order)

        if min_size > 1:
            matrix_bin = matrix_bin.filter_size(min=min_size)

        if min_common_haplotye > 1:
            matrix_bin = matrix_bin.filter_min_common_haplotype(min_common_haplotye)

        #or min_frequency > 1 or min_common_haplotye > 1:
        #    matrix_bin = self.__construct_binary_matrix(min_size=min_size, min_frequency=min_frequency,
        #                                                min_common_haplotye=min_common_haplotye)
        
        if size_tolerance > 0:
            matrix_bin = matrix_bin.apply_tolerance(size_tolerance=size_tolerance)

        matrix_bin_data = matrix_bin.to_df()

        plt.figure(figsize=figsize)
        ax = sns.heatmap(matrix_bin_data, cmap="YlGnBu", cbar=False)
        plt.show()

    def distance(self, method='jaccard', squareform=False):
        binary_data = self.__binary_matrix_df.transpose() # Required rows as Samples, columns as Haplotypes
        pdist_vector = 1 - distance.pdist(binary_data, method)

        if squareform:
            return distance.squareform(pdist_vector)

        return pdist_vector