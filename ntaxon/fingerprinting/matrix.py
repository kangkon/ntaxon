import pandas as pd
import numpy as np

class BinaryMatrix(pd.DataFrame):
    """
    BinaryMatrix od 0 and 1 mostly used for gel scoring

    Attributes
    ----------
    locus_name  : str
        Name of the locus or Name of the marker i.e. current primer or restriction enzyme

    Methods
    -------
    get___binary_matrix() : DataFrame
        Binary Matrix of 0 for band absence and 1 for band presence

    """
    __locus_name = "Unknown"
    __sample_count = 0
    __raw_data = None
    __size_matrix = None
    __binary_matrix = None
    __boolean_matrix = None

    def __init__(self, locus_name, data: pd.DataFrame, size_column=None):
        # Convert NA to o
        data = data.fillna(0)

        #TODO: check data to be in proper format
        self.__raw_data = data
        self.__locus_name = locus_name

        if size_column:
            __binary_matrix = data.loc[:, data.columns != size_column]
            self.__binary_matrix = __binary_matrix
            self.__boolean_matrix = __binary_matrix.astype('bool')
            self.__size_matrix = data[size_column]
        else:
            self.__binary_matrix = data
            self.__boolean_matrix = data.astype('bool')

        self.__sample_count = self.__binary_matrix.shape[1]

    def get_binary_matrix(self):
        return self.__binary_matrix

    def get_boolean_matrix(self):
        return self.__boolean_matrix
    
    def get_size_matrix(self):
        return self.__size_matrix

    def get_locus_name(self):
        return self.__locus_name

    def get_sample_count(self):
        return self.__sample_count

    def set_band_sizes(self, sizes, match_index_label=True):
        # TODO Validate sizes matrix input
        self.__size_matrix = sizes

    def join(self, matrix):
        """
        Join Two Binary Matrix and returns a combined a BinaryMatrix

        Parameters
        ----------
        matrix : BinaryMatrix
            BinaryMatrix Instance
        """
        matrix_a = self.__binary_matrix
        matrix_b = matrix.get_binary_matrix()

        if self.__size_matrix is not None:
            matrix_a['size'] = self.__size_matrix

        if matrix.get_size_matrix() is not None:
            matrix_b['size'] = matrix.get_size_matrix()

        matrix_a = matrix_a.set_index(f'{self.__locus_name}_' + matrix_a.index.astype(str))
        matrix_b = matrix_b.set_index(f'{matrix.get_locus_name()}_' + matrix_b.index.astype(str))

        if matrix.get_size_matrix() is not None or self.__size_matrix is not None:
            return BinaryMatrix(locus_name=self.__locus_name + '_' + matrix.get_locus_name(),
                                data=matrix_a.append(matrix_b), size_column="size")
        else:
            return BinaryMatrix(locus_name=self.__locus_name + '_' + matrix.get_locus_name(),
                                data=matrix_a.append(matrix_b))

    # Get the band or haplotype size
    def sizeof(self, index: int = None, index_label: str = None):
        """
        Get the size of a band
        :param index: int
        :param index_label: str
        :return: int
        """
        if self.__size_matrix is None:
            return 0

        if index:
            return self.__size_matrix.iloc[index]
        else:
            return self.__size_matrix.loc[index_label]

    def sample_names(self):
        """
        Get list of sample names

        :return: list
        """
        return self.__binary_matrix.columns

    # returns: total_alleles, sample_count
    def shape(self):
        return self.__binary_matrix.shape

    def band_count(self, presence=True, samples=None, samplewise=False):
        """
        Counts bands in a binary matrix. Returns a series with index equivalent to binary matrix

        Parameters
        ----------
        presence : bool
            Allele present or absent i.e allele 0 or allele 1
        samples : ndarray, optional
            The array of sample or column names.
        samplewise : bool, optional
            Sample wise Band Count
        """
        matrix = self.__binary_matrix
        if samples:
            matrix = self.__binary_matrix[samples]

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

    def total_bands(self, unique=True, samples = None):
        """
        Counts total bands present in a binary matrix.
        Parameters
        ----------
        unique : bool
            Sum of Unique bands or Repeating bands
        samples : ndarray
            The array of sample or column names.
        """
        if not unique:
            return self.band_count(samples=samples).sum()
        return self.__binary_matrix.shape[0]

    def band_frequency(self, singular=True, samples=None):
        """
        Calculate Allele frequency

        Parameters
        ----------
        singular : bool
            True: Calculate frequency of individual band based on total occurrence for a particular bands
            False: Calculate frequency of individual band based on total occurrence of all bands
        samples : ndarray
            The array of sample or column names.
        """

        # TODO; Evaluate the frequency for singular param

        #if singular:
        #    return


        return self.band_count(presence=True, samples=samples)/self.total_bands(unique=False, samples=samples)

    def monomorphic_bands(self, samples=None):
        # Row wise Total band count == sample count
        if samples:
            return self.__binary_matrix.loc[self.band_count(presence=True, samples=samples) == len(samples)][samples]
        return self.__binary_matrix.loc[self.band_count(presence=True) == self.shape()[1]]

    def monomorphic_band_count(self, samples=None):
        # Row wise Total band count == sample count
        if samples:
            return len(self.__binary_matrix.loc[self.band_count(presence=True, samples=samples) == len(samples)].index)
        return len(self.__binary_matrix.loc[self.band_count(presence=True) == self.shape()[1]].index)

    def polymorphic_bands(self, samples=None):
        if samples:
            return self.__binary_matrix.loc[self.band_count(presence=True, samples=samples) != len(samples)][samples]
        return self.__binary_matrix.loc[self.band_count(presence=True) != self.shape()[1]]

    def polymorphic_band_count(self, samples=None):
        # Row wise Total band count == sample count
        if samples:
            return len(self.__binary_matrix.loc[self.band_count(presence=True, samples=samples) != len(samples)].index)
        return len(self.__binary_matrix.loc[self.band_count(presence=True) != self.shape()[1]].index)

    def percent_polymorphism(self, samples=None):
        return self.polymorphic_band_count(samples=samples) * 100/self.total_bands(samples=samples)

    def get_marker_frequency(self):
        return self.total_bands(unique=False) / (self.total_bands(unique=True) * self.__sample_count)

    def pic_dominant_1(self):
        """
        Calculate Polymorphic Information Content for a Dominant marker (Roldan-Ruiz et al., 2000)
        PIC = 2f(1-f)

        :return: A list of PIC for individual bands/allele or average PIC value
        """
        marker_frequency = self.get_marker_frequency()
        return 2 * marker_frequency * (1 - marker_frequency)

    def pic_dominant_2(self):
        """
        Calculate Polymorphic Information Content for a Dominant marker (De Riek et al., 2001)
        PIC = 1 - [f^2 + (1 - f)^2]

        :return: A list of PIC for individual bands/allele or average PIC value
        """
        marker_frequency = self.get_marker_frequency()
        return 1 - (np.square(marker_frequency) + np.square((1 - marker_frequency)))

    def pic_dominant_3(self, average=True):
        """
        Calculate Polymorphic Information Content for a Dominant marker
        PIC = 1 - (p^2 + q^2), where p = frequency of presence and q = frequency of absence in a particular allele

        :param average: (boolean) Output The average PIC

        :return: A list of PIC for individual bands/allele or average PIC value
        """

