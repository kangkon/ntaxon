import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ntaxon.nucleotide.sequence import Sequence


class RestrictionDigestionAlignment:
    aligned_fragments = None

    def __init__(self, digestion_table: pd.DataFrame):
        self.aligned_fragments = None

    def plot(self):
        print('plot')


class RestrictionDigestion:
    accessions = None
    digestion_table = None
    digestion_profile = None

    def __init__(self, enzyme, accessions: pd.DataFrame, label_col: str = 'isolate', sequence_col: str = 'sequence'):
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
        digestion['haplotype'] = digestion['size'].astype(str) + '_' + digestion['haplotype'].astype(str)

        self.digestion_table = digestion

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

        self.digestion_profile = digestion_bin

    def plot_histogram(self, figsize=(14, 8)):
        plt.figure(figsize=figsize)
        ax = sns.histplot(data=self.digestion_table, x="size")
        plt.show()

    def plot_electrophoretic_diagram(self, figsize=(16, 10)):
        plt.figure(figsize=figsize)
        ax = sns.heatmap(self.digestion_profile.drop('size', axis=1), cmap="YlGnBu", cbar=False)
        plt.show()

    def get_fragment_alignment(self):
        print('Alignment')

