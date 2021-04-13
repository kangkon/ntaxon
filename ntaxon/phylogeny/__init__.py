from Bio.Phylo.TreeConstruction import DistanceMatrix as Dm
from scipy.spatial.distance import squareform

class DistanceMatrix(Dm):
    def squareform(self):
        print(self.matrix)