from Bio.Seq import Seq
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


class Sequence(Seq):
    def restriction_search(self, enzyme, linear=True):
        #amb_seq = Seq(self, IUPACAmbiguousDNA())
        #return enzyme.search(amb_seq)
        return enzyme.search(self, linear)

    def haplotypes(self, enzyme, linear=True, min_size=50):
        fragments = enzyme.catalyse(self, linear)
        lengths = [len(f) for f in fragments]
        haplotypes = list()
        for i, sz in enumerate(lengths):
            if sz > min_size:
                haplotypes.append([sz, fragments[i]])
        return haplotypes
