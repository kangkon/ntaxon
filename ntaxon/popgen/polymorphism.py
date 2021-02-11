from ntaxon.fingerprinting.matrix import BinaryMatrix
class Polymorphism:
    __locus_list = []
    def __init__(self, locus):
        # TODO: Validate locus list
        self.__locus_list = locus
