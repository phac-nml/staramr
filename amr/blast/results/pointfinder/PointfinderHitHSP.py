from amr.blast.results.AMRHitHSP import AMRHitHSP


class PointfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        super().__init__(file, blast_record, hit, hsp)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
