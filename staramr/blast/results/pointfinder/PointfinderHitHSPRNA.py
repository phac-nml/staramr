from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP

class PointfinderHitHSPRNA(PointfinderHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        """
        Creates a new PointfinderHitHSPRNA.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param hit: The particular Bio.Blast.Record.Alignment.
        :param hsp: The particular Bio.Blast.Record.HSP.
        """
        super().__init__(file, blast_record, hit, hsp)