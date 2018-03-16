from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


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

    def _get_mutation_positions(self, start, database_frame, query_frame):
        return [NucleotideMutationPosition(i, self.hsp.sbjct, self.hsp.query, start, database_frame, query_frame) for i
                in self._get_match_positions()]
