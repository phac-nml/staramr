import logging

from amr.blast.results.AMRHitHSP import AMRHitHSP
from amr.blast.results.pointfinder.NucleotideMutationPosition import NucleotideMutationPosition

logger = logging.getLogger('PointfinderHitHSP')


class PointfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        super().__init__(file, blast_record, hit, hsp)

    def get_gene(self):
        return self.hit.hit_id

    def _get_hsp_frame(self, index):
        frame = self.hsp.frame[index]
        if frame not in [1, -1]:
            raise Exception("frame=" + str(frame) + ", is unexpected")
        else:
            return frame

    def get_database_frame(self):
        return self._get_hsp_frame(1)

    def get_query_frame(self):
        return self._get_hsp_frame(0)

    def _get_match_positions(self):
        return [i for i, c in enumerate(self.hsp.match) if c == ' ']

    def _get_nucleotide_mutation_positions(self, start, database_frame, query_frame):
        return [NucleotideMutationPosition(i, self.hsp.sbjct, self.hsp.query, start, database_frame, query_frame) for i
                in self._get_match_positions()]

    def get_nucleotide_mutations(self):
        start = self.hsp.sbjct_start
        database_frame = self.get_database_frame()
        query_frame = self.get_query_frame()
        return self._get_nucleotide_mutation_positions(start, database_frame, query_frame)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
