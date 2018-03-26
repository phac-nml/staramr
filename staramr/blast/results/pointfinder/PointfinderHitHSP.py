import logging

from staramr.blast.results.AMRHitHSP import AMRHitHSP
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition

logger = logging.getLogger('PointfinderHitHSP')

"""
Class storing a PointFinder BLAST hit/hsp.
"""


class PointfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        """
        Creates a new PointfinderHitHSP.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param hit: The particular Bio.Blast.Record.Alignment.
        :param hsp: The particular Bio.Blast.Record.HSP.
        """
        super().__init__(file, blast_record, hit, hsp)

    def get_gene(self):
        """
        Gets the particular gene name for the PointFinder hit.
        :return: The gene name.
        """
        return self.hit.hit_id

    def _get_match_positions(self):
        return [i for i, c in enumerate(self.hsp.match) if c == ' ']

    def _get_mutation_positions(self, start, database_frame, query_frame):
        return [CodonMutationPosition(i, self.hsp.sbjct, self.hsp.query, start, database_frame, query_frame) for i
                in self._get_match_positions()]

    def get_mutations(self):
        """
        Gets a list of NucleotideMutationPosition for the individual mutations.
        :return: A list of NucleotideMutationPosition.
        """
        start = self.hsp.sbjct_start
        database_frame = self.get_database_frame()
        query_frame = self.get_query_frame()
        return self._get_mutation_positions(start, database_frame, query_frame)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
