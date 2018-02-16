from amr.blast.results.AMRHitHSP import AMRHitHSP

import logging
import pprint
import math


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

    def _get_nucleotide_mutation_positions(self, start, frame):
        mutations = []
        if frame == 1:
            mutations = [start + i for i, c in enumerate(self.hsp.match) if c == ' ']
        else:
            mutations = [start - i for i, c in enumerate(self.hsp.match) if c == ' ']
        return mutations

    def get_database_nucleotide_mutation_positions(self):
        start = self.hsp.sbjct_start
        frame = self.get_database_frame()
        return self._get_nucleotide_mutation_positions(start,frame)

    def get_query_nucleotide_mutation_positions(self):
        start = self.hsp.query_start
        frame = self.get_query_frame()
        return self._get_nucleotide_mutation_positions(start,frame)

    def get_codon_mutation_positions_at(self, nucleotide_mutation_positions):
        codon_mutation_positions = [math.ceil(x/3) for x in nucleotide_mutation_positions]
        logger.info("codon_mutations: "+str(codon_mutation_positions))
        return codon_mutation_positions

    def _get_nucleotide_codons_at(self, nucleotide_string, nucleotide_mutation_positions):
        nucleotide_string = self.hsp.sbjct
        codon_position_starts_index_0 = [(x - 1) - ((x - 1) % 3) for x in nucleotide_mutation_positions]
        return [nucleotide_string[x:(x+3)].upper() for x in codon_position_starts_index_0]

    def get_database_nucleotide_codons_at(self, nucleotide_mutation_positions):
        return self._get_nucleotide_codons_at(self.hsp.sbjct, nucleotide_mutation_positions)

    def get_query_nucleotide_codons_at(self, nucleotide_mutation_positions):
        return self._get_nucleotide_codons_at(self.hsp.query, nucleotide_mutation_positions)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
