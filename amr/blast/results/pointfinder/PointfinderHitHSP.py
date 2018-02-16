from amr.blast.results.AMRHitHSP import AMRHitHSP

import logging
import pprint
import math


logger = logging.getLogger('PointfinderHitHSP')

class PointfinderHitHSP(AMRHitHSP):

    pp = pprint.PrettyPrinter(indent=4)

    def __init__(self, file, blast_record, hit, hsp):
        super().__init__(file, blast_record, hit, hsp)

    def get_gene(self):
        return self.hit.hit_id

    def get_database_frame(self):
        frame = self.hsp.frame[1]
        if frame not in [1, -1]:
            raise Exception("frame=" + str(frame) + ", is unexpected")
        else:
            return frame

    def get_nucleotide_mutation_positions(self):
        start = self.hsp.sbjct_start
        mutations = []
        if self.get_database_frame() == 1:
            mutations = [start + i for i, c in enumerate(self.hsp.match) if c == ' ']
        else:
            mutations = [start - i for i, c in enumerate(self.hsp.match) if c == ' ']
        logger.info("")
        logger.info("gene: "+self.hit.title)
        logger.info("sbjct_start: " + str(self.hsp.sbjct_start))
        logger.info("sbjct_start: " + str(self.hsp.sbjct_end))
        logger.info("nuc_mutations: "+str(mutations))
        #logger.info("hsp: "+self.pp.pprint(vars(self.hsp)))
        return mutations

    def get_codon_mutation_positions_at(self, nucleotide_mutation_positions):
        codon_mutation_positions = [math.ceil(x/3) for x in nucleotide_mutation_positions]
        logger.info("codon_mutations: "+str(codon_mutation_positions))
        return codon_mutation_positions

    def get_nucleotide_codons_at(self, nucleotide_mutation_positions):
        codon_position_remainders_index_0 = [(x-1)%3 for x in nucleotide_mutation_positions]
        return "X"

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
