import os

from amr.blast.results.BlastResultsParser import BlastResultsParser
from amr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP

class BlastResultsParserPointfinder(BlastResultsParser):

    def __init__(self, file_blast_map, pid_threshold, plength_threshold):
        super().__init__(file_blast_map, pid_threshold, plength_threshold)

    def _create_hit(self, file, alignment, hsp):
        return PointfinderHitHSP(file, alignment, hsp)

    def _print_hit(self, hit, file_handle):
        file_handle.write("%s\t%s\t%0.2f\t%0.2f\t%d/%d\n" % (
            hit.get_file(), hit.get_hit_id(), hit.get_pid(), hit.get_plength(), hit.get_hsp_alignment_length(),
            hit.get_alignment_length()))