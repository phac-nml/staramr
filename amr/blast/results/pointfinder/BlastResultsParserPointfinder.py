import os
import pandas

from amr.blast.results.BlastResultsParser import BlastResultsParser
from amr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP

class BlastResultsParserPointfinder(BlastResultsParser):

    def __init__(self, file_blast_map, pid_threshold, plength_threshold):
        super().__init__(file_blast_map, pid_threshold, plength_threshold)

    def _create_hit(self, file, alignment, hsp):
        return PointfinderHitHSP(file, alignment, hsp)

    def _create_data_frame(self, results):
        return pandas.DataFrame(results, columns=('File', 'Resistance gene', '% Identity', '% Overlap', 'HSP/Alignment'))

    def _append_results_to(self, hit, results):
        results.append([hit.get_file(),
                        hit.get_hit_id(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                       ])