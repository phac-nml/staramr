import os
import pandas

from amr.blast.results.BlastResultsParser import BlastResultsParser
from amr.blast.results.resfinder.ResfinderHitHSP import ResfinderHitHSP

class BlastResultsParserResfinder(BlastResultsParser):

    def __init__(self, file_blast_map, pid_threshold, plength_threshold):
        super().__init__(file_blast_map, pid_threshold, plength_threshold)

    def _create_hit(self, file, blast_record, alignment, hsp):
        return ResfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        return pandas.DataFrame(results, columns=('FILE', 'GENE', '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP',
                                                  'CONTIG', 'START', 'END', 'ACCESSION'))

    def _append_results_to(self, hit, results):
        results.append([hit.get_file(),
                        hit.get_gene(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length()),
                        hit.get_contig(),
                        hit.get_contig_start(),
                        hit.get_contig_end(),
                        hit.get_accession()
                       ])