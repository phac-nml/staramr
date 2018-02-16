import os
import pandas

from amr.blast.results.BlastResultsParser import BlastResultsParser
from amr.blast.results.resfinder.ResfinderHitHSP import ResfinderHitHSP

class BlastResultsParserResfinder(BlastResultsParser):

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold):
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold)

    def _create_hit(self, file, blast_record, alignment, hsp):
        return ResfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        return pandas.DataFrame(results, columns=('FILE', 'GENE', 'RESFINDER_PHENOTYPE', '%IDENTITY', '%OVERLAP',
                                                  'DB_SEQ_LENGTH/QUERY_HSP', 'CONTIG', 'START', 'END', 'ACCESSION'))

    def _append_results_to(self, hit, results):
        phenotype=self._blast_database.get_phenotype(hit.get_gene())

        results.append([hit.get_file(),
                        hit.get_gene(),
                        phenotype,
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length()),
                        hit.get_contig(),
                        hit.get_contig_start(),
                        hit.get_contig_end(),
                        hit.get_accession()
                       ])