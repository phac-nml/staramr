import pandas

from staramr.blast.results.BlastResultsParser import BlastResultsParser
from staramr.blast.results.resfinder.ResfinderHitHSP import ResfinderHitHSP

"""
Class used to parse out BLAST results for ResFinder.
"""


class BlastResultsParserResfinder(BlastResultsParser):

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold):
        """
        Creates a new BlastResultsParserResfinder.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold)

    def _create_hit(self, file, database_name, blast_record, alignment, hsp):
        return ResfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results, columns=('FILE', 'GENE', '%IDENTITY', '%OVERLAP',
                                                'DB_SEQ_LENGTH/QUERY_HSP', 'CONTIG', 'START', 'END', 'ACCESSION'))
        return df.set_index('FILE')

    def _append_results_to(self, hit, database_name, results):

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
