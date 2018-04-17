from os import path

import pandas

from staramr.blast.results.BlastResultsParser import BlastResultsParser
from staramr.blast.results.resfinder.ResfinderHitHSP import ResfinderHitHSP

"""
Class used to parse out BLAST results for ResFinder.
"""


class BlastResultsParserResfinder(BlastResultsParser):

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold, report_all=False,
                 output_dir=None):
        """
        Creates a new BlastResultsParserResfinder.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all,
                         output_dir=output_dir)

    def _create_hit(self, file, database_name, blast_record, alignment, hsp):
        return ResfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results, columns=('Isolate ID', 'Gene', '%Identity', '%Overlap',
                                                'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession'))
        return df.set_index('Isolate ID')

    def _append_results_to(self, hit, database_name, results, seq_records):
        self._append_seqrecords_to(hit, seq_records)
        results.append([hit.get_isolate_id(),
                        hit.get_gene(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length()),
                        hit.get_contig(),
                        hit.get_contig_start(),
                        hit.get_contig_end(),
                        hit.get_accession()
                        ])

    def _get_out_file_name(self, in_file):
        if self._output_dir:
            return path.join(self._output_dir, 'resfinder_' + path.basename(in_file))
        else:
            raise Exception("output_dir is None")
