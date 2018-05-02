import pandas

from staramr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder

"""
Class used to parse out BLAST results for ResFinder, including phenotyhpes/resistances.
"""


class BlastResultsParserResfinderResistance(BlastResultsParserResfinder):

    def __init__(self, file_blast_map, arg_drug_table, blast_database, pid_threshold, plength_threshold,
                 report_all=False, output_dir=None):
        """
        Creates a new BlastResultsParserResfinderResistance.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param arg_drug_table: A table mapping the resistance gene to a specific drug resistance.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all,
                         output_dir=output_dir)
        self._arg_drug_table = arg_drug_table

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results, columns=('Isolate ID', 'Gene', 'Predicted Phenotype', '%Identity', '%Overlap',
                                                'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession'))
        return df.set_index('Isolate ID')

    def _append_results_to(self, hit, database_name, results, seq_records):
        self._append_seqrecords_to(hit, seq_records)
        drug = self._arg_drug_table.get_drug(database_name, hit.get_gene_with_variant(), hit.get_accession())

        if drug is None:
            drug = 'unknown['+hit.get_gene_variant_accession()+']'

        results.append([hit.get_isolate_id(),
                        hit.get_gene(),
                        drug,
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length()),
                        hit.get_contig(),
                        hit.get_contig_start(),
                        hit.get_contig_end(),
                        hit.get_accession()
                        ])
