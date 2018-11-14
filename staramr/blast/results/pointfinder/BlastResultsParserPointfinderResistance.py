import logging

from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder

"""
Class used to parse out BLAST results for PointFinder, including phenotyhpes/resistances.
"""

logger = logging.getLogger('BlastResultsParserPointfinderResistance')


class BlastResultsParserPointfinderResistance(BlastResultsParserPointfinder):
    COLUMNS = [x.strip() for x in '''
    Isolate ID
    Gene
    Predicted Phenotype
    Type
    Position
    Mutation
    %Identity
    %Overlap
    HSP Length/Total Length
    Contig
    Start
    End
    '''.strip().split('\n')]

    def __init__(self, file_blast_map, arg_drug_table, blast_database, pid_threshold, plength_threshold,
                 report_all=False, output_dir=None, genes_to_exclude=[]):
        """
        Creates a new BlastResultsParserPointfinderResistance.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param arg_drug_table: A table mapping the resistance gene to a specific drug resistance.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all,
                         output_dir=output_dir, genes_to_exclude=genes_to_exclude)
        self._arg_drug_table = arg_drug_table

    def _get_result(self, hit, db_mutation):
        drug = self._arg_drug_table.get_drug(self._blast_database.get_organism(), hit.get_amr_gene_id(),
                                             db_mutation.get_mutation_position())
        gene_name = hit.get_amr_gene_id() + " (" + db_mutation.get_mutation_string_short() + ")"

        if drug is None:
            drug = 'unknown[' + gene_name + ']'

        return [hit.get_genome_id(),
                gene_name,
                drug,
                db_mutation.get_type(),
                db_mutation.get_mutation_position(),
                db_mutation.get_mutation_string(),
                hit.get_pid(),
                hit.get_plength(),
                str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                hit.get_genome_contig_id(),
                hit.get_genome_contig_start(),
                hit.get_genome_contig_end()
                ]
