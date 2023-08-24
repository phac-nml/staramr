from staramr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder

"""
Class used to parse out BLAST results for ResFinder, including phenotyhpes/resistances.
"""


class BlastResultsParserResfinderResistance(BlastResultsParserResfinder):
    COLUMNS = [x.strip() for x in '''
    Isolate ID
    Gene
    Predicted Phenotype
    CGE Predicted Phenotype
    %Identity
    %Overlap
    HSP Length/Total Length
    Contig
    Start
    End
    Accession
    Sequence
    CGE Notes
    '''.strip().split('\n')]

    def __init__(self, file_blast_map, arg_drug_table, cge_drug_table, blast_database, pid_threshold,
                 plength_threshold, report_all=False, output_dir=None, genes_to_exclude=[]):
        """
        Creates a new BlastResultsParserResfinderResistance.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param arg_drug_table: A table mapping the resistance gene to a specific drug resistance.
        :param cge_drug_table: A table mapping the CGE-predicted resistance gene to a specific drug resistance.
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
        self._cge_drug_table = cge_drug_table

    def _get_result_rows(self, hit, database_name):
        arg_drug = self._arg_drug_table.get_drug(database_name, hit.get_amr_gene_name_with_variant(),
                                                 hit.get_amr_gene_accession())

        cge_drug = self._cge_drug_table.get_drug(database_name, hit.get_amr_gene_name_with_variant(),
                                                 hit.get_amr_gene_accession())

        if arg_drug is None:
            arg_drug = 'unknown[' + hit.get_amr_gene_variant_accession() + ']'

        if cge_drug is None:
            cge_drug = 'unknown[' + hit.get_amr_gene_variant_accession() + ']'

        return [[hit.get_genome_id(),
                 hit.get_amr_gene_name(),
                 arg_drug,
                 cge_drug,
                 hit.get_pid(),
                 hit.get_plength(),
                 str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                 hit.get_genome_contig_id(),
                 hit.get_genome_contig_start(),
                 hit.get_genome_contig_end(),
                 hit.get_amr_gene_accession(),
                 hit.get_genome_contig_hsp_seq(),
                 self._cge_drug_table.get_cge_notes(hit.get_amr_gene_name_with_variant(),
                                                hit.get_amr_gene_accession())
                 ]]
