from staramr.blast.results.plasmidfinder.BlastResultsParserPlasmidfinder import BlastResultsParserPlasmidfinder
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from typing import Dict

"""
Class used to parse out BLAST results for PlasmidFinder, including phenotyhpes/resistances.
"""


class BlastResultsParserPlasmidfinderResistance(BlastResultsParserPlasmidfinder):
    COLUMNS = [x.strip() for x in '''
    Isolate ID
    Gene
    %Identity
    %Overlap
    HSP Length/Total Length
    Contig
    Start
    End
    Accession
    '''.strip().split('\n')]

    def __init__(self, file_blast_map: Dict[str, BlastResultsParserPlasmidfinder], blast_database: PlasmidfinderBlastDatabase, pid_threshold: int, plength_threshold: int,
                 report_all=False, output_dir=None, genes_to_exclude=[]) -> None:
        """
        Creates a new BlastResultsParserPlasmidfinderResistance.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        :param report_all: Whether or not to report all blast hits.
        :param output_dir: The directory where output files are being written.
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all,
                         output_dir=output_dir, genes_to_exclude=genes_to_exclude)

    def _get_result_rows(self, hit: list, database_name: str) -> list:

        return [[hit.get_genome_id(),
                 hit.get_amr_gene_name(),
                 hit.get_pid(),
                 hit.get_plength(),
                 str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                 hit.get_genome_contig_id(),
                 hit.get_genome_contig_start(),
                 hit.get_genome_contig_end(),
                 hit.get_amr_gene_accession()
                 ]]
