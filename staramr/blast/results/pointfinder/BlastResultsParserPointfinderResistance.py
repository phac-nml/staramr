import logging
import pandas

from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder
from staramr.blast.results.pointfinder.codon.CodonInsertionPosition import CodonInsertionPosition

"""
Class used to parse out BLAST results for PointFinder, including phenotyhpes/resistances.
"""

logger = logging.getLogger('BlastResultsParserPointfinderResistance')


class BlastResultsParserPointfinderResistance(BlastResultsParserPointfinder):
    COLUMNS = [x.strip() for x in '''
    Isolate ID
    Gene
    Predicted Phenotype
    CGE Predicted Phenotype
    Type
    Position
    Mutation
    %Identity
    %Overlap
    HSP Length/Total Length
    Contig
    Start
    End
    Pointfinder Position
    CGE Notes
    CGE Required Mutation
    CGE Mechanism
    CGE PMID
    '''.strip().split('\n')]

    def __init__(self, file_blast_map, arg_drug_table, blast_database, pid_threshold, plength_threshold,
                 report_all=False, output_dir=None, genes_to_exclude=[], complex_mutations=None):
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
        :param complex_mutations: An object mapping a set of multiple point mutations to a single phenotype.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all,
                         output_dir=output_dir, genes_to_exclude=genes_to_exclude)
        self._arg_drug_table = arg_drug_table
        self._complex_mutations = complex_mutations

    def _get_result(self, hit, db_mutation):

        # We need to correct for Pointfinder codon insertions being off by 1.
        if type(db_mutation) is CodonInsertionPosition:
            mutation_position = db_mutation.get_mutation_position() + 1
        else:
            mutation_position = db_mutation.get_mutation_position()

        arg_drug = self._arg_drug_table.get_drug(self._blast_database.get_organism(), hit.get_amr_gene_name(),
                                             mutation_position)
        
        cge_drug = self._blast_database.get_cge_phenotype(hit.get_amr_gene_name(), db_mutation)

        gene_name = hit.get_amr_gene_name() + " (" + db_mutation.get_mutation_string_short() + ")"

        if arg_drug is None:
            arg_drug = 'unknown[' + gene_name + ']'

        if cge_drug is None:
            cge_drug = 'unknown[' + gene_name + ']'

        result = [hit.get_genome_id(),
                gene_name,
                arg_drug,
                cge_drug,
                db_mutation.get_type(),
                db_mutation.get_mutation_position(),
                db_mutation.get_mutation_string(),
                hit.get_pid(),
                hit.get_plength(),
                str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                hit.get_genome_contig_id(),
                hit.get_genome_contig_start(),
                hit.get_genome_contig_end(),
                db_mutation.get_pointfinder_mutation_string(),
                self._blast_database.get_cge_notes(hit.get_amr_gene_name(), db_mutation),
                self._blast_database.get_cge_required_mutation(hit.get_amr_gene_name(), db_mutation),
                self._blast_database.get_cge_mechanism(hit.get_amr_gene_name(), db_mutation),
                self._blast_database.get_cge_pmid(hit.get_amr_gene_name(), db_mutation),
                ]

        return result
    
    def _get_result_rows(self, hit, database_name):
        results = super()._get_result_rows(hit, database_name)

        # We ought to handle complex mutations after results are processed, in case we need to account
        # for any Pointfinder position corrections (i.e. indels).
        if self._complex_mutations:
            results_table = pandas.DataFrame(columns=self.COLUMNS, data=results)
            matches = self._complex_mutations.get_matches(results_table, hit)
            
            if results is not None:
                results.extend(matches)
            else:
                results = matches

        return results

