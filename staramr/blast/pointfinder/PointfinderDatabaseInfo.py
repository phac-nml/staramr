import logging
from os import path

import pandas as pd
import Bio.Seq
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition
from staramr.blast.results.pointfinder.codon.CodonInsertionPosition import CodonInsertionPosition

from staramr.exceptions.GenotypePhenotypeMatchException import GenotypePhenotypeMatchException

"""
A Class storing information about the specific PointFinder database.
"""

logger = logging.getLogger('PointfinderDatabaseInfo')


class PointfinderDatabaseInfo:

    def __init__(self, database_info_dataframe, file=None):
        """
        Creates a new PointfinderDatabaseInfo.
        :param database_info_dataframe: A pd.DataFrame containing the information in PointFinder.
        :param file: The file where the pointfinder database info originates from.
        """
        self._pointfinder_info = database_info_dataframe
        self._file = file

        self._resistance_table_hacks(self._pointfinder_info)

    @classmethod
    def from_file(cls, file):
        """
        Builds a new PointfinderDatabaseInfo from the passed file containing PointFinder information on drug resistance
        mutations.
        :param file: The file containing drug resistance mutations.
        :return: A new PointfinderDatabaseInfo.
        """

        with open(file) as f:
            line = f.readline()
        
        line = line.lstrip("#")
        column_names = line.split()

        pointfinder_info = pd.read_csv(file, sep='\t', index_col=False, comment='#', header=None, names=column_names)

        return cls(pointfinder_info, file)

    @classmethod
    def from_pandas_table(cls, database_info_dataframe):
        """
        Builds a new PointfinderDatabaseInfo from the passed pd.DataFrame.
        :param database_info_dataframe: A pd.DataFrame containing the information in PointFinder.
        :return: A new PointfinderDatabaseInfo.
        """
        return cls(database_info_dataframe)
    
    @staticmethod
    def to_codons(regex_match):
        # Sometimes, the regex will match a string with a comma and return multiple matches.
        # Ex: TTCATGGAC,TTC
        # We need to make sure we handle these commas and multiple matches.
        entries = regex_match.string.split(",")

        for i in range(0, len(entries)):
            entries[i] = Bio.Seq.translate(entries[i], table='Standard').upper()

        return ",".join(entries)

    def _resistance_table_hacks(self, table):
        """
        A function implementing some hacks to try and fix mismatched strings in the pointfinder databases.
        These should be removed when the underlying database is corrected.
        :param table: The pointfinder resistance table to fix.
        :return: None, but modifies the passed table in place.
        """
        if self._file and 'salmonella' in str(self._file) and path.exists(
                path.join(path.dirname(self._file), '16S_rrsD.fsa')):
            logger.debug("Replacing [16S] with [16S_rrsD] for pointfinder organism [salmonella]")
            table[['Gene_ID']] = table[['Gene_ID']].replace('16S', '16S_rrsD')

        # We need to fix table entries where codon deletions are listed,
        # but no reference nucleotides are listed in the corresponding entry.
        # This will mean that this specific entry (for some reason)
        # is using 0-based nucleotide coordinates and nucleotides for the Res_codon
        # instead of 1-based codon coordinates and codons for the Res_codon.
        # However, we need make an exception for entries where "promoter", "16S", or "23S"
        # are found in the 'Gene_ID'.

        table.loc[(table['Ref_nuc'] == '-') & (table['Ref_codon'] == 'del')
        & (table['Gene_ID'].str.contains('promoter') == False)
        & (table['Gene_ID'].str.contains('16S') == False)
        & (table['Gene_ID'].str.contains('23S') == False), "Codon_pos"] += 3
        # We increment by 3 because the table is using base-0 coordinates
        # and we have base-1 coordindates for everything codon-related.
        # However, what's listed is nucleotide coordinates for the codon
        # mutation, so we actually need to bump it by 3 nucleotides.

        # Furthermore, we need to convert 'Res_codon' entries related to the above problem
        # (codon deletions) from nucleotides into codons. We choose to convert nucleotides into
        # codons instead of the other way around so that can more accurately compare observed mutations
        # to the table entries.
        table.loc[(table['Ref_nuc'] == '-') & (table['Ref_codon'] == 'del')
        & (table['Gene_ID'].str.contains('promoter') == False)
        & (table['Gene_ID'].str.contains('16S') == False)
        & (table['Gene_ID'].str.contains('23S') == False), "Res_codon"] = table.loc[(table['Ref_nuc'] == '-') & (table['Ref_codon'] == 'del')
        & (table['Gene_ID'].str.contains('promoter') == False)
        & (table['Gene_ID'].str.contains('16S') == False)
        & (table['Gene_ID'].str.contains('23S') == False), "Res_codon"].str.replace('[A-Z,]+', self.to_codons, regex=True)

    def _get_resistance_codon_match(self, gene, codon_mutation):
        table = self._pointfinder_info

        # We need to handle codon deletions as a special case:
        if(type(codon_mutation) is CodonMutationPosition
        and codon_mutation.get_input_genome_amino_acid() == 'del'):
            matches = table[(table['Gene_ID'] == gene)
                    # Codon deletions are listed using nucleotide coordinates
                    & (table['Codon_pos'] == codon_mutation.get_mutation_position() * 3)
                        # Such coords are denoted in nucleotide coordinates in the reference table for some reason,
                        # so we need to convert to nucleotide coordinates before making the comparison.
                    & (table['Ref_codon'] == codon_mutation.get_database_amr_gene_mutation())
                    & (table['Res_codon'].str.contains(codon_mutation.get_input_genome_mutation(), regex=False))]
        
        # We need to handle codon insertions as a special case:
        # Pointfinder mis-reports the position of codon insertions. For example:
        # ref:     ACG --- ACG
        # query:   ACG GGG ACG
        # ref_pos: 1       2
        # Pointfinder is incorrectly reporting the insertion as 2_3insG instead of the correct 1_2insG,
        # which is to say, the insertion happens between reference codon coordinates 1 and 2.
        # We need to shift by 1 in our interpretation.
        elif type(codon_mutation) is CodonInsertionPosition:
            matches = table[(table['Gene_ID'] == gene)
                    # Codon inerstions need to be shifted by 1:
                    & (table['Codon_pos'] == codon_mutation.get_mutation_position() + 1)
                    & (table['Ref_codon'] == codon_mutation.get_database_amr_gene_mutation())
                    & (table['Res_codon'].str.contains(codon_mutation.get_input_genome_mutation(), regex=False))]

        # Normal cases:
        else:
            matches = table[(table['Gene_ID'] == gene)
                            & (table['Codon_pos'] == codon_mutation.get_mutation_position())
                            & (table['Ref_codon'] == codon_mutation.get_database_amr_gene_mutation())
                            & (table['Res_codon'].str.contains(codon_mutation.get_input_genome_mutation(), regex=False))]

        if len(matches.index) > 1:
            # If more then one match, try to match Res_codon exactly to subselect
            matches_subset = matches[matches['Res_codon'] == codon_mutation.get_input_genome_mutation()]

            if len(matches_subset.index) >= 1:
                matches = matches_subset

        return matches

    def _get_resistance_nucleotide_match(self, gene, nucleotide_mutations):
        return self._get_resistance_codon_match(gene, nucleotide_mutations)

    def get_phenotype(self, gene, codon_mutation):
        """
        Gets the phenotype for a given gene and codon mutation from PointFinder.
        :param gene: The gene.
        :param codon_mutation: The codon mutation.
        :return: A string describing the phenotype.
        """
        match = self._get_resistance_codon_match(gene, codon_mutation)

        if len(match.index) > 0:
            return match['Resistance'].iloc[0]
        else:
            raise GenotypePhenotypeMatchException("Error, no match for gene=" + str(gene) + ", codon_mutation=" + str(codon_mutation))

    def get_resistance_codons(self, gene, codon_mutations):
        """
        Gets a list of resistance codons from the given gene and codon mutations.
        :param gene: The gene.
        :param codon_mutations: The codon mutations.
        :return: The resistance codons.
        """
        resistance_mutations = []

        for codon_mutation in codon_mutations:
            match = self._get_resistance_codon_match(gene, codon_mutation)
            if len(match.index) > 0:
                resistance_mutations.append(codon_mutation)

        return resistance_mutations

    def get_resistance_nucleotides(self, gene, nucleotide_mutations):
        """
        Gets a list of resistance nucleotides from the given gene and nucleotide mutations.
        :param gene: The gene.
        :param nucleotide_mutations: The nucleotide mutations.
        :return: The resistance nucleotides.
        """
        resistance_mutations = []

        for nucleotide_mutation in nucleotide_mutations:
            match = self._get_resistance_nucleotide_match(gene, nucleotide_mutation)
            if len(match.index) > 0:
                resistance_mutations.append(nucleotide_mutation)

        return resistance_mutations
