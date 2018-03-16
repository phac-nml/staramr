import pandas

"""
A Class storing information about the specific PointFinder database.
"""


class PointfinderDatabaseInfo:

    def __init__(self, database_info_dataframe):
        """
        Creates a new PointfinderDatabaseInfo.
        :param database_info_dataframe: A pandas.DataFrame containing the information in PointFinder.
        """
        self._pointfinder_info = database_info_dataframe

    @classmethod
    def from_file(cls, file):
        """
        Builds a new PointfinderDatabaseInfo from the passed file containing PointFinder information on drug resistance
        mutations.
        :param file: The file containing drug resistance mutations.
        :return: A new PointfinderDatabaseInfo.
        """
        pointfinder_info = pandas.read_csv(file, sep="\t", index_col=False)
        return cls(pointfinder_info)

    @classmethod
    def from_pandas_table(cls, database_info_dataframe):
        """
        Builds a new PointfinderDatabaseInfo from the passed pandas.DataFrame.
        :param database_info_dataframe: A pandas.DataFrame containing the information in PointFinder.
        :return: A new PointfinderDatabaseInfo.
        """
        return cls(database_info_dataframe)

    def _get_resistance_codon_match(self, gene, codon_mutation):
        table = self._pointfinder_info

        matches = table[(table['#Gene_ID'] == gene)
                        & (table['Codon_pos'] == codon_mutation.get_mutation_position())
                        & (table['Ref_codon'] == codon_mutation.get_database_mutation())
                        & (table['Res_codon'].str.contains(codon_mutation.get_query_mutation(), regex=False))]

        if len(matches.index) > 1:
            raise Exception("Error, multiple matches for gene=" + str(gene) + ", codon_mutation=" + str(codon_mutation))
        else:
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
            raise Exception("Error, no match for gene=" + str(gene) + ", codon_mutation=" + str(codon_mutation))

    def get_resistance_codons(self, gene, codon_mutations):
        """
        Gets a list of resistance codons from the given gene and codon mutations.
        :param gene: The gene.
        :param codon_mutations: The codon mutations.
        :return: The resistance codons.
        """
        resistance_mutations = []

        table = self._pointfinder_info
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

        table = self._pointfinder_info
        for nucleotide_mutation in nucleotide_mutations:
            match = self._get_resistance_nucleotide_match(gene, nucleotide_mutation)
            if len(match.index) > 0:
                resistance_mutations.append(nucleotide_mutation)

        return resistance_mutations
