import pandas


class PointfinderDatabaseInfo:

    def __init__(self, table):
        self._pointfinder_info = table

    @classmethod
    def from_file(cls, file):
        pointfinder_info = pandas.read_csv(file, sep="\t", index_col=False)
        return cls(pointfinder_info)

    @classmethod
    def from_pandas_table(cls, table):
        return cls(table)

    def _get_resistance_codon_match(self, gene, codon_mutation):
        table = self._pointfinder_info

        matches = table[(table['#Gene_ID'] == gene)
                        & (table['Codon_pos'] == codon_mutation.get_codon_start())
                        & (table['Ref_codon'] == codon_mutation.get_database_amino_acid())
                        & (table['Res_codon'].str.contains(codon_mutation.get_query_amino_acid().upper(), regex=False))]

        if len(matches.index) > 1:
            raise Exception("Error, multiple matches for gene=" + str(gene) + ", codon_mutation=" + str(codon_mutation))
        else:
            return matches

    def get_phenotype(self, gene, codon_mutation):
        match = self._get_resistance_codon_match(gene, codon_mutation)

        if len(match.index) > 0:
            return match['Resistance'].iloc[0]
        else:
            raise Exception("Error, no match for gene=" + str(gene) + ", codon_mutation=" + str(codon_mutation))

    def get_resistance_codons(self, gene, codon_mutations):
        resistance_mutations = []

        table = self._pointfinder_info
        for codon_mutation in codon_mutations:
            match = self._get_resistance_codon_match(gene, codon_mutation)
            if len(match.index) > 0:
                resistance_mutations.append(codon_mutation)

        return resistance_mutations
