from os import path, listdir

import pandas as pd

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase
from staramr.blast.pointfinder.PointfinderDatabaseInfo import PointfinderDatabaseInfo

"""
A Class for a PointFinder database for a particular organism.
"""


class PointfinderBlastDatabase(AbstractBlastDatabase):
    def __init__(self, database_dir, organism):
        """
        Creates a new PointfinderBlastDatabase.
        :param database_dir: The root directory for the PointFinder organisms.
        :param organism: The particular organism.
        """
        super().__init__(database_dir)
        self.organism = organism
        self.pointfinder_database_dir = path.join(self.database_dir, organism)

        if (not path.isdir(self.pointfinder_database_dir)):
            raise Exception(
                "Error, pointfinder organism [" + organism + "] is either incorrect or pointfinder database not installed properly")
        elif organism not in PointfinderBlastDatabase.get_organisms(database_dir):
            raise Exception("Pointfinder organism [" + organism + "] is not valid")

        self._pointfinder_info = PointfinderDatabaseInfo.from_file(
            path.join(self.pointfinder_database_dir, "resistens-overview.txt"))

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in listdir(self.pointfinder_database_dir) if
                (path.isfile(path.join(self.pointfinder_database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return path.join(self.pointfinder_database_dir, database_name + self.fasta_suffix)

    def is_validated(self):
        """
        Whether or not this particular PointFinder organism is part of the validated set.
        :return: True if this PointFinder organism/database is validated, False otherwise.
        """
        return self.organism in self.get_available_organisms()

    def get_resistance_codons(self, gene, codon_mutations):
        """
        Gets a list of resistance codons from the given gene and codon mutations.
        :param gene: The gene.
        :param codon_mutations: The codon mutations.
        :return: The resistance codons.
        """
        return self._pointfinder_info.get_resistance_codons(gene, codon_mutations)

    def get_phenotype(self, gene, codon_mutation):
        """
        Gets the phenotype for a given gene and codon mutation from PointFinder.
        :param gene: The gene.
        :param codon_mutation: The codon mutation.
        :return: A string describing the phenotype.
        """
        return self._pointfinder_info.get_phenotype(gene, codon_mutation)

    def get_resistance_nucleotides(self, gene, nucleotide_mutations):
        """
        Gets a list of resistance nucleotides from the given gene and nucleotide mutations.
        :param gene: The gene.
        :param nucleotide_mutations: The nucleotide mutations.
        :return: The resistance nucleotides.
        """
        return self._pointfinder_info.get_resistance_nucleotides(gene, nucleotide_mutations)

    def get_resistance_promoter(self, gene, nucleotide_mutations):
        """
        Gets a list of resistance nucleotides and codons located within the promoter, derived from the list
        of nucleotide mutations.
        :param gene: The name of the gene.
        :param nucleotide_mutations: The positions of the nucleotide mutations.
        :return: The resistance positions (both nucleotide and codon positions).
        """

        # Get the mutations in the nucleotide (non-gene) part:
        # Filter the list for negative coordinate positions.
        nucleotide_part = list(filter(lambda x: (x._nucleotide_position_amr_gene < 0), nucleotide_mutations))
        resistance_nucleotides = self._pointfinder_info.get_resistance_nucleotides(gene, nucleotide_part)

        # Get the mutations in the coding part:
        # Filter the list for non-negative coordinate positions.
        codon_part = list(filter(lambda x: (x._nucleotide_position_amr_gene >= 0), nucleotide_mutations))
        resistance_codons = self._pointfinder_info.get_resistance_codons(gene, codon_part)

        # Combine and return the results:
        return resistance_nucleotides + resistance_codons

    def get_organism(self):
        """
        Gets the particular organism of this database.
        :return: The organism.
        """
        return self.organism

    def get_name(self):
        return 'pointfinder'

    @classmethod
    def get_available_organisms(cls):
        """
        A Class Method to get a list of organisms that are currently supported by staramr.
        :return: The list of organisms currently supported by staramr.
        """
        return ['salmonella', 'campylobacter', 'enterococcus_faecalis', 'enterococcus_faecium',
                'escherichia_coli', 'helicobacter_pylori']

    @classmethod
    def get_organisms(cls, database_dir):
        """
        A Class Method to get the list of organisms from the PointFinder database root directory.
        :param database_dir: The PointFinder database root directory.
        :return: A list of organisms.
        """
        config = pd.read_csv(path.join(database_dir, 'config'), sep='\t', comment='#', header=None,
                             names=['db_prefix', 'name', 'description'])
        return config['db_prefix'].tolist()

    @classmethod
    def build_databases(cls, database_dir):
        """
        A Class Method to build a list of PointfinderBlastDatabase for all organisms in PointFinder.
        :param database_dir: The root PointFinder database directory.
        :return: A list of PointfinderBlastDatabase objects for all organisms in PointFinder.
        """
        return [PointfinderBlastDatabase(database_dir, organism) for organism in cls.get_organisms(database_dir)]
