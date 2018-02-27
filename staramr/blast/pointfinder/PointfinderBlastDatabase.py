from os import path, listdir

import pandas

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase
from staramr.blast.pointfinder.PointfinderDatabaseInfo import PointfinderDatabaseInfo


class PointfinderBlastDatabase(AbstractBlastDatabase):
    def __init__(self, database_dir, organism):
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

    def get_resistance_codons(self, gene, codon_mutations):
        return self._pointfinder_info.get_resistance_codons(gene, codon_mutations)

    def get_phenotype(self, gene, codon_mutation):
        return self._pointfinder_info.get_phenotype(gene, codon_mutation)

    def get_organism(self):
        return self.organism

    @classmethod
    def get_available_organisms(cls):
        return ['salmonella']

    @classmethod
    def get_organisms(cls, database_dir):
        config = pandas.read_csv(path.join(database_dir, 'config'), sep='\t', comment='#', header=None,
                                 names=['db_prefix', 'name', 'description'])
        return config['db_prefix'].tolist()

    @classmethod
    def build_databases(cls, database_dir):
        return [PointfinderBlastDatabase(database_dir, organism) for organism in cls.get_organisms(database_dir)]
